#include "online_evaluator.h"

#include <array>

#include "quadsquad/helpers.h"

namespace quadsquad {
OnlineEvaluator::OnlineEvaluator(int id,
                                 std::shared_ptr<io::NetIOMP<4>> network,
                                 PreprocCircuit<Ring> preproc,
                                 utils::LevelOrderedCircuit circ,
                                 int security_param, int threads, int seed)
    : id_(id),
      security_param_(security_param),
      rgen_(id, seed),
      network_(std::move(network)),
      preproc_(std::move(preproc)),
      circ_(std::move(circ)),
      wires_(circ.num_gates),
      jump_(id),
      msb_circ_(
          utils::Circuit<BoolRing>::generatePPAMSB().orderGatesByLevel()) {
  tpool_ = std::make_shared<ThreadPool>(threads);
}

OnlineEvaluator::OnlineEvaluator(int id,
                                 std::shared_ptr<io::NetIOMP<4>> network,
                                 PreprocCircuit<Ring> preproc,
                                 utils::LevelOrderedCircuit circ,
                                 int security_param,
                                 std::shared_ptr<ThreadPool> tpool, int seed)
    : id_(id),
      security_param_(security_param),
      rgen_(id, seed),
      network_(std::move(network)),
      preproc_(std::move(preproc)),
      circ_(std::move(circ)),
      tpool_(std::move(tpool)),
      wires_(circ.num_gates),
      jump_(id) {}

void OnlineEvaluator::setInputs(
    const std::unordered_map<utils::wire_t, Ring>& inputs) {
  // Input gates have depth 0.
  std::vector<Ring> my_betas;
  std::vector<size_t> num_inp_pid(4, 0);

  for (auto& g : circ_.gates_by_level[0]) {
    if (g->type == utils::GateType::kInp) {
      auto* pre_input =
          static_cast<PreprocInput<Ring>*>(preproc_.gates[g->out].get());
      auto pid = pre_input->pid;

      num_inp_pid[pid]++;
      if (pid == id_) {
        my_betas.push_back(pre_input->mask_value + inputs.at(g->out));
      }
    }
  }

  int prev_pid = pidFromOffset(id_, -1);
  std::vector<Ring> prev_betas(num_inp_pid[prev_pid]);
  std::vector<std::future<void>> res;

  // Send betas to next party.
  if (!my_betas.empty()) {
    res.push_back(tpool_->enqueue([&]() {
      network_->sendRelative(1, my_betas.data(),
                             my_betas.size() * sizeof(Ring));
      network_->flush(pidFromOffset(id_, 1));
    }));
  }

  // Receive betas from previous party.
  if (num_inp_pid[prev_pid] != 0) {
    res.push_back(tpool_->enqueue([&]() {
      network_->recv(prev_pid, prev_betas.data(),
                     prev_betas.size() * sizeof(Ring));
    }));
  }

  for (auto& f : res) {
    f.get();
  }

  // Run jump to receive all inputs.
  for (int pid = 0; pid < 4; ++pid) {
    if (num_inp_pid[pid] != 0) {
      if (pid == id_) {
        auto next_pid = pidFromOffset(id_, 1);
        jump_.jumpUpdate(id_, next_pid, pidFromOffset(id_, 2),
                         my_betas.size() * sizeof(Ring), my_betas.data());
        jump_.jumpUpdate(id_, next_pid, pidFromOffset(id_, 3),
                         my_betas.size() * sizeof(Ring), my_betas.data());
      } else if (pid == prev_pid) {
        jump_.jumpUpdate(id_, prev_pid, pidFromOffset(prev_pid, 2),
                         prev_betas.size() * sizeof(Ring), prev_betas.data());
        jump_.jumpUpdate(id_, prev_pid, pidFromOffset(prev_pid, 3),
                         prev_betas.size() * sizeof(Ring), prev_betas.data());
      } else {
        jump_.jumpUpdate(pid, pidFromOffset(pid, 1), id_,
                         num_inp_pid[pid] * sizeof(Ring));
      }
    }
  }
  jump_.communicate(*network_, *tpool_);

  std::vector<size_t> pid_inp_idx(4, 0);
  for (auto& g : circ_.gates_by_level[0]) {
    if (g->type == utils::GateType::kInp) {
      auto* pre_input =
          static_cast<PreprocInput<Ring>*>(preproc_.gates[g->out].get());
      auto pid = pre_input->pid;

      if (pid == id_) {
        wires_[g->out] = my_betas[pid_inp_idx[pid]];
      } else if (pid == prev_pid) {
        wires_[g->out] = prev_betas[pid_inp_idx[pid]];
      } else {
        const auto* values = reinterpret_cast<const Ring*>(
            jump_.getValues(pid, pidFromOffset(pid, 1)).data());
        wires_[g->out] = values[pid_inp_idx[pid]];
      }
      pid_inp_idx[pid]++;
    }
  }
  jump_.reset();
}

void OnlineEvaluator::setRandomInputs() {
  // Input gates have depth 0.
  for (auto& g : circ_.gates_by_level[0]) {
    if (g->type == utils::GateType::kInp) {
      rgen_.all().random_data(&wires_[g->out], sizeof(Ring));
    }
  }
}

std::array<std::vector<Ring>, 3> OnlineEvaluator::reluEvaluate(
    const std::vector<utils::FIn1Gate>& relu_gates) {
  auto num_relu_gates = relu_gates.size();
  std::vector<preprocg_ptr_t<BoolRing>*> vpreproc(num_relu_gates);

  // Iterate through preproc_ and extract info of relu gates.
  std::vector<utils::wire_t> win(num_relu_gates);
  for (size_t i = 0; i < num_relu_gates; ++i) {
    auto* pre_relu = static_cast<PreprocReluGate<Ring>*>(
        preproc_.gates[relu_gates[i].out].get());
    vpreproc[i] = pre_relu->msb_gates.data();
  }

  BoolEvaluator bool_eval(id_, vpreproc, msb_circ_);

  // Set the inputs.
  for (size_t i = 0; i < num_relu_gates; ++i) {
    auto val = wires_[relu_gates[i].in];

    auto val_bits = bitDecompose(val);
    for (size_t j = 0; j < msb_circ_.gates_by_level[0].size(); ++j) {
      const auto& gate = msb_circ_.gates_by_level[0][j];

      if (gate->type == utils::GateType::kInp) {
        bool_eval.vwires[i][gate->out] = 0;
        if (gate->out > 63) {
          bool_eval.vwires[i][gate->out] = val_bits[j - 64];
        }
      }
    }
  }

  bool_eval.evaluateAllLevels(*network_, jump_, *tpool_);
  auto output_shares = bool_eval.getOutputShares();

  std::vector<Ring> output_share_val(num_relu_gates);
  for (size_t i = 0; i < num_relu_gates; ++i) {
    if (output_shares[i][0].val()) {
      output_share_val[i] = 0;
    } else {
      output_share_val[i] = 1;
    }
  }

  // Bit to A.
  std::array<std::vector<Ring>, 3> relu_recon_shares;
  for (size_t i = 0; i < num_relu_gates; ++i) {
    auto* pre_relu = static_cast<PreprocReluGate<Ring>*>(
        preproc_.gates[relu_gates[i].out].get());
    auto beta_w = pre_relu->mask_w + (pre_relu->mask_msb * output_share_val[i]);
    for (int j = 0; j < 3; ++j) {
      relu_recon_shares[j].push_back(beta_w[j]);
    }
  }

  auto vw = reconstruct(relu_recon_shares);
  std::vector<Ring> vbtoa(num_relu_gates);
  for (size_t i = 0; i < num_relu_gates; ++i) {
    vbtoa[i] = vw[i] * static_cast<Ring>(-2) + output_share_val[i];
  }

  // Bit injection.
  std::vector<ReplicatedShare<Ring>> beta_inj(num_relu_gates);
  for (size_t i = 0; i < num_relu_gates; ++i) {
    auto* pre_relu = static_cast<PreprocReluGate<Ring>*>(
        preproc_.gates[relu_gates[i].out].get());

    auto win = relu_gates[i].in;
    auto mask_in = preproc_.gates[win]->mask;
    auto beta_v = wires_[win];

    beta_inj[i] = pre_relu->mask_binj + pre_relu->mask;
    beta_inj[i] -= (pre_relu->mask_btoa * beta_v + mask_in * vbtoa[i]);
    beta_inj[i].add(vbtoa[i] * beta_v, id_);
  }

  std::array<std::vector<Ring>, 3> outputs;
  for (const auto& share : beta_inj) {
    for (int i = 0; i < 3; ++i) {
      outputs[i].push_back(share[i]);
    }
  }

  return outputs;
}

std::array<std::vector<Ring>, 3> OnlineEvaluator::msbEvaluate(
    const std::vector<utils::FIn1Gate>& msb_gates) {
  auto num_msb_gates = msb_gates.size();
  std::vector<preprocg_ptr_t<BoolRing>*> vpreproc(num_msb_gates);

  // Iterate through preproc_ and extract info of msb gates.
  std::vector<utils::wire_t> win(num_msb_gates);
  for (size_t i = 0; i < num_msb_gates; ++i) {
    auto* pre_msb = static_cast<PreprocMsbGate<Ring>*>(
        preproc_.gates[msb_gates[i].out].get());
    vpreproc[i] = pre_msb->msb_gates.data();
  }

  BoolEvaluator bool_eval(id_, vpreproc, msb_circ_);

  // Set the inputs.
  for (size_t i = 0; i < num_msb_gates; ++i) {
    auto val = wires_[msb_gates[i].in];

    auto val_bits = bitDecompose(val);
    for (size_t j = 0; j < msb_circ_.gates_by_level[0].size(); ++j) {
      const auto& gate = msb_circ_.gates_by_level[0][j];

      if (gate->type == utils::GateType::kInp) {
        bool_eval.vwires[i][gate->out] = 0;
        if (gate->out > 63) {
          bool_eval.vwires[i][gate->out] = val_bits[j - 64];
        }
      }
    }
  }

  bool_eval.evaluateAllLevels(*network_, jump_, *tpool_);
  auto output_shares = bool_eval.getOutputShares();

  std::vector<Ring> output_share_val(num_msb_gates);
  for (size_t i = 0; i < num_msb_gates; ++i) {
    if (output_shares[i][0].val()) {
      output_share_val[i] = 1;
    } else {
      output_share_val[i] = 0;
    }
  }

  // Bit to A.
  std::array<std::vector<Ring>, 3> outputs;
  for (size_t i = 0; i < num_msb_gates; ++i) {
    auto* pre_msb = static_cast<PreprocMsbGate<Ring>*>(
        preproc_.gates[msb_gates[i].out].get());
    auto beta_w = pre_msb->mask_w + (pre_msb->mask_msb * output_share_val[i]);
    beta_w *= static_cast<Ring>(-2);
    beta_w.add(output_share_val[i], id_);

    for (int j = 0; j < 3; ++j) {
      outputs[j].push_back(beta_w[j]);
    }
  }

  return outputs;
}

std::vector<Ring> OnlineEvaluator::reconstruct(
    const std::array<std::vector<Ring>, 3>& recon_shares) {
  // All vectors in recon_shares should have same size.
  size_t num = recon_shares[0].size();
  size_t nbytes = sizeof(Ring) * num;

  if (nbytes == 0) {
    return {};
  }

  std::vector<Ring> vres(num);
  switch (id_) {
    case 0: {
      // Round 1.
      jump_.jumpUpdate(0, 1, 2, nbytes,
                       recon_shares[offsetFromPid(id_, 1) - 1].data());
      jump_.communicate(*network_, *tpool_);
      jump_.reset();

      // Round 2.
      jump_.jumpUpdate(0, 3, 2, nbytes,
                       recon_shares[offsetFromPid(id_, 3) - 1].data());
      jump_.jumpUpdate(0, 2, 3, nbytes,
                       recon_shares[offsetFromPid(id_, 2) - 1].data());
      jump_.communicate(*network_, *tpool_);
      jump_.reset();

      // Round 3.
      jump_.jumpUpdate(2, 3, 0, nbytes);
      jump_.communicate(*network_, *tpool_);
      const auto* vresj =
          reinterpret_cast<const Ring*>(jump_.getValues(2, 3).data());
      std::copy(vresj, vresj + num, vres.begin());
      jump_.reset();
      break;
    }

    case 1: {
      // Round 1.
      jump_.jumpUpdate(0, 1, 2, nbytes,
                       recon_shares[offsetFromPid(id_, 0) - 1].data());
      jump_.communicate(*network_, *tpool_);
      jump_.reset();

      const auto& v01 = recon_shares[offsetFromPid(id_, 0) - 1];
      const auto& v12 = recon_shares[offsetFromPid(id_, 2) - 1];
      for (int i = 0; i < num; ++i) {
        vres[i] = v01[i] + v12[i];
      }

      // Round 2.
      jump_.jumpUpdate(1, 3, 2, nbytes,
                       recon_shares[offsetFromPid(id_, 3) - 1].data());
      jump_.jumpUpdate(1, 2, 3, nbytes, vres.data());
      jump_.communicate(*network_, *tpool_);
      jump_.reset();

      // Round 3.
      jump_.jumpUpdate(2, 3, 1, nbytes);
      jump_.communicate(*network_, *tpool_);
      const auto* vresj =
          reinterpret_cast<const Ring*>(jump_.getValues(2, 3).data());
      std::copy(vresj, vresj + num, vres.begin());
      jump_.reset();
      break;
    }

    case 2: {
      // Round 1.
      jump_.jumpUpdate(0, 1, 2, nbytes);
      jump_.communicate(*network_, *tpool_);

      const auto* v01 =
          reinterpret_cast<const Ring*>(jump_.getValues(0, 1).data());
      const auto& v12 = recon_shares[offsetFromPid(id_, 1) - 1];
      for (int i = 0; i < num; ++i) {
        vres[i] = v01[i] + v12[i];
      }
      jump_.reset();

      // Round 2.
      jump_.jumpUpdate(0, 2, 3, nbytes,
                       recon_shares[offsetFromPid(id_, 0) - 1].data());
      jump_.jumpUpdate(1, 2, 3, nbytes, vres.data());
      jump_.jumpUpdate(0, 3, 2, nbytes);
      jump_.jumpUpdate(1, 3, 2, nbytes);
      jump_.communicate(*network_, *tpool_);

      const auto& v02 = recon_shares[offsetFromPid(id_, 0) - 1];
      const auto& v23 = recon_shares[offsetFromPid(id_, 3) - 1];
      const auto* v03 =
          reinterpret_cast<const Ring*>(jump_.getValues(0, 3).data());
      const auto* v13 =
          reinterpret_cast<const Ring*>(jump_.getValues(1, 3).data());
      for (int i = 0; i < num; ++i) {
        vres[i] += v02[i] + v23[i] + v03[i] + v13[i];
      }
      jump_.reset();

      // Round 3.
      // Jump 4.
      jump_.jumpUpdate(2, 3, 0, nbytes, vres.data());
      jump_.jumpUpdate(2, 3, 1, nbytes, vres.data());
      jump_.communicate(*network_, *tpool_);
      jump_.reset();
      break;
    }

    case 3: {
      jump_.communicate(*network_, *tpool_);
      jump_.reset();

      // Round 2.
      jump_.jumpUpdate(0, 3, 2, nbytes,
                       recon_shares[offsetFromPid(id_, 0) - 1].data());
      jump_.jumpUpdate(1, 3, 2, nbytes,
                       recon_shares[offsetFromPid(id_, 1) - 1].data());
      jump_.jumpUpdate(0, 2, 3, nbytes);
      jump_.jumpUpdate(1, 2, 3, nbytes);
      jump_.communicate(*network_, *tpool_);

      const auto& v03 = recon_shares[offsetFromPid(id_, 0) - 1];
      const auto& v13 = recon_shares[offsetFromPid(id_, 1) - 1];
      const auto& v23 = recon_shares[offsetFromPid(id_, 2) - 1];
      const auto* v01_12 =
          reinterpret_cast<const Ring*>(jump_.getValues(1, 2).data());
      const auto* v02 =
          reinterpret_cast<const Ring*>(jump_.getValues(0, 2).data());
      for (int i = 0; i < num; ++i) {
        vres[i] = v03[i] + v13[i] + v23[i] + v01_12[i] + v02[i];
      }
      jump_.reset();

      // Round 3.
      // Jump 4.
      jump_.jumpUpdate(2, 3, 0, nbytes, vres.data());
      jump_.jumpUpdate(2, 3, 1, nbytes, vres.data());
      jump_.communicate(*network_, *tpool_);
      jump_.reset();
      break;
    }
  }

  return vres;
}

void OnlineEvaluator::evaluateGatesAtDepth(size_t depth) {
  std::array<std::vector<Ring>, 3> recon_shares;
  std::vector<utils::FIn1Gate> relu_gates;
  std::vector<utils::FIn1Gate> msb_gates;

  for (auto& gate : circ_.gates_by_level[depth]) {
    switch (gate->type) {
      case utils::GateType::kMul: {
        auto* g = static_cast<utils::FIn2Gate*>(gate.get());
        auto& m_in1 = preproc_.gates[g->in1]->mask;
        auto& m_in2 = preproc_.gates[g->in2]->mask;
        auto* pre_out =
            static_cast<PreprocMultGate<Ring>*>(preproc_.gates[g->out].get());

        auto rec_share = pre_out->mask + pre_out->mask_prod -
                         m_in1 * wires_[g->in2] - m_in2 * wires_[g->in1];
        rec_share.add(wires_[g->in1] * wires_[g->in2], id_);

        for (int i = 0; i < 3; ++i) {
          recon_shares[i].push_back(rec_share[i]);
        }
        break;
      }

      case utils::GateType::kDotprod: {
        auto* g = static_cast<utils::SIMDGate*>(gate.get());
        auto* pre_out =
            static_cast<PreprocDotpGate<Ring>*>(preproc_.gates[g->out].get());

        auto rec_share = pre_out->mask + pre_out->mask_prod;
        for (size_t i = 0; i < g->in1.size(); i++) {
          auto win1 = g->in1[i];
          auto win2 = g->in2[i];
          auto& m_in1 = preproc_.gates[win1]->mask;
          auto& m_in2 = preproc_.gates[win2]->mask;

          rec_share -= m_in1 * wires_[win2] + m_in2 * wires_[win1];
          rec_share.add(wires_[win1] * wires_[win2], id_);
        }

        for (int i = 0; i < 3; i++) {
          recon_shares[i].push_back(rec_share[i]);
        }
        break;
      }

      case utils::GateType::kTrdotp: {
        auto* g = static_cast<utils::SIMDGate*>(gate.get());
        auto* pre_out =
            static_cast<PreprocTrDotpGate<Ring>*>(preproc_.gates[g->out].get());

        auto rec_share = pre_out->mask_prod + pre_out->mask_d;
        for (size_t i = 0; i < g->in1.size(); i++) {
          auto win1 = g->in1[i];
          auto win2 = g->in2[i];
          auto& m_in1 = preproc_.gates[win1]->mask;
          auto& m_in2 = preproc_.gates[win2]->mask;

          rec_share -= (m_in1 * wires_[win2] + m_in2 * wires_[win1]);
          rec_share.add(wires_[win1] * wires_[win2], id_);
        }

        for (int i = 0; i < 3; i++) {
          recon_shares[i].push_back(rec_share[i]);
        }

        break;
      }

      case utils::GateType::kRelu: {
        auto* g = static_cast<utils::FIn1Gate*>(gate.get());
        relu_gates.push_back(*g);
        break;
      }

      case utils::GateType::kMsb: {
        auto* g = static_cast<utils::FIn1Gate*>(gate.get());
        msb_gates.push_back(*g);
        break;
      }

      default:
        break;
    }
  }

  size_t non_relu_recon = recon_shares[0].size();
  if (!relu_gates.empty()) {
    auto shares = reluEvaluate(relu_gates);
    for (size_t i = 0; i < 3; ++i) {
      recon_shares[i].insert(recon_shares[i].end(), shares[i].begin(),
                             shares[i].end());
    }
  }

  if (!msb_gates.empty()) {
    auto shares = msbEvaluate(msb_gates);
    for (size_t i = 0; i < 3; ++i) {
      recon_shares[i].insert(recon_shares[i].end(), shares[i].begin(),
                             shares[i].end());
    }
  }

  auto vres = reconstruct(recon_shares);

  size_t idx = 0;
  size_t relu_idx = 0;
  size_t msb_idx = 0;
  for (auto& gate : circ_.gates_by_level[depth]) {
    switch (gate->type) {
      case utils::GateType::kAdd: {
        auto* g = static_cast<utils::FIn2Gate*>(gate.get());
        wires_[g->out] = wires_[g->in1] + wires_[g->in2];
        break;
      }

      case utils::GateType::kSub: {
        auto* g = static_cast<utils::FIn2Gate*>(gate.get());
        wires_[g->out] = wires_[g->in1] - wires_[g->in2];
        break;
      }

      case utils::GateType::kMul:
      case utils::GateType::kDotprod: {
        wires_[gate->out] = vres[idx++];
        break;
      }

      case utils::GateType::kRelu: {
        wires_[gate->out] = vres[non_relu_recon + relu_idx];
        relu_idx++;
        break;
      }

      case utils::GateType::kMsb: {
        wires_[gate->out] = vres[non_relu_recon + relu_idx + msb_idx];
        msb_idx++;
        break;
      }

      case utils::GateType::kTrdotp: {
        wires_[gate->out] = vres[idx++] >> FRACTION;
        break;
      }

      case utils::GateType::kConstAdd: {
        auto* g = static_cast<utils::ConstOpGate<Ring>*>(gate.get());
        wires_[g->out] = wires_[g->in] + g->cval;
        break;
      }

      case utils::GateType::kConstMul: {
        auto* g = static_cast<utils::ConstOpGate<Ring>*>(gate.get());
        wires_[g->out] = wires_[g->in] * g->cval;
        break;
      }

      default:
        break;
    }
  }
}

std::vector<Ring> OnlineEvaluator::reconstruct(
    const std::vector<ReplicatedShare<Ring>>& shares) {
  std::array<std::vector<Ring>, 3> recon_shares;
  for (const auto& s : shares) {
    for (size_t i = 0; i < 3; ++i) {
      recon_shares[i].push_back(s[i]);
    }
  }

  return reconstruct(recon_shares);
}

std::vector<Ring> OnlineEvaluator::getOutputs() {
  std::vector<Ring> outvals(circ_.outputs.size());
  if (circ_.outputs.empty()) {
    return outvals;
  }

  std::vector<std::future<void>> res;

  // Send openings to share commitments.
  for (int receiver = 0; receiver < 4; ++receiver) {
    if (receiver == id_) {
      continue;
    }

    res.push_back(tpool_->enqueue([&, receiver]() {
      for (int other = 0; other < 4; ++other) {
        if (other == id_ || other == receiver) {
          continue;
        }

        for (size_t i = 0; i < circ_.outputs.size(); ++i) {
          auto wout = circ_.outputs[i];
          Ring share_element =
              preproc_.gates[wout]->mask.commonValueWithParty(id_, other);
          const auto& opening =
              preproc_.output[i].openings[offsetFromPid(id_, other) - 1];
          network_->send(receiver, &share_element, sizeof(Ring));
          network_->send(receiver, opening.data(), opening.size());
        }
      }

      network_->flush(receiver);
    }));
  }

  // Receive share openings.
  std::vector<std::array<std::array<Ring, 3>, 3>> recv_shares(
      circ_.outputs.size());
  for (int sender = 0; sender < 4; ++sender) {
    if (sender == id_) {
      continue;
    }

    res.push_back(tpool_->enqueue([&, sender]() {
      std::vector<uint8_t> opening(security_param_ / 8);
      std::array<char, emp::Hash::DIGEST_SIZE> digest{};
      emp::Hash hash;

      for (int other = 0; other < 4; ++other) {
        if (other == id_ || other == sender) {
          continue;
        }

        for (size_t i = 0; i < circ_.outputs.size(); ++i) {
          Ring share_element = 0;
          network_->recv(sender, &share_element, sizeof(Ring));
          network_->recv(sender, opening.data(), opening.size());

          hash.reset();
          hash.put(&share_element, sizeof(Ring));
          hash.put(opening.data(), opening.size());
          hash.digest(digest.data());

          auto sender_fs = offsetFromPid(id_, sender) - 1;
          auto other_fs = offsetFromPid(id_, other) - 1;
          const auto& commitment =
              preproc_.output[i].commitments[sender_fs + other_fs - 1];
          if (commitment == digest) {
            recv_shares[i][sender_fs][other_fs] = share_element;
          } else {
            recv_shares[i][sender_fs][other_fs] = 0;
          }
        }
      }
    }));
  }

  for (auto& f : res) {
    f.get();
  }

  for (size_t i = 0; i < outvals.size(); ++i) {
    auto wout = circ_.outputs[i];
    outvals[i] = wires_[wout] - preproc_.gates[wout]->mask.sum();

    for (int p1 = 0; p1 < 3; ++p1) {
      for (int p2 = p1 + 1; p2 < 3; ++p2) {
        outvals[i] -= recv_shares[i][p1][p2];
        if (recv_shares[i][p1][p2] == 0) {
          outvals[i] -= recv_shares[i][p2][p1];
        }
      }
    }
  }

  return outvals;
}

std::vector<Ring> OnlineEvaluator::evaluateCircuit(
    const std::unordered_map<utils::wire_t, Ring>& inputs) {
  setInputs(inputs);
  for (size_t i = 0; i < circ_.gates_by_level.size(); ++i) {
    evaluateGatesAtDepth(i);
  }
  return getOutputs();
}

BoolEvaluator::BoolEvaluator(int my_id,
                             std::vector<preprocg_ptr_t<BoolRing>*> vpreproc,
                             utils::LevelOrderedCircuit circ)
    : id(my_id),
      vwires(vpreproc.size(), std::vector<BoolRing>(circ.num_gates)),
      vpreproc(std::move(vpreproc)),
      circ(std::move(circ)) {}

std::vector<BoolRing> BoolEvaluator::reconstruct(
    int id, const std::array<std::vector<BoolRing>, 3>& recon_shares,
    io::NetIOMP<4>& network, JumpProvider& jump, ThreadPool& tpool) {
  size_t num = recon_shares[0].size();

  std::array<std::vector<uint8_t>, 3> packed_recon_shares;
  for (int i = 0; i < 3; ++i) {
    packed_recon_shares[i] = BoolRing::pack(recon_shares[i].data(), num);
  }
  size_t nbytes = sizeof(uint8_t) * packed_recon_shares[0].size();

  if (nbytes == 0) {
    return {};
  }

  std::vector<BoolRing> vres;
  switch (id) {
    case 0: {
      // Round 1.
      jump.jumpUpdate(0, 1, 2, nbytes,
                      packed_recon_shares[offsetFromPid(id, 1) - 1].data());
      jump.communicate(network, tpool);
      jump.reset();

      // Round 2.
      jump.jumpUpdate(0, 3, 2, nbytes,
                      packed_recon_shares[offsetFromPid(id, 3) - 1].data());
      jump.jumpUpdate(0, 2, 3, nbytes,
                      packed_recon_shares[offsetFromPid(id, 2) - 1].data());
      jump.communicate(network, tpool);
      jump.reset();

      // Round 3.
      jump.jumpUpdate(2, 3, 0, nbytes);
      jump.communicate(network, tpool);
      vres = BoolRing::unpack(jump.getValues(2, 3).data(), num);
      jump.reset();
      break;
    }

    case 1: {
      // Round 1.
      jump.jumpUpdate(0, 1, 2, nbytes,
                      packed_recon_shares[offsetFromPid(id, 0) - 1].data());
      jump.communicate(network, tpool);
      jump.reset();

      std::vector<uint8_t> vtemp(packed_recon_shares[0].size());
      const auto& v01 = packed_recon_shares[offsetFromPid(id, 0) - 1];
      const auto& v12 = packed_recon_shares[offsetFromPid(id, 2) - 1];
      for (int i = 0; i < vtemp.size(); ++i) {
        vtemp[i] = v01[i] ^ v12[i];
      }

      // Round 2.
      jump.jumpUpdate(1, 3, 2, nbytes,
                      packed_recon_shares[offsetFromPid(id, 3) - 1].data());
      jump.jumpUpdate(1, 2, 3, nbytes, vtemp.data());
      jump.communicate(network, tpool);
      jump.reset();

      // Round 3.
      jump.jumpUpdate(2, 3, 1, nbytes);
      jump.communicate(network, tpool);
      vres = BoolRing::unpack(jump.getValues(2, 3).data(), num);
      jump.reset();
      break;
    }

    case 2: {
      // Round 1.
      jump.jumpUpdate(0, 1, 2, nbytes);
      jump.communicate(network, tpool);

      std::vector<uint8_t> vtemp(packed_recon_shares[0].size());
      const auto& v01 = jump.getValues(0, 1);
      const auto& v12 = packed_recon_shares[offsetFromPid(id, 1) - 1];
      for (int i = 0; i < vtemp.size(); ++i) {
        vtemp[i] = v01[i] ^ v12[i];
      }
      jump.reset();

      // Round 2.
      jump.jumpUpdate(0, 2, 3, nbytes,
                      packed_recon_shares[offsetFromPid(id, 0) - 1].data());
      jump.jumpUpdate(1, 2, 3, nbytes, vtemp.data());
      jump.jumpUpdate(0, 3, 2, nbytes);
      jump.jumpUpdate(1, 3, 2, nbytes);
      jump.communicate(network, tpool);

      const auto& v02 = packed_recon_shares[offsetFromPid(id, 0) - 1];
      const auto& v23 = packed_recon_shares[offsetFromPid(id, 3) - 1];
      const auto& v03 = jump.getValues(0, 3).data();
      const auto& v13 = jump.getValues(1, 3).data();
      for (int i = 0; i < vtemp.size(); ++i) {
        vtemp[i] = vtemp[i] ^ v02[i] ^ v23[i] ^ v03[i] ^ v13[i];
      }
      jump.reset();

      // Round 3.
      // Jump 4.
      jump.jumpUpdate(2, 3, 0, nbytes, vtemp.data());
      jump.jumpUpdate(2, 3, 1, nbytes, vtemp.data());
      jump.communicate(network, tpool);
      jump.reset();
      vres = BoolRing::unpack(vtemp.data(), num);
      break;
    }

    case 3: {
      jump.communicate(network, tpool);
      jump.reset();

      // Round 2.
      jump.jumpUpdate(0, 3, 2, nbytes,
                      packed_recon_shares[offsetFromPid(id, 0) - 1].data());
      jump.jumpUpdate(1, 3, 2, nbytes,
                      packed_recon_shares[offsetFromPid(id, 1) - 1].data());
      jump.jumpUpdate(0, 2, 3, nbytes);
      jump.jumpUpdate(1, 2, 3, nbytes);
      jump.communicate(network, tpool);

      std::vector<uint8_t> vtemp(packed_recon_shares[0].size());
      const auto& v03 = packed_recon_shares[offsetFromPid(id, 0) - 1];
      const auto& v13 = packed_recon_shares[offsetFromPid(id, 1) - 1];
      const auto& v23 = packed_recon_shares[offsetFromPid(id, 2) - 1];
      const auto& v01_12 = jump.getValues(1, 2).data();
      const auto& v02 = jump.getValues(0, 2).data();
      for (int i = 0; i < vtemp.size(); ++i) {
        vtemp[i] = v03[i] ^ v13[i] ^ v23[i] ^ v01_12[i] ^ v02[i];
      }
      jump.reset();

      // Round 3.
      // Jump 4.
      jump.jumpUpdate(2, 3, 0, nbytes, vtemp.data());
      jump.jumpUpdate(2, 3, 1, nbytes, vtemp.data());
      jump.communicate(network, tpool);
      jump.reset();
      vres = BoolRing::unpack(vtemp.data(), num);
      break;
    }
  }

  return vres;
}

void BoolEvaluator::evaluateGatesAtDepth(size_t depth, io::NetIOMP<4>& network,
                                         JumpProvider& jump,
                                         ThreadPool& tpool) {
  std::array<std::vector<BoolRing>, 3> recon_shares;
  for (size_t i = 0; i < vwires.size(); ++i) {
    const auto& preproc = vpreproc[i];
    auto& wires = vwires[i];

    for (auto& gate : circ.gates_by_level[depth]) {
      switch (gate->type) {
        case utils::GateType::kMul: {
          auto* g = static_cast<utils::FIn2Gate*>(gate.get());
          auto& m_in1 = preproc[g->in1]->mask;
          auto& m_in2 = preproc[g->in2]->mask;
          auto* pre_out =
              static_cast<PreprocMultGate<BoolRing>*>(preproc[g->out].get());

          auto rec_share = pre_out->mask + pre_out->mask_prod -
                           m_in1 * wires[g->in2] - m_in2 * wires[g->in1];
          rec_share.add(wires[g->in1] * wires[g->in2], id);

          for (int i = 0; i < 3; ++i) {
            recon_shares[i].push_back(rec_share[i]);
          }
          break;
        }

        default:
          break;
      }
    }
  }

  auto vres = reconstruct(id, recon_shares, network, jump, tpool);

  // Update mult gate values.
  size_t idx = 0;
  for (auto& wires : vwires) {
    for (auto& gate : circ.gates_by_level[depth]) {
      switch (gate->type) {
        case utils::GateType::kAdd: {
          auto* g = static_cast<utils::FIn2Gate*>(gate.get());
          wires[g->out] = wires[g->in1] + wires[g->in2];
          break;
        }

        case utils::GateType::kSub: {
          auto* g = static_cast<utils::FIn2Gate*>(gate.get());
          wires[g->out] = wires[g->in1] - wires[g->in2];
          break;
        }

        case utils::GateType::kMul: {
          wires[gate->out] = vres[idx++];
          break;
        }

        case utils::GateType::kConstAdd: {
          auto* g = static_cast<utils::ConstOpGate<BoolRing>*>(gate.get());
          wires[g->out] = wires[g->in] + g->cval;
          break;
        }

        case utils::GateType::kConstMul: {
          auto* g = static_cast<utils::ConstOpGate<BoolRing>*>(gate.get());
          wires[g->out] = wires[g->in] * g->cval;
          break;
        }

        default:
          break;
      }
    }
  }
}

void BoolEvaluator::evaluateAllLevels(io::NetIOMP<4>& network,
                                      JumpProvider& jump, ThreadPool& tpool) {
  for (size_t i = 0; i < circ.gates_by_level.size(); ++i) {
    evaluateGatesAtDepth(i, network, jump, tpool);
  }
}

std::vector<std::vector<BoolRing>> BoolEvaluator::getOutputShares() {
  std::vector<std::vector<BoolRing>> outputs(
      vwires.size(), std::vector<BoolRing>(circ.outputs.size()));

  for (size_t i = 0; i < vwires.size(); ++i) {
    const auto& wires = vwires[i];
    for (size_t j = 0; j < circ.outputs.size(); ++j) {
      outputs[i][j] = wires[circ.outputs[j]];
    }
  }

  return outputs;
}
};  // namespace quadsquad
