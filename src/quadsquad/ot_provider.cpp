#include "ot_provider.h"

#include <boost/format.hpp>
#include <memory>

namespace quadsquad {
constexpr int64_t ot_bsize = emp::ot_bsize;

OTProvider::OTProvider(int my_id, int other_id, emp::NetIO* send_io,
                       emp::NetIO* recv_io)
    : send_ios_{send_io}, recv_ios_{recv_io} {
  std::string send_filename = boost::str(
      boost::format("./ot_data/s%1%_r%2%_sender.bin") % my_id % other_id);

  std::string recv_filename = boost::str(
      boost::format("./ot_data/s%1%_r%2%_receiver.bin") % other_id % my_id);

  if (my_id < other_id) {
    send_ot_ = std::make_unique<FerretCOT<NetIO>>(ALICE, 1, send_ios_.data(),
                                                  true, true, send_filename);
    recv_ot_ = std::make_unique<FerretCOT<NetIO>>(BOB, 1, recv_ios_.data(),
                                                  true, true, recv_filename);
  } else {
    recv_ot_ = std::make_unique<FerretCOT<NetIO>>(BOB, 1, recv_ios_.data(),
                                                  true, true, recv_filename);
    send_ot_ = std::make_unique<FerretCOT<NetIO>>(ALICE, 1, send_ios_.data(),
                                                  true, true, send_filename);
  }
}

void OTProvider::send(const uint64_t* data0, const uint64_t* data1,
                      int64_t length) {
  auto* data = new emp::block[length];
  send_ot_->send_cot(data, length);
  emp::block s;
  send_ot_->prg.random_block(&s, 1);
  send_ios_[0]->send_block(&s, 1);
  send_ot_->mitccrh.setS(s);
  send_ios_[0]->flush();

  block pad[2 * ot_bsize];
  uint64_t upad[2 * ot_bsize];
  auto* tpad = reinterpret_cast<uint64_t*>(pad);
  for (int64_t i = 0; i < length; i += ot_bsize) {
    for (int64_t j = i; j < min(i + ot_bsize, length); ++j) {
      pad[2 * (j - i)] = data[j];
      pad[2 * (j - i) + 1] = data[j] ^ send_ot_->Delta;
    }
    send_ot_->mitccrh.hash<ot_bsize, 2>(pad);
    for (int64_t j = i; j < min(i + ot_bsize, length); ++j) {
      upad[2 * (j - i)] = tpad[4 * (j - i)] ^ data0[j];
      upad[2 * (j - i) + 1] = tpad[4 * (j - i) + 2] ^ data1[j];
    }
    send_ios_[0]->send_data(upad,
                            2 * sizeof(uint64_t) * min(ot_bsize, length - i));
  }
  delete[] data;
}

void OTProvider::recv(uint64_t* udata, const bool* r, int64_t length) {
  auto* data = new emp::block[length];
  recv_ot_->recv_cot(data, r, length);
  emp::block s;
  recv_ios_[0]->recv_block(&s, 1);
  recv_ot_->mitccrh.setS(s);
  recv_ios_[0]->flush();

  uint64_t res[2 * ot_bsize];
  block pad[ot_bsize];
  auto* tpad = reinterpret_cast<uint64_t*>(pad);
  for (int64_t i = 0; i < length; i += ot_bsize) {
    memcpy(pad, data + i, min(ot_bsize, length - i) * sizeof(block));
    recv_ot_->mitccrh.hash<ot_bsize, 1>(pad);
    recv_ios_[0]->recv_data(res,
                            2 * sizeof(uint64_t) * min(ot_bsize, length - i));
    for (int64_t j = 0; j < ot_bsize and j < length - i; ++j) {
      udata[i + j] = res[2 * j + r[i + j]] ^ tpad[2 * j];
    }
  }
  delete[] data;
}

std::vector<Ring> OTProvider::multiplySend(const std::vector<Ring>& inputs,
                                           emp::PRG& prg) {
  size_t num_bits = sizeof(Ring) * 8;
  size_t num_blocks = num_bits * inputs.size();

  std::vector<Ring> vrand(num_blocks);
  prg.random_data(vrand.data(), sizeof(Ring) * vrand.size());

  std::vector<uint64_t> inp_0(num_blocks);
  std::vector<uint64_t> inp_1(num_blocks);
  std::vector<Ring> shares(inputs.size());
  size_t idx = 0;
  for (size_t i = 0; i < inputs.size(); ++i) {
    const auto& input = inputs[i];
    auto& share = shares[i];
    for (size_t j = 0; j < num_bits; ++j) {
      auto val = vrand[idx];
      share -= val;
      inp_0[idx] = val;
      inp_1[idx] = (input << j) + val;
      idx++;
    }
  }

  send(inp_0.data(), inp_1.data(), num_blocks);
  send_ios_[0]->flush();
  return shares;
}

std::vector<Ring> OTProvider::multiplyReceive(const std::vector<Ring>& inputs) {
  size_t num_bits = sizeof(Ring) * 8;
  size_t num_blocks = num_bits * inputs.size();

  std::unique_ptr<bool[]> choice_bits(new bool[num_blocks]);
  size_t idx = 0;
  for (auto input : inputs) {
    for (size_t j = 0; j < num_bits; ++j) {
      bool lsb = ((input >> j) & 1U) == 1;
      choice_bits[idx++] = lsb;
    }
  }

  std::vector<uint64_t> recv_blocks(num_blocks);
  recv(recv_blocks.data(), choice_bits.get(), num_blocks);

  std::vector<Ring> shares(inputs.size(), 0);
  idx = 0;
  for (size_t i = 0; i < inputs.size(); ++i) {
    for (size_t j = 0; j < num_bits; ++j) {
      shares[i] += recv_blocks[idx];
      idx++;
    }
  }

  return shares;
}
};  // namespace quadsquad
