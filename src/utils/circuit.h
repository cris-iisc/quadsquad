#pragma once

#include <algorithm>
#include <array>
#include <boost/format.hpp>
#include <cmath>
#include <iostream>
#include <memory>
#include <optional>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include "../quadsquad/helpers.h"
#include "../quadsquad/types.h"

namespace quadsquad::utils {

using wire_t = size_t;

enum GateType {
  kInp,
  kAdd,
  kMul,
  kSub,
  kConstAdd,
  kConstMul,
  kRelu,
  kMsb,
  kDotprod,
  kTrdotp,
  kInvalid,
  NumGates
};

std::ostream& operator<<(std::ostream& os, GateType type);

// Gates represent primitive operations.
// All gates have one output.
struct Gate {
  GateType type{GateType::kInvalid};
  wire_t out;

  Gate() = default;
  Gate(GateType type, wire_t out);

  virtual ~Gate() = default;
};

// Represents a gate with fan-in 2.
struct FIn2Gate : public Gate {
  wire_t in1{0};
  wire_t in2{0};

  FIn2Gate() = default;
  FIn2Gate(GateType type, wire_t in1, wire_t in2, wire_t out);
};

// Represents a gate with fan-in 1.
struct FIn1Gate : public Gate {
  wire_t in{0};

  FIn1Gate() = default;
  FIn1Gate(GateType type, wire_t in, wire_t out);
};

// Represents a gate used to denote SIMD operations.
// These type is used to represent operations that take vectors of inputs but
// might not necessarily be SIMD e.g., dot product.
struct SIMDGate : public Gate {
  std::vector<wire_t> in1{0};
  std::vector<wire_t> in2{0};

  SIMDGate() = default;
  SIMDGate(GateType type, std::vector<wire_t> in1, std::vector<wire_t> in2,
           wire_t out);
};

// Represents gates where one input is a constant.
template <class R>
struct ConstOpGate : public Gate {
  wire_t in{0};
  R cval;

  ConstOpGate() = default;
  ConstOpGate(GateType type, wire_t in, R cval, wire_t out)
      : Gate(type, out), in(in), cval(std::move(cval)) {}
};

using gate_ptr_t = std::shared_ptr<Gate>;

// Gates ordered by multiplicative depth.
//
// Addition gates are not considered to increase the depth.
// Moreover, if gates_by_level[l][i]'s output is input to gates_by_level[l][j]
// then i < j.
struct LevelOrderedCircuit {
  size_t num_gates;
  std::array<uint64_t, GateType::NumGates> count;
  std::vector<wire_t> outputs;
  std::vector<std::vector<gate_ptr_t>> gates_by_level;

  friend std::ostream& operator<<(std::ostream& os,
                                  const LevelOrderedCircuit& circ);
};

// Represents an arithmetic circuit.
template <class R>
class Circuit {
  std::vector<wire_t> outputs_;
  std::vector<gate_ptr_t> gates_;

  bool isWireValid(wire_t wid) { return wid < gates_.size(); }

 public:
  Circuit() = default;

  // Methods to manually build a circuit.
  wire_t newInputWire() {
    wire_t wid = gates_.size();
    gates_.push_back(std::make_shared<Gate>(GateType::kInp, wid));
    return wid;
  }

  void setAsOutput(wire_t wid) {
    if (!isWireValid(wid)) {
      throw std::invalid_argument("Invalid wire ID.");
    }

    outputs_.push_back(wid);
  }

  // Function to add a gate with fan-in 2.
  wire_t addGate(GateType type, wire_t input1, wire_t input2) {
    if (type != GateType::kAdd && type != GateType::kMul &&
        type != GateType::kSub) {
      throw std::invalid_argument("Invalid gate type.");
    }

    if (!isWireValid(input1) || !isWireValid(input2)) {
      throw std::invalid_argument("Invalid wire ID.");
    }

    wire_t output = gates_.size();
    gates_.push_back(std::make_shared<FIn2Gate>(type, input1, input2, output));

    return output;
  }

  // Function to add a gate with one input from a wire and a second constant
  // input.
  wire_t addConstOpGate(GateType type, wire_t wid, R cval) {
    if (type != kConstAdd && type != kConstMul) {
      throw std::invalid_argument("Invalid gate type.");
    }

    if (!isWireValid(wid)) {
      throw std::invalid_argument("Invalid wire ID.");
    }

    wire_t output = gates_.size();
    gates_.push_back(std::make_shared<ConstOpGate<R>>(type, wid, cval, output));

    return output;
  }

  // Function to add a single input gate.
  wire_t addGate(GateType type, wire_t input) {
    if (type != GateType::kRelu && type != GateType::kMsb) {
      throw std::invalid_argument("Invalid gate type.");
    }

    if (!isWireValid(input)) {
      throw std::invalid_argument("Invalid wire ID.");
    }

    wire_t output = gates_.size();
    gates_.push_back(std::make_shared<FIn1Gate>(type, input, output));

    return output;
  }

  // Function to add a multiple fan-in gate.
  wire_t addGate(GateType type, const std::vector<wire_t>& input1,
                 const std::vector<wire_t>& input2) {
    if (type != GateType::kDotprod && type != GateType::kTrdotp) {
      throw std::invalid_argument("Invalid gate type.");
    }

    if (input1.size() != input2.size()) {
      throw std::invalid_argument("Expected same length inputs.");
    }

    for (size_t i = 0; i < input1.size(); ++i) {
      if (!isWireValid(input1[i]) || !isWireValid(input2[i])) {
        throw std::invalid_argument("Invalid wire ID.");
      }
    }

    wire_t output = gates_.size();
    gates_.push_back(std::make_shared<SIMDGate>(type, input1, input2, output));
    return output;
  }

  // Level ordered gates are helpful for evaluation.
  [[nodiscard]] LevelOrderedCircuit orderGatesByLevel() const {
    LevelOrderedCircuit res;
    res.outputs = outputs_;
    res.num_gates = gates_.size();

    // Map from output wire id to multiplicative depth/level.
    // Input gates have a depth of 0.
    std::vector<size_t> gate_level(res.num_gates, 0);
    size_t depth = 0;

    // This assumes that if gates_[i]'s output is input to gates_[j] then
    // i < j.
    for (const auto& gate : gates_) {
      switch (gate->type) {
        case GateType::kAdd:
        case GateType::kSub: {
          const auto* g = static_cast<FIn2Gate*>(gate.get());
          gate_level[g->out] = std::max(gate_level[g->in1], gate_level[g->in2]);
          break;
        }

        case GateType::kMul: {
          const auto* g = static_cast<FIn2Gate*>(gate.get());
          gate_level[g->out] =
              std::max(gate_level[g->in1], gate_level[g->in2]) + 1;
          break;
        }

        case GateType::kConstAdd:
        case GateType::kConstMul: {
          const auto* g = static_cast<ConstOpGate<R>*>(gate.get());
          gate_level[g->out] = gate_level[g->in];
          break;
        }

        case GateType::kRelu: {
          const auto* g = static_cast<FIn1Gate*>(gate.get());
          gate_level[g->out] = gate_level[g->in] + 1;
          break;
        }

        case GateType::kMsb: {
          const auto* g = static_cast<FIn1Gate*>(gate.get());
          gate_level[g->out] = gate_level[g->in] + 1;
          break;
        }

        case GateType::kDotprod:
        case GateType::kTrdotp: {
          const auto* g = static_cast<SIMDGate*>(gate.get());
          size_t gate_depth = 0;
          for (size_t i = 0; i < g->in1.size(); ++i) {
            gate_depth = std::max(
                {gate_level[g->in1[i]], gate_level[g->in2[i]], gate_depth});
          }
          gate_level[g->out] = gate_depth + 1;
          break;
        }

        default:
          break;
      }

      depth = std::max(depth, gate_level[gate->out]);
    }

    std::fill(res.count.begin(), res.count.end(), 0);

    std::vector<std::vector<gate_ptr_t>> gates_by_level(depth + 1);
    for (const auto& gate : gates_) {
      res.count[gate->type]++;
      gates_by_level[gate_level[gate->out]].push_back(gate);
    }

    res.gates_by_level = std::move(gates_by_level);

    return res;
  }

  // Evaluate circuit on plaintext inputs.
  [[nodiscard]] std::vector<R> evaluate(
      const std::unordered_map<wire_t, R>& inputs) const {
    auto level_circ = orderGatesByLevel();
    std::vector<R> wires(level_circ.num_gates);

    auto num_inp_gates = level_circ.count[GateType::kInp];
    if (inputs.size() != num_inp_gates) {
      throw std::invalid_argument(boost::str(
          boost::format("Expected %1% inputs but received %2% inputs.") %
          num_inp_gates % inputs.size()));
    }

    for (const auto& level : level_circ.gates_by_level) {
      for (const auto& gate : level) {
        switch (gate->type) {
          case GateType::kInp: {
            wires[gate->out] = inputs.at(gate->out);
            break;
          }

          case GateType::kMul: {
            auto* g = static_cast<FIn2Gate*>(gate.get());
            wires[g->out] = wires[g->in1] * wires[g->in2];
            break;
          }

          case GateType::kAdd: {
            auto* g = static_cast<FIn2Gate*>(gate.get());
            wires[g->out] = wires[g->in1] + wires[g->in2];
            break;
          }

          case GateType::kSub: {
            auto* g = static_cast<FIn2Gate*>(gate.get());
            wires[g->out] = wires[g->in1] - wires[g->in2];
            break;
          }

          case GateType::kConstAdd: {
            auto* g = static_cast<ConstOpGate<R>*>(gate.get());
            wires[g->out] = wires[g->in] + g->cval;
            break;
          }

          case GateType::kConstMul: {
            auto* g = static_cast<ConstOpGate<R>*>(gate.get());
            wires[g->out] = wires[g->in] * g->cval;
            break;
          }

          case GateType::kRelu: {
            // ReLU gates don't make sense for boolean rings.
            if constexpr (std::is_same_v<R, BoolRing>) {
              throw std::runtime_error("ReLU gates are invalid for BoolRing.");
            } else {
              auto* g = static_cast<FIn1Gate*>(gate.get());
              std::vector<BoolRing> bin = bitDecompose(wires[g->in]);

              if (bin[63].val())
                wires[g->out] = 0;
              else
                wires[g->out] = wires[g->in];
            }
            break;
          }

          case GateType::kMsb: {
            auto* g = static_cast<FIn1Gate*>(gate.get());

            if constexpr (std::is_same_v<R, BoolRing>) {
              wires[g->out] = wires[g->in];
            } else {
              std::vector<BoolRing> bin = bitDecompose(wires[g->in]);
              wires[g->out] = bin[63].val();
            }
            break;
          }

          case GateType::kDotprod: {
            auto* g = static_cast<SIMDGate*>(gate.get());
            for (size_t i = 0; i < g->in1.size(); i++) {
              wires[g->out] += wires[g->in1.at(i)] * wires[g->in2.at(i)];
            }
            break;
          }

          case GateType::kTrdotp: {
            // Truncation makes sense only for non-boolean rings.
            if constexpr (std::is_same_v<R, BoolRing>) {
              throw std::runtime_error(
                  "Truncation gates are invalid for BoolRing.");
            } else {
              auto* g = static_cast<SIMDGate*>(gate.get());
              for (size_t i = 0; i < g->in1.size(); i++) {
                auto temp = wires[g->in1.at(i)] * wires[g->in2.at(i)];
                wires[g->out] += temp;
              }

              wires[g->out] = wires[g->out] >> FRACTION;
            }
            break;
          }

          default: {
            throw std::runtime_error("Invalid gate type.");
          }
        }
      }
    }

    std::vector<R> outputs;
    for (auto i : level_circ.outputs) {
      outputs.push_back(wires[i]);
    }

    return outputs;
  }

  static Circuit generatePPA() {
    Circuit circ;
    std::vector<wire_t> input_a(64);
    std::vector<wire_t> input_b(64);

    std::vector<wire_t> loc_p, loc_g;
    for (int i = 0; i < 64; i++) {
      input_a[i] = circ.newInputWire();
    }

    for (int i = 0; i < 64; i++) {
      input_b[i] = circ.newInputWire();
    }

    // input_a[0] stores the lsb.
    for (int i = 0; i < 64; i++) {
      auto p_id = circ.addGate(GateType::kAdd, input_a[i], input_b[i]);
      loc_p.push_back(p_id);

      auto g_id = circ.addGate(GateType::kMul, input_a[i], input_b[i]);
      loc_g.push_back(g_id);
    }

    for (int level = 1; level <= 6; level++) {
      for (int count = 1; count <= 64 / std::pow(2, level); count++) {
        int temp =
            std::pow(2, level - 1) + (count - 1) * std::pow(2, level) - 1;
        for (int i = 0; i < std::pow(2, level - 1); i++) {
          auto w1 =
              circ.addGate(GateType::kMul, loc_p[temp + i + 1], loc_g[temp]);

          auto w2 = circ.addGate(GateType::kAdd, loc_g[temp + i + 1], w1);

          loc_g[temp + i + 1] = w2;
          auto w3 =
              circ.addGate(GateType::kMul, loc_p[temp + i + 1], loc_p[temp]);

          loc_p[temp + i + 1] = w3;
        }
      }
    }

    std::vector<wire_t> S;

    S.push_back(circ.addGate(GateType::kAdd, input_a[0], input_b[0]));
    for (int i = 1; i < 64; i++) {
      auto w = circ.addGate(GateType::kAdd, input_a[i], input_b[i]);
      S.push_back(circ.addGate(GateType::kAdd, w, loc_g[i - 1]));
    }

    for (int i = 0; i < 64; i++) {
      circ.setAsOutput(S[i]);
    }
    return circ;
  }

  static Circuit generatePPAMSB() {
    Circuit circ;
    std::vector<wire_t> input_a(64);
    std::vector<wire_t> input_b(64);

    std::vector<wire_t> loc_p, loc_g;
    for (int i = 0; i < 64; i++) {
      input_a[i] = circ.newInputWire();
    }

    for (int i = 0; i < 64; i++) {
      input_b[i] = circ.newInputWire();
    }

    // input_a[0] stores the lsb.
    for (int i = 0; i < 64; i++) {
      auto p_id = circ.addGate(GateType::kAdd, input_a[i], input_b[i]);
      loc_p.push_back(p_id);
      auto g_id = circ.addGate(GateType::kMul, input_a[i], input_b[i]);
      loc_g.push_back(g_id);
    }

    for (int level = 1; level <= 6; level++) {
      for (int count = 1; count <= 64 / std::pow(2, level); count++) {
        int temp =
            std::pow(2, level - 1) + (count - 1) * std::pow(2, level) - 1;
        int offset = std::pow(2, level - 1);
        if (count < 64 / std::pow(2, level)) {
          auto w1 =
              circ.addGate(GateType::kMul, loc_p[temp + offset], loc_g[temp]);
          auto w2 = circ.addGate(GateType::kAdd, loc_g[temp + offset], w1);
          loc_g[temp + offset] = w2;

          auto w3 =
              circ.addGate(GateType::kMul, loc_p[temp + offset], loc_p[temp]);
          loc_p[temp + offset] = w3;
        } else {
          if (level != 1) {
            auto w1 = circ.addGate(GateType::kMul, loc_p[62], loc_g[temp]);
            auto w2 = circ.addGate(GateType::kAdd, loc_g[62], w1);
            loc_g[62] = w2;

            auto w3 = circ.addGate(GateType::kMul, loc_p[62], loc_p[temp]);
            loc_p[62] = w3;
          }
        }
      }
    }
    auto w = circ.addGate(GateType::kAdd, input_a[63], input_b[63]);

    auto msb = circ.addGate(GateType::kAdd, w, loc_g[62]);

    circ.setAsOutput(msb);
    return circ;
  }
};
};  // namespace quadsquad::utils
