#pragma once

#include "../utils/circuit.h"
#include "sharing.h"
#include "types.h"

namespace quadsquad {
// Preprocessed data for a gate.
template <class R>
struct PreprocGate {
  // Secret shared mask for the output wire of the gate.
  ReplicatedShare<R> mask{};

  PreprocGate() = default;

  explicit PreprocGate(const ReplicatedShare<R>& mask) : mask(mask) {}

  virtual ~PreprocGate() = default;
};

template <class R>
using preprocg_ptr_t = std::unique_ptr<PreprocGate<R>>;

template <class R>
struct PreprocInput : public PreprocGate<R> {
  // ID of party providing input on wire.
  int pid{};
  // Plaintext value of mask on input wire. Non-zero for all parties except
  // party with id 'pid'.
  R mask_value{};

  PreprocInput() = default;
  PreprocInput(const ReplicatedShare<R>& mask, int pid, R mask_value = 0)
      : PreprocGate<R>(mask), pid(pid), mask_value(mask_value) {}
};

template <class R>
struct PreprocMultGate : public PreprocGate<R> {
  // Secret shared product of inputs masks.
  ReplicatedShare<R> mask_prod{};

  PreprocMultGate() = default;
  PreprocMultGate(const ReplicatedShare<R>& mask,
                  const ReplicatedShare<R>& mask_prod)
      : PreprocGate<R>(mask), mask_prod(mask_prod) {}
};

template <class R>
struct PreprocDotpGate : public PreprocGate<R> {
  ReplicatedShare<Ring> mask_prod{};

  PreprocDotpGate() = default;
  PreprocDotpGate(const ReplicatedShare<Ring>& mask,
                  const ReplicatedShare<Ring>& mask_prod)
      : PreprocGate<R>(mask), mask_prod(mask_prod) {}
};

template <class R>
struct PreprocTrDotpGate : public PreprocGate<R> {
  ReplicatedShare<Ring> mask_prod{};
  ReplicatedShare<Ring> mask_d{};

  PreprocTrDotpGate() = default;
  PreprocTrDotpGate(const ReplicatedShare<Ring>& mask,
                    const ReplicatedShare<Ring>& mask_prod,
                    const ReplicatedShare<Ring>& mask_d)
      : PreprocGate<R>(mask), mask_prod(mask_prod), mask_d(mask_d) {}
};

template <class R>
struct PreprocReluGate : public PreprocGate<R> {
  std::vector<preprocg_ptr_t<BoolRing>> msb_gates;
  ReplicatedShare<R> mask_msb;
  ReplicatedShare<R> mask_w;
  ReplicatedShare<R> mask_btoa;
  ReplicatedShare<R> mask_binj;

  PreprocReluGate() = default;
  PreprocReluGate(ReplicatedShare<R> mask,
                  std::vector<preprocg_ptr_t<BoolRing>> msb_gates,
                  ReplicatedShare<R> mask_msb, ReplicatedShare<R> mask_w,
                  ReplicatedShare<R> mask_btoa, ReplicatedShare<R> mask_binj)
      : PreprocGate<R>(mask),
        msb_gates(std::move(msb_gates)),
        mask_msb(mask_msb),
        mask_w(mask_w),
        mask_btoa(mask_btoa),
        mask_binj(mask_binj) {}
};

template <class R>
struct PreprocMsbGate : public PreprocGate<R> {
  std::vector<preprocg_ptr_t<BoolRing>> msb_gates;
  ReplicatedShare<R> mask_msb;
  ReplicatedShare<R> mask_w;

  PreprocMsbGate() = default;
  PreprocMsbGate(ReplicatedShare<R> mask,
                 std::vector<preprocg_ptr_t<BoolRing>> msb_gates,
                 ReplicatedShare<R> mask_msb, ReplicatedShare<R> mask_w)
      : PreprocGate<R>(mask),
        msb_gates(std::move(msb_gates)),
        mask_msb(mask_msb),
        mask_w(mask_w) {}
};

// Preprocessed data for output wires.
struct PreprocOutput {
  // Commitment corresponding to share elements not available with the party
  // for the output wire. If party's ID is 'i' then array is of the form
  // {s[i+1, i+2], s[i+1, i+3], s[i+2, i+3]}.
  std::array<std::array<char, emp::Hash::DIGEST_SIZE>, 3> commitments{};

  // Opening info for commitments to party's output shares.
  // If party's ID is 'i' then array is of the form
  // {o[i, i+1], o[i, i+2], o[i, i+3]} where o[i, j] is the opening info for
  // share common to parties with ID i and j.
  std::array<std::vector<uint8_t>, 3> openings;
};

// Preprocessed data for the circuit.
template <class R>
struct PreprocCircuit {
  std::vector<preprocg_ptr_t<R>> gates;
  std::vector<PreprocOutput> output;

  PreprocCircuit() = default;
  PreprocCircuit(size_t num_gates, size_t num_output)
      : gates(num_gates), output(num_output) {}
};
};  // namespace quadsquad
