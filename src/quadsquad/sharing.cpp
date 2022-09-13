#include "sharing.h"

namespace quadsquad {
template <>
void ReplicatedShare<BoolRing>::randomize(emp::PRG& prg) {
  bool values[3];
  prg.random_bool(static_cast<bool*>(values), 3);
  for (size_t i = 0; i < 3; ++i) {
    values_[i] = values[i];
  }
}

template <>
DummyShare<BoolRing>::DummyShare(BoolRing secret, emp::PRG& prg) {
  bool values[5];
  prg.random_bool(static_cast<bool*>(values), 5);

  BoolRing sum;
  for (size_t i = 0; i < 5; ++i) {
    share_elements[i] = values[i];
    sum += share_elements[i];
  }
  share_elements[5] = secret - sum;
}

template <>
void DummyShare<BoolRing>::randomize(emp::PRG& prg) {
  bool values[6];
  prg.random_bool(static_cast<bool*>(values), 6);

  for (size_t i = 0; i < 6; ++i) {
    share_elements[i] = values[i];
  }
}
};  // namespace quadsquad
