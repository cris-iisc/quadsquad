#pragma once

#include <emp-tool/emp-tool.h>

#include <array>
#include <vector>

#include "helpers.h"
#include "types.h"

namespace quadsquad {

template <class R>
class ReplicatedShare {
  // values_[i] will denote element common with party having my_id + i + 1.
  std::array<R, 3> values_;

 public:
  ReplicatedShare() = default;
  explicit ReplicatedShare(std::array<R, 3> values)
      : values_{std::move(values)} {}

  void randomize(emp::PRG& prg) {
    prg.random_data(values_.data(), sizeof(R) * 3);
  }

  // Access share elements.
  // idx = i retreives value common with party having my_id + i + 1.
  R& operator[](size_t idx) { return values_.at(idx); }

  R operator[](size_t idx) const { return values_.at(idx); }

  R& commonValueWithParty(int my_id, int pid) {
    int idx = pid - my_id;
    if (idx < 0) {
      idx += 4;
    }
    idx -= 1;

    return values_.at(idx);
  }

  [[nodiscard]] R commonValueWithParty(int my_id, int pid) const {
    int idx = pid - my_id;
    if (idx < 0) {
      idx += 4;
    }
    idx -= 1;

    return values_.at(idx);
  }

  [[nodiscard]] R sum() const { return values_[0] + values_[1] + values_[2]; }

  // Arithmetic operators.
  ReplicatedShare<R>& operator+=(const ReplicatedShare<R>& rhs) {
    values_[0] += rhs.values_[0];
    values_[1] += rhs.values_[1];
    values_[2] += rhs.values_[2];
    return *this;
  }

  friend ReplicatedShare<R> operator+(ReplicatedShare<R> lhs,
                                      const ReplicatedShare<R>& rhs) {
    lhs += rhs;
    return lhs;
  }

  ReplicatedShare<R>& operator-=(const ReplicatedShare<R>& rhs) {
    (*this) += (rhs * -1);
    return *this;
  }

  friend ReplicatedShare<R> operator-(ReplicatedShare<R> lhs,
                                      const ReplicatedShare<R>& rhs) {
    lhs -= rhs;
    return lhs;
  }

  ReplicatedShare<R>& operator*=(const R& rhs) {
    values_[0] *= rhs;
    values_[1] *= rhs;
    values_[2] *= rhs;
    return *this;
  }

  friend ReplicatedShare<R> operator*(ReplicatedShare<R> lhs, const R& rhs) {
    lhs *= rhs;
    return lhs;
  }

  ReplicatedShare<R>& add(R val, int pid) {
    if (pid == 0) {
      values_[0] += val;
    } else if (pid == 1) {
      values_[2] += val;
    }

    return *this;
  }
};

// Contains all elements of a secret sharing. Used only for generating dummy
// preprocessing data.
template <class R>
struct DummyShare {
  std::array<R, 6> share_elements;

  DummyShare() = default;

  explicit DummyShare(std::array<R, 6> share_elements)
      : share_elements(std::move(share_elements)) {}

  DummyShare(R secret, emp::PRG& prg) {
    prg.random_data(share_elements.data(), sizeof(R) * 5);

    R sum = share_elements[0];
    for (int i = 1; i < 5; ++i) {
      sum += share_elements[i];
    }
    share_elements[5] = secret - sum;
  }

  void randomize(emp::PRG& prg) {
    prg.random_data(share_elements.data(), sizeof(R) * 6);
  }

  [[nodiscard]] R secret() const {
    R sum = share_elements[0];
    for (size_t i = 1; i < 6; ++i) {
      sum += share_elements[i];
    }

    return sum;
  }

  DummyShare<R>& operator+=(const DummyShare<R>& rhs) {
    for (size_t i = 0; i < 6; ++i) {
      share_elements[i] += rhs.share_elements[i];
    }

    return *this;
  }

  friend DummyShare<R> operator+(DummyShare<R> lhs, const DummyShare<R>& rhs) {
    lhs += rhs;
    return lhs;
  }

  DummyShare<R>& operator-=(const DummyShare<R>& rhs) {
    for (size_t i = 0; i < 6; ++i) {
      share_elements[i] -= rhs.share_elements[i];
    }

    return *this;
  }

  friend DummyShare<R> operator-(DummyShare<R> lhs, const DummyShare<R>& rhs) {
    lhs -= rhs;
    return lhs;
  }

  DummyShare<R>& operator*=(const R& rhs) {
    for (size_t i = 0; i < 6; ++i) {
      share_elements[i] *= rhs;
    }

    return *this;
  }

  friend DummyShare<R> operator*(DummyShare<R> lhs, const R& rhs) {
    lhs *= rhs;
    return lhs;
  }

  friend DummyShare<R> operator*(const R& lhs, DummyShare<R> rhs) {
    // Assumes abelian ring.
    rhs *= lhs;
    return rhs;
  }

  ReplicatedShare<R> getRSS(size_t pid) {
    return ReplicatedShare<R>({getShareElement(pid, pidFromOffset(pid, 1)),
                               getShareElement(pid, pidFromOffset(pid, 2)),
                               getShareElement(pid, pidFromOffset(pid, 3))});
  }

  R getShareElement(size_t i, size_t j) {
    return share_elements.at(upperTriangularToArray(i, j));
  }
};

template <>
void ReplicatedShare<BoolRing>::randomize(emp::PRG& prg);

template <>
DummyShare<BoolRing>::DummyShare(BoolRing secret, emp::PRG& prg);

template <>
void DummyShare<BoolRing>::randomize(emp::PRG& prg);
};  // namespace quadsquad
