#include "rand_gen_pool.h"

#include <algorithm>

#include "helpers.h"

namespace quadsquad {

RandGenPool::RandGenPool(int my_id, uint64_t seed) : id_{my_id} {
  auto seed_block = emp::makeBlock(seed, 0);
  auto encode = [](std::vector<int> parties) {
    // Ordering of parties should not impact encoding.
    std::sort(parties.begin(), parties.end());

    int res = 0;
    int pow = 1;
    for (int party : parties) {
      res += pow * party;
      pow *= 4;
    }

    return res;
  };

  for (int i = 0; i < 4; ++i) {
    v_rgen_.emplace_back(&seed_block, encode({id_, i}));
  }

  for (int i = 0; i < 4; ++i) {
    std::vector<int> parties = {id_};
    for (int j = 0; j < 4; ++j) {
      if (j == i || j == id_) {
        continue;
      }

      parties.push_back(j);
    }
    v_rgen_.emplace_back(&seed_block, encode(parties));
  }
}

emp::PRG& RandGenPool::self() { return v_rgen_[id_]; }

emp::PRG& RandGenPool::all() { return v_rgen_[id_ + 4]; }

emp::PRG& RandGenPool::get(int pid) { return v_rgen_.at(pid); }

emp::PRG& RandGenPool::getRelative(int offset) {
  return v_rgen_.at(pidFromOffset(id_, offset));
}

emp::PRG& RandGenPool::get(int pid1, int pid2) {
  return v_rgen_.at(4 + (pid1 ^ pid2 ^ id_));
}

emp::PRG& RandGenPool::getRelative(int offset1, int offset2) {
  return get(pidFromOffset(id_, offset1), pidFromOffset(id_, offset2));
}

emp::PRG& RandGenPool::getComplement(int pid) { return v_rgen_.at(4 + pid); }

emp::PRG& RandGenPool::getComplementRelative(int offset) {
  return getComplement(pidFromOffset(id_, offset));
}

}  // namespace quadsquad
