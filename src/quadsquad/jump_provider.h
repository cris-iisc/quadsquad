#pragma once

#include <emp-tool/emp-tool.h>

#include <array>

#include "../io/netmp.h"

namespace quadsquad {

// Manages instances of jump.
class JumpProvider {
  int id_;
  std::array<std::array<bool, 4>, 4> send_;
  std::array<std::array<emp::Hash, 4>, 4> send_hash_;
  std::array<std::array<std::vector<uint8_t>, 4>, 4> send_values_;
  std::array<std::array<size_t, 4>, 4> recv_lengths_;
  std::array<std::array<std::vector<uint8_t>, 4>, 4> recv_values_;

  static bool isHashSender(int sender, int other_sender, int receiver);

 public:
  explicit JumpProvider(int my_id);

  void reset();

  void jumpUpdate(int sender1, int sender2, int receiver, size_t nbytes,
                  const void* data = nullptr);
  void communicate(io::NetIOMP<4>& network, ThreadPool& tpool);
  const std::vector<uint8_t>& getValues(int sender1, int sender2);
};

};  // namespace quadsquad
