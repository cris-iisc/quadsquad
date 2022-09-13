#include "jump_provider.h"

#include <algorithm>
#include <boost/format.hpp>
#include <future>
#include <stdexcept>

#include "helpers.h"

namespace quadsquad {
JumpProvider::JumpProvider(int my_id) : id_(my_id), recv_lengths_{}, send_{} {}

void JumpProvider::reset() {
  for (size_t i = 0; i < 4; ++i) {
    for (size_t j = 0; j < 4; ++j) {
      send_[i][j] = false;
      send_hash_[i][j].reset();
      send_values_[i][j].clear();
      recv_lengths_[i][j] = 0;
      recv_values_[i][j].clear();
    }
  }
}

bool JumpProvider::isHashSender(int sender, int other_sender, int receiver) {
  if (other_sender == pidFromOffset(sender, 1)) {
    return receiver == pidFromOffset(sender, 3);
  }
  return receiver == pidFromOffset(sender, 1);
}

void JumpProvider::jumpUpdate(int sender1, int sender2, int receiver,
                              size_t nbytes, const void* data) {
  if (sender1 == sender2 || sender1 == receiver || sender2 == receiver) {
    throw std::invalid_argument(boost::str(
        boost::format(
            "ID, other_sender and receiver must be distinct for Jump3 but "
            "received sender1=%1%, sender2=%2% and receiver=%3%") %
        sender1 % sender2 % receiver));
  }

  if (id_ == receiver) {
    recv_lengths_[std::min(sender1, sender2)][std::max(sender1, sender2)] +=
        nbytes;
    return;
  }
  if (id_ != sender1 && id_ != sender2) {
    return;
  }

  auto other_sender = id_ ^ sender1 ^ sender2;
  if (isHashSender(id_, other_sender, receiver)) {
    send_hash_[other_sender][receiver].put(data, nbytes);
  } else {
    const auto* temp = static_cast<const uint8_t*>(data);
    auto& values = send_values_[other_sender][receiver];
    values.insert(values.end(), temp, temp + nbytes);
  }
  send_[other_sender][receiver] = true;
}

void JumpProvider::communicate(io::NetIOMP<4>& network, ThreadPool& tpool) {
  std::vector<std::future<void>> res;

  // Send data.
  for (int receiver = 0; receiver < 4; ++receiver) {
    if (receiver == id_) {
      continue;
    }

    res.push_back(tpool.enqueue([&, receiver]() {
      for (int other_sender = 0; other_sender < 4; ++other_sender) {
        if (other_sender == receiver || other_sender == id_) {
          continue;
        }

        auto& hash = send_hash_[other_sender][receiver];
        auto& values = send_values_[other_sender][receiver];
        bool should_send = send_[other_sender][receiver];
        if (should_send && isHashSender(id_, other_sender, receiver)) {
          std::array<char, emp::Hash::DIGEST_SIZE> digest{};
          hash.digest(digest.data());
          network.send(receiver, digest.data(), digest.size());
        } else if (should_send) {
          std::array<char, emp::Hash::DIGEST_SIZE> digest{};
          network.send(receiver, values.data(), values.size());
        }
      }

      network.flush(receiver);
    }));
  }

  // Receive data.
  std::array<std::array<std::array<char, emp::Hash::DIGEST_SIZE>, 4>, 4>
      recv_hash{};
  for (int sender = 0; sender < 4; ++sender) {
    if (sender == id_) {
      continue;
    }

    res.push_back(tpool.enqueue([&, sender]() {
      for (int other_sender = 0; other_sender < 4; ++other_sender) {
        if (other_sender == sender || other_sender == id_) {
          continue;
        }

        auto min_sender = std::min(sender, other_sender);
        auto max_sender = std::max(sender, other_sender);
        auto nbytes = recv_lengths_[min_sender][max_sender];

        if (isHashSender(sender, other_sender, id_) && (nbytes != 0)) {
          network.recv(sender, recv_hash[min_sender][max_sender].data(),
                       emp::Hash::DIGEST_SIZE);
        } else if (nbytes != 0) {
          auto& values = recv_values_[min_sender][max_sender];
          values.resize(values.size() + nbytes);
          network.recv(sender, values.data() + values.size() - nbytes, nbytes);
        }
      }
    }));
  }

  for (auto& f : res) {
    f.get();
  }

  // Verify.
  emp::Hash hash;
  std::array<char, emp::Hash::DIGEST_SIZE> digest{};
  for (int sender1 = 0; sender1 < 4; ++sender1) {
    for (int sender2 = sender1 + 1; sender2 < 4; ++sender2) {
      if (sender1 == id_ || sender2 == id_) {
        continue;
      }

      auto& values = recv_values_[sender1][sender2];
      if (!values.empty()) {
        hash.put(values.data(), values.size());
        hash.digest(digest.data());
        if (digest != recv_hash[sender1][sender2]) {
          throw std::runtime_error(boost::str(
              boost::format("Malicious behaviour detected in jump with "
                            "sender1=%1%, sender2=%2% and receiver=%3%") %
              sender1 % sender2 % id_));
        }
        hash.reset();
      }
    }
  }
}

const std::vector<uint8_t>& JumpProvider::getValues(int sender1, int sender2) {
  return recv_values_[std::min(sender1, sender2)][std::max(sender1, sender2)];
}
};  // namespace quadsquad
