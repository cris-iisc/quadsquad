#pragma once

#include <emp-ot/emp-ot.h>
#include <emp-tool/emp-tool.h>

#include <array>
#include <vector>

#include "types.h"

namespace quadsquad {

class OTProvider {
  std::array<NetIO*, 1> send_ios_;
  std::array<NetIO*, 1> recv_ios_;
  std::unique_ptr<emp::FerretCOT<emp::NetIO>> send_ot_;
  std::unique_ptr<emp::FerretCOT<emp::NetIO>> recv_ot_;

 public:
  OTProvider(int my_id, int other_id, emp::NetIO* send_io, emp::NetIO* recv_io);

  std::vector<Ring> multiplySend(const std::vector<Ring>& inputs,
                                 emp::PRG& prg);
  std::vector<Ring> multiplyReceive(const std::vector<Ring>& inputs);

  void send(const uint64_t* data0, const uint64_t* data1, int64_t length);
  void recv(uint64_t* data, const bool* r, int64_t length);
};
};  // namespace quadsquad
