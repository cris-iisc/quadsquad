// MIT License
//
// Copyright (c) 2018 Xiao Wang (wangxiao@gmail.com)
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

// The following code has been adopted from
// https://github.com/emp-toolkit/emp-agmpc. It has been modified to define the
// class within a namespace and add additional methods (sendRelative,
// recvRelative).

#pragma once

#include <emp-tool/emp-tool.h>

namespace io {
using namespace emp;

template <int nP>
class NetIOMP {
 public:
  NetIO* ios[nP];
  NetIO* ios2[nP];
  int party;
  bool sent[nP];

  NetIOMP(int party, int port, char* IP[], bool localhost = false) {
    this->party = party;
    memset(sent, false, nP);
    for (int i = 0; i < nP; ++i) {
      for (int j = i + 1; j < nP; ++j) {
        if (i == party) {
          usleep(1000);
          if (localhost) {
            ios[j] = new NetIO("127.0.0.1", port + 2 * (i * nP + j), true);
          } else {
            ios[j] = new NetIO(IP[j], port + 2 * (i), true);
          }
          ios[j]->set_nodelay();

          usleep(1000);
          if (localhost) {
            ios2[j] = new NetIO(nullptr, port + 2 * (i * nP + j) + 1, true);
          } else {
            ios2[j] = new NetIO(nullptr, port + 2 * (j) + 1, true);
          }
          ios2[j]->set_nodelay();
        } else if (j == party) {
          usleep(1000);
          if (localhost) {
            ios[i] = new NetIO(nullptr, port + 2 * (i * nP + j), true);
          } else {
            ios[i] = new NetIO(nullptr, port + 2 * (i), true);
          }
          ios[i]->set_nodelay();

          usleep(1000);
          if (localhost) {
            ios2[i] = new NetIO("127.0.0.1", port + 2 * (i * nP + j) + 1, true);
          } else {
            ios2[i] = new NetIO(IP[i], port + 2 * (j) + 1, true);
          }
        }
      }
    }
  }

  int64_t count() {
    int64_t res = 0;
    for (int i = 0; i < nP; ++i)
      if (i != party) {
        res += ios[i]->counter;
        res += ios2[i]->counter;
      }
    return res;
  }

  void resetStats() {
    for (int i = 0; i < nP; ++i) {
      if (i != party) {
        ios[i]->counter = 0;
        ios2[i]->counter = 0;
      }
    }
  }

  ~NetIOMP() {
    for (int i = 0; i < nP; ++i)
      if (i != party) {
        delete ios[i];
        delete ios2[i];
      }
  }

  void send(int dst, const void* data, size_t len) {
    if (dst != -1 and dst != party) {
      if (party < dst)
        ios[dst]->send_data(data, len);
      else
        ios2[dst]->send_data(data, len);
      sent[dst] = true;
    }
#ifdef __clang__
    flush(dst);
#endif
  }

  void sendRelative(int offset, const void* data, size_t len) {
    int dst = (party + offset) % nP;
    if (dst < 0) {
      dst += nP;
    }
    send(dst, data, len);
  }

  void sendBool(int dst, const bool* data, size_t len) {
    for (int i = 0; i < len;) {
      uint64_t tmp = 0;
      for (int j = 0; j < 64 && i < len; ++i, ++j) {
        if (data[i]) {
          tmp |= (0x1ULL << j);
        }
      }
      send(dst, &tmp, 8);
    }
  }

  void sendBoolRelative(int offset, const bool* data, size_t len) {
    int dst = (party + offset) % nP;
    if (dst < 0) {
      dst += nP;
    }
    sendBool(dst, data, len);
  }

  void recv(int src, void* data, size_t len) {
    if (src != -1 && src != party) {
      if (sent[src]) flush(src);
      if (src < party)
        ios[src]->recv_data(data, len);
      else
        ios2[src]->recv_data(data, len);
    }
  }

  void recvRelative(int offset, void* data, size_t len) {
    int src = (party + offset) % nP;
    if (src < 0) {
      src += nP;
    }
    recv(src, data, len);
  }

  void recvBool(int src, bool* data, size_t len) {
    for (int i = 0; i < len;) {
      uint64_t tmp = 0;
      recv(src, &tmp, 8);
      for (int j = 63; j >= 0 && i < len; ++i, --j) {
        data[i] = (tmp & 0x1) == 0x1;
        tmp >>= 1;
      }
    }
  }

  void recvRelative(int offset, bool* data, size_t len) {
    int src = (party + offset) % nP;
    if (src < 0) {
      src += nP;
    }
    recvBool(src, data, len);
  }

  NetIO*& get(size_t idx, bool b = false) {
    if (b)
      return ios[idx];
    else
      return ios2[idx];
  }

  NetIO*& getSendChannel(size_t idx) {
    if (party < idx) {
      return ios[idx];
    }

    return ios2[idx];
  }

  NetIO*& getRecvChannel(size_t idx) {
    if (idx < party) {
      return ios[idx];
    }

    return ios2[idx];
  }

  void flush(int idx = -1) {
    if (idx == -1) {
      for (int i = 0; i < nP; ++i) {
        if (i != party) {
          ios[i]->flush();
          ios2[i]->flush();
        }
      }
    } else {
      if (party < idx) {
        ios[idx]->flush();
      } else {
        ios2[idx]->flush();
      }
    }
  }

  void sync() {
    for (int i = 0; i < nP; ++i) {
      for (int j = 0; j < nP; ++j) {
        if (i < j) {
          if (i == party) {
            ios[j]->sync();
            ios2[j]->sync();
          } else if (j == party) {
            ios[i]->sync();
            ios2[i]->sync();
          }
        }
      }
    }
  }
};
};  // namespace io
