#define BOOST_TEST_MODULE jump
#include <emp-tool/emp-tool.h>
#include <io/netmp.h>
#include <quadsquad/jump_provider.h>

#include <boost/test/included/unit_test.hpp>
#include <future>
#include <string>
#include <vector>

using namespace quadsquad;

BOOST_AUTO_TEST_SUITE(jump_provider)

BOOST_AUTO_TEST_CASE(receiver_offset_2) {
  std::string message("A test string.");
  std::vector<uint8_t> input(message.begin(), message.end());

  std::vector<std::future<void>> parties;
  for (int i = 0; i < 4; ++i) {
    parties.push_back(std::async(std::launch::async, [&, i]() {
      io::NetIOMP<4> network(i, 10000, nullptr, true);
      JumpProvider jump(i);
      ThreadPool tpool(1);

      jump.jumpUpdate(1, 2, 3, input.size(), input.data());
      jump.communicate(network, tpool);

      if (i == 3) {
        BOOST_TEST(jump.getValues(1, 2) == input);
      }
    }));
  }

  for (auto& p : parties) {
    p.wait();
  }
}

BOOST_AUTO_TEST_CASE(receiver_offset_3) {
  std::string message("A test string.");
  std::vector<uint8_t> input(message.begin(), message.end());

  std::vector<std::future<void>> parties;
  for (int i = 0; i < 4; ++i) {
    parties.push_back(std::async(std::launch::async, [&, i]() {
      io::NetIOMP<4> network(i, 10000, nullptr, true);
      JumpProvider jump(i);
      ThreadPool tpool(1);

      jump.jumpUpdate(0, 1, 3, input.size(), input.data());
      jump.communicate(network, tpool);

      if (i == 3) {
        BOOST_TEST(jump.getValues(0, 1) == input);
      }
    }));
  }

  for (auto& p : parties) {
    p.wait();
  }
}

BOOST_AUTO_TEST_CASE(all_combinations) {
  std::string message("A test string.");
  std::vector<uint8_t> input(message.begin(), message.end());

  std::vector<std::future<void>> parties;
  for (int i = 0; i < 4; ++i) {
    parties.push_back(std::async(std::launch::async, [&, i]() {
      io::NetIOMP<4> network(i, 10000, nullptr, true);
      JumpProvider jump(i);
      ThreadPool tpool(1);

      for (int sender1 = 0; sender1 < 4; ++sender1) {
        for (int sender2 = sender1 + 1; sender2 < 4; ++sender2) {
          for (int receiver = 0; receiver < 4; ++receiver) {
            if (receiver == sender1 || receiver == sender2) {
              continue;
            }
            jump.jumpUpdate(sender1, sender2, receiver, input.size(),
                            input.data());
          }
        }
      }

      jump.communicate(network, tpool);

      for (int sender1 = 0; sender1 < 4; ++sender1) {
        for (int sender2 = sender1 + 1; sender2 < 4; ++sender2) {
          if (i == sender1 || i == sender2) {
            continue;
          }
          BOOST_TEST(jump.getValues(sender1, sender2) == input);
        }
      }
    }));
  }

  for (auto& p : parties) {
    p.wait();
  }
}

BOOST_AUTO_TEST_SUITE_END()
