#define BOOST_TEST_MODULE utils
#include <emp-tool/emp-tool.h>
#include <quadsquad/rand_gen_pool.h>
#include <utils/circuit.h>
#include <utils/liquidity_matching.h>
#include <utils/neural_network.h>

#include <boost/test/data/monomorphic.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/included/unit_test.hpp>
#include <map>
#include <random>

using namespace quadsquad;
using namespace quadsquad::utils;
namespace bdata = boost::unit_test::data;

constexpr int TEST_DATA_MAX_VAL = 1000;

BOOST_AUTO_TEST_SUITE(circuit)

BOOST_DATA_TEST_CASE(no_op_circuit,
                     bdata::random(0, TEST_DATA_MAX_VAL) ^ bdata::xrange(1),
                     input, idx) {
  Circuit<int> circ;
  auto wa = circ.newInputWire();
  circ.setAsOutput(wa);

  auto output = circ.evaluate({{wa, input}});

  BOOST_TEST(output[0] == input);

  auto level_circ = circ.orderGatesByLevel();
  BOOST_TEST(level_circ.gates_by_level.size() == 1);
  BOOST_TEST(level_circ.num_gates == 1);
  BOOST_TEST(level_circ.count[GateType::kInp] == 1);
  BOOST_TEST(level_circ.count[GateType::kAdd] == 0);
  BOOST_TEST(level_circ.count[GateType::kMul] == 0);
  BOOST_TEST(level_circ.count[GateType::kSub] == 0);
  BOOST_TEST(level_circ.count[GateType::kConstAdd] == 0);
  BOOST_TEST(level_circ.count[GateType::kConstMul] == 0);
  BOOST_TEST(level_circ.count[GateType::kDotprod] == 0);
  BOOST_TEST(level_circ.count[GateType::kTrdotp] == 0);
  BOOST_TEST(level_circ.count[GateType::kRelu] == 0);
}

BOOST_DATA_TEST_CASE(add_gate,
                     bdata::random(0, TEST_DATA_MAX_VAL) ^
                         bdata::random(0, TEST_DATA_MAX_VAL) ^ bdata::xrange(1),
                     input_a, input_b, idx) {
  Circuit<int> circ;
  auto wa = circ.newInputWire();
  auto wb = circ.newInputWire();
  auto wsum = circ.addGate(GateType::kAdd, wa, wb);
  circ.setAsOutput(wsum);

  auto output = circ.evaluate({{wa, input_a}, {wb, input_b}});

  BOOST_TEST(output[0] == input_a + input_b);

  auto level_circ = circ.orderGatesByLevel();
  BOOST_TEST(level_circ.gates_by_level.size() == 1);
  BOOST_TEST(level_circ.num_gates == 3);
  BOOST_TEST(level_circ.count[GateType::kInp] == 2);
  BOOST_TEST(level_circ.count[GateType::kAdd] == 1);
  BOOST_TEST(level_circ.count[GateType::kMul] == 0);
  BOOST_TEST(level_circ.count[GateType::kSub] == 0);
  BOOST_TEST(level_circ.count[GateType::kConstAdd] == 0);
  BOOST_TEST(level_circ.count[GateType::kConstMul] == 0);
  BOOST_TEST(level_circ.count[GateType::kDotprod] == 0);
  BOOST_TEST(level_circ.count[GateType::kTrdotp] == 0);
  BOOST_TEST(level_circ.count[GateType::kRelu] == 0);
}

BOOST_DATA_TEST_CASE(sub_gate,
                     bdata::random(0, TEST_DATA_MAX_VAL) ^
                         bdata::random(0, TEST_DATA_MAX_VAL) ^ bdata::xrange(1),
                     input_a, input_b, idx) {
  Circuit<int> circ;
  auto wa = circ.newInputWire();
  auto wb = circ.newInputWire();
  auto wdiff = circ.addGate(GateType::kSub, wa, wb);
  circ.setAsOutput(wdiff);

  auto output = circ.evaluate({{wa, input_a}, {wb, input_b}});

  BOOST_TEST(output[0] == input_a - input_b);

  auto level_circ = circ.orderGatesByLevel();
  BOOST_TEST(level_circ.gates_by_level.size() == 1);
  BOOST_TEST(level_circ.num_gates == 3);
  BOOST_TEST(level_circ.count[GateType::kInp] == 2);
  BOOST_TEST(level_circ.count[GateType::kAdd] == 0);
  BOOST_TEST(level_circ.count[GateType::kMul] == 0);
  BOOST_TEST(level_circ.count[GateType::kSub] == 1);
  BOOST_TEST(level_circ.count[GateType::kConstAdd] == 0);
  BOOST_TEST(level_circ.count[GateType::kConstMul] == 0);
  BOOST_TEST(level_circ.count[GateType::kDotprod] == 0);
  BOOST_TEST(level_circ.count[GateType::kTrdotp] == 0);
  BOOST_TEST(level_circ.count[GateType::kRelu] == 0);
}

BOOST_DATA_TEST_CASE(mul_gate,
                     bdata::random(0, TEST_DATA_MAX_VAL) ^
                         bdata::random(0, TEST_DATA_MAX_VAL) ^ bdata::xrange(1),
                     input_a, input_b, idx) {
  Circuit<int> circ;
  auto wa = circ.newInputWire();
  auto wb = circ.newInputWire();
  auto wprod = circ.addGate(GateType::kMul, wa, wb);
  circ.setAsOutput(wprod);

  auto output = circ.evaluate({{wa, input_a}, {wb, input_b}});

  BOOST_TEST(output[0] == input_a * input_b);

  auto level_circ = circ.orderGatesByLevel();
  BOOST_TEST(level_circ.gates_by_level.size() == 2);
  BOOST_TEST(level_circ.num_gates == 3);
  BOOST_TEST(level_circ.count[GateType::kInp] == 2);
  BOOST_TEST(level_circ.count[GateType::kAdd] == 0);
  BOOST_TEST(level_circ.count[GateType::kMul] == 1);
  BOOST_TEST(level_circ.count[GateType::kSub] == 0);
  BOOST_TEST(level_circ.count[GateType::kConstAdd] == 0);
  BOOST_TEST(level_circ.count[GateType::kConstMul] == 0);
  BOOST_TEST(level_circ.count[GateType::kDotprod] == 0);
  BOOST_TEST(level_circ.count[GateType::kTrdotp] == 0);
  BOOST_TEST(level_circ.count[GateType::kRelu] == 0);
}

BOOST_DATA_TEST_CASE(const_add_gate,
                     bdata::random(0, TEST_DATA_MAX_VAL) ^
                         bdata::random(0, TEST_DATA_MAX_VAL) ^ bdata::xrange(1),
                     input_a, input_b, idx) {
  Circuit<int> circ;
  auto wa = circ.newInputWire();
  auto wsum = circ.addConstOpGate(GateType::kConstAdd, wa, input_b);
  circ.setAsOutput(wsum);

  auto output = circ.evaluate({{wa, input_a}});

  BOOST_TEST(output[0] == input_a + input_b);

  auto level_circ = circ.orderGatesByLevel();
  BOOST_TEST(level_circ.gates_by_level.size() == 1);
  BOOST_TEST(level_circ.num_gates == 2);
  BOOST_TEST(level_circ.count[GateType::kInp] == 1);
  BOOST_TEST(level_circ.count[GateType::kAdd] == 0);
  BOOST_TEST(level_circ.count[GateType::kMul] == 0);
  BOOST_TEST(level_circ.count[GateType::kSub] == 0);
  BOOST_TEST(level_circ.count[GateType::kConstAdd] == 1);
  BOOST_TEST(level_circ.count[GateType::kConstMul] == 0);
  BOOST_TEST(level_circ.count[GateType::kDotprod] == 0);
  BOOST_TEST(level_circ.count[GateType::kTrdotp] == 0);
  BOOST_TEST(level_circ.count[GateType::kRelu] == 0);
}

BOOST_DATA_TEST_CASE(const_mul_gate,
                     bdata::random(0, TEST_DATA_MAX_VAL) ^
                         bdata::random(0, TEST_DATA_MAX_VAL) ^ bdata::xrange(1),
                     input_a, input_b, idx) {
  Circuit<int> circ;
  auto wa = circ.newInputWire();
  auto wprod = circ.addConstOpGate(GateType::kConstMul, wa, input_b);
  circ.setAsOutput(wprod);

  auto output = circ.evaluate({{wa, input_a}});

  BOOST_TEST(output[0] == input_a * input_b);

  auto level_circ = circ.orderGatesByLevel();
  BOOST_TEST(level_circ.gates_by_level.size() == 1);
  BOOST_TEST(level_circ.num_gates == 2);
  BOOST_TEST(level_circ.count[GateType::kInp] == 1);
  BOOST_TEST(level_circ.count[GateType::kAdd] == 0);
  BOOST_TEST(level_circ.count[GateType::kMul] == 0);
  BOOST_TEST(level_circ.count[GateType::kSub] == 0);
  BOOST_TEST(level_circ.count[GateType::kConstAdd] == 0);
  BOOST_TEST(level_circ.count[GateType::kConstMul] == 1);
  BOOST_TEST(level_circ.count[GateType::kDotprod] == 0);
  BOOST_TEST(level_circ.count[GateType::kTrdotp] == 0);
  BOOST_TEST(level_circ.count[GateType::kRelu] == 0);
}

BOOST_DATA_TEST_CASE(depth_2_circuit,
                     bdata::random(0, TEST_DATA_MAX_VAL) ^
                         bdata::random(0, TEST_DATA_MAX_VAL) ^
                         bdata::random(0, TEST_DATA_MAX_VAL) ^
                         bdata::random(0, TEST_DATA_MAX_VAL) ^ bdata::xrange(1),
                     input_a, input_b, input_c, input_d, idx) {
  std::vector<int> inputs = {input_a, input_b, input_c, input_d};

  Circuit<int> circ;
  std::vector<wire_t> input_wires;
  for (size_t i = 0; i < inputs.size(); ++i) {
    input_wires.push_back(circ.newInputWire());
  }

  auto w_aab = circ.addGate(GateType::kAdd, input_wires[0], input_wires[1]);
  auto v_aab = inputs[0] + inputs[1];

  auto w_cmd = circ.addGate(GateType::kMul, input_wires[2], input_wires[3]);
  auto v_cmd = inputs[2] * inputs[3];

  auto w_mout = circ.addGate(GateType::kMul, w_aab, w_cmd);
  auto v_mout = v_aab * v_cmd;

  auto w_aout = circ.addGate(GateType::kAdd, w_aab, w_cmd);
  auto v_aout = v_aab + v_cmd;

  circ.setAsOutput(w_mout);
  circ.setAsOutput(w_aout);

  std::unordered_map<wire_t, int> input_map;
  for (size_t i = 0; i < inputs.size(); ++i) {
    input_map[input_wires[i]] = inputs[i];
  }
  auto outputs = circ.evaluate(input_map);

  BOOST_TEST(outputs[0] == v_mout);
  BOOST_TEST(outputs[1] == v_aout);

  auto level_circ = circ.orderGatesByLevel();
  BOOST_TEST(level_circ.gates_by_level.size() == 3);
  BOOST_TEST(level_circ.num_gates == 8);
  BOOST_TEST(level_circ.count[GateType::kInp] == 4);
  BOOST_TEST(level_circ.count[GateType::kAdd] == 2);
  BOOST_TEST(level_circ.count[GateType::kMul] == 2);
  BOOST_TEST(level_circ.count[GateType::kSub] == 0);
  BOOST_TEST(level_circ.count[GateType::kConstAdd] == 0);
  BOOST_TEST(level_circ.count[GateType::kConstMul] == 0);
  BOOST_TEST(level_circ.count[GateType::kDotprod] == 0);
  BOOST_TEST(level_circ.count[GateType::kTrdotp] == 0);
  BOOST_TEST(level_circ.count[GateType::kRelu] == 0);
}

BOOST_DATA_TEST_CASE(double_relu_gate,
                     bdata::random(0, TEST_DATA_MAX_VAL) ^
                         bdata::random(0, TEST_DATA_MAX_VAL) ^ bdata::xrange(1),
                     input_a, input_b, idx) {
  Circuit<Ring> circ;
  auto wa = circ.newInputWire();
  auto wb = circ.newInputWire();
  auto wrelu_a = circ.addGate(GateType::kRelu, wa);
  auto wprod = circ.addGate(GateType::kMul, wrelu_a, wb);
  auto wrelu_prod = circ.addGate(GateType::kRelu, wprod);

  circ.setAsOutput(wrelu_a);
  circ.setAsOutput(wrelu_prod);

  auto output = circ.evaluate({{wa, input_a}, {wb, input_b}});
  BOOST_TEST(output[0] == input_a);
  BOOST_TEST(output[1] == input_a * input_b);
  output = circ.evaluate({{wa, -1 * input_a}, {wb, input_b}});
  BOOST_TEST(output[0] == 0);
  BOOST_TEST(output[1] == 0);

  auto level_circ = circ.orderGatesByLevel();
  BOOST_TEST(level_circ.gates_by_level.size() == 4);
  BOOST_TEST(level_circ.num_gates == 5);
  BOOST_TEST(level_circ.count[GateType::kInp] == 2);
  BOOST_TEST(level_circ.count[GateType::kAdd] == 0);
  BOOST_TEST(level_circ.count[GateType::kMul] == 1);
  BOOST_TEST(level_circ.count[GateType::kSub] == 0);
  BOOST_TEST(level_circ.count[GateType::kConstAdd] == 0);
  BOOST_TEST(level_circ.count[GateType::kConstMul] == 0);
  BOOST_TEST(level_circ.count[GateType::kDotprod] == 0);
  BOOST_TEST(level_circ.count[GateType::kTrdotp] == 0);
  BOOST_TEST(level_circ.count[GateType::kRelu] == 2);
}

BOOST_AUTO_TEST_CASE(dotp_gate) {
  int nf = 10;

  Circuit<int> circ;
  std::vector<wire_t> vwa(nf), vwb(nf);
  std::vector<int> vx(nf), vy(nf);
  int result = 0;

  std::mt19937 gen(200);
  std::uniform_int_distribution<int> distrib(0, TEST_DATA_MAX_VAL);
  for (int i = 0; i < nf; i++) {
    vwa[i] = circ.newInputWire();
    vwb[i] = circ.newInputWire();
    vx[i] = distrib(gen);
    vy[i] = distrib(gen);
    result += vx[i] * vy[i];
  }

  auto wdotp = circ.addGate(GateType::kDotprod, vwa, vwb);
  circ.setAsOutput(wdotp);

  std::unordered_map<wire_t, int> input_map;
  for (size_t i = 0; i < nf; ++i) {
    input_map[vwa[i]] = vx[i];
    input_map[vwb[i]] = vy[i];
  }
  auto output = circ.evaluate(input_map);
  BOOST_TEST(output[0] == result);

  auto level_circ = circ.orderGatesByLevel();
  BOOST_TEST(level_circ.gates_by_level.size() == 2);
  BOOST_TEST(level_circ.num_gates == 2 * nf + 1);
  BOOST_TEST(level_circ.count[GateType::kInp] == 2 * nf);
  BOOST_TEST(level_circ.count[GateType::kAdd] == 0);
  BOOST_TEST(level_circ.count[GateType::kMul] == 0);
  BOOST_TEST(level_circ.count[GateType::kSub] == 0);
  BOOST_TEST(level_circ.count[GateType::kConstAdd] == 0);
  BOOST_TEST(level_circ.count[GateType::kConstMul] == 0);
  BOOST_TEST(level_circ.count[GateType::kDotprod] == 1);
  BOOST_TEST(level_circ.count[GateType::kTrdotp] == 0);
  BOOST_TEST(level_circ.count[GateType::kRelu] == 0);
}

BOOST_AUTO_TEST_CASE(msb_circuit) {
  Circuit circ = Circuit<BoolRing>::generatePPAMSB();
  uint64_t a, b;
  a = 4;
  b = -5;
  auto aval = static_cast<uint64_t>(a);
  auto bval = static_cast<uint64_t>(b);
  auto sum = aval + bval;

  std::vector<BoolRing> aval_bits(64);
  std::vector<BoolRing> bval_bits(64);
  std::vector<BoolRing> sum_bits(64);
  for (size_t i = 0; i < 64; ++i) {
    aval_bits[i] = ((aval >> i) & 1) == 1;
    bval_bits[i] = ((bval >> i) & 1) == 1;
    sum_bits[i] = ((sum >> i) & 1) == 1;
  }

  std::unordered_map<wire_t, BoolRing> input_map;
  for (size_t i = 0; i < 64; ++i) {
    input_map[i] = aval_bits[i];
    input_map[64 + i] = bval_bits[i];
  }

  auto output = circ.evaluate(input_map);

  BOOST_TEST(output[0] == sum_bits[63]);
}

BOOST_DATA_TEST_CASE(ppa_circuit,
                     bdata::random(0, TEST_DATA_MAX_VAL) ^
                         bdata::random(0, TEST_DATA_MAX_VAL) ^ bdata::xrange(1),
                     input_a, input_b, idx) {
  Circuit circ = Circuit<BoolRing>::generatePPA();
  uint64_t a, b;
  a = 12;
  b = -20;
  auto aval = static_cast<uint64_t>(a);
  auto bval = static_cast<uint64_t>(b);
  auto sum = aval + bval;
  std::vector<BoolRing> aval_bits(64);
  std::vector<BoolRing> bval_bits(64);
  std::vector<BoolRing> sum_bits(64);
  for (size_t i = 0; i < 64; ++i) {
    aval_bits[i] = ((aval >> i) & 1) == 1;
    bval_bits[i] = ((bval >> i) & 1) == 1;
    sum_bits[i] = ((sum >> i) & 1) == 1;
  }
  std::unordered_map<wire_t, BoolRing> input_map;
  for (size_t i = 0; i < 64; ++i) {
    input_map[i] = aval_bits[i];
    input_map[64 + i] = bval_bits[i];
  }

  auto output = circ.evaluate(input_map);
  for (int i = 0; i < 64; ++i) {
  }
  BOOST_TEST(output == sum_bits);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(rand_gen_pool)

BOOST_AUTO_TEST_CASE(matching_output) {
  const int num_tests = 10;
  const uint64_t seed = 200;

  for (int i = 0; i < 4; ++i) {
    for (int j = i + 1; j < 4; ++j) {
      auto rpool_i = RandGenPool(i, seed);
      auto rpool_j = RandGenPool(j, seed);

      uint64_t bi;
      uint64_t bj;
      for (int t = 0; t < num_tests; ++t) {
        rpool_i.all().random_data(&bi, sizeof(uint64_t));
        rpool_j.all().random_data(&bj, sizeof(uint64_t));
        BOOST_TEST(bi == bj);

        rpool_i.get(j).random_data(&bi, sizeof(uint64_t));
        rpool_j.get(i).random_data(&bj, sizeof(uint64_t));
        BOOST_TEST(bi == bj);

        for (int k = 0; k < 4; ++k) {
          if (k != i && k != j) {
            emp::block bk;
            rpool_i.get(k, j).random_data(&bi, sizeof(uint64_t));
            rpool_j.get(k, i).random_data(&bj, sizeof(uint64_t));
            BOOST_TEST(bi == bj);
          }
        }
      }
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(bool_ring)

BOOST_AUTO_TEST_CASE(operations) {
  for (int i = 0; i < 1; ++i) {
    for (int j = 0; j < 1; ++j) {
      BoolRing ar(i);
      BoolRing br(j);
      bool a = (i != 0);
      bool b = (j != 0);

      BOOST_TEST((ar + br).val() == (a != b));
      BOOST_TEST((ar - br).val() == (a != b));
      BOOST_TEST((ar * br).val() == (a && b));
    }
  }

  BoolRing a = 0;
  BoolRing b = 1;

  BOOST_TEST(!a.val());
  BOOST_TEST(b.val());
}

BOOST_AUTO_TEST_CASE(add_gate) {
  Circuit<BoolRing> circ;
  auto wa = circ.newInputWire();
  auto wb = circ.newInputWire();
  auto wsum = circ.addGate(GateType::kAdd, wa, wb);
  circ.setAsOutput(wsum);

  for (int i = 0; i < 1; ++i) {
    for (int j = 0; j < 1; ++j) {
      BoolRing a(i);
      BoolRing b(j);
      auto output = circ.evaluate({{wa, a}, {wb, b}});
      BOOST_TEST(output[0] == a * b);
    }
  }
}

BOOST_AUTO_TEST_CASE(mul_gate) {
  Circuit<BoolRing> circ;
  auto wa = circ.newInputWire();
  auto wb = circ.newInputWire();
  auto wprod = circ.addGate(GateType::kMul, wa, wb);
  circ.setAsOutput(wprod);

  for (int i = 0; i < 1; ++i) {
    for (int j = 0; j < 1; ++j) {
      BoolRing a(i);
      BoolRing b(j);
      auto output = circ.evaluate({{wa, a}, {wb, b}});
      BOOST_TEST(output[0] == a * b);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(neural_network)

BOOST_AUTO_TEST_CASE(linear) {
  const size_t batch_size = 2;
  const size_t inp_len = 5;
  const size_t out_len = 2;

  NeuralNetwork<Ring> nn;
  auto input = nn.newInput<2>({batch_size, inp_len});
  auto output = nn.linear(input, out_len);
  nn.setOutput(output);
  const auto& circ = nn.getCircuit();

  auto inp_shape = input.shape();
  BOOST_TEST(inp_shape[0] == batch_size);
  BOOST_TEST(inp_shape[1] == inp_len);

  std::unordered_map<wire_t, Ring> input_map;
  for (size_t b = 0; b < batch_size; ++b) {
    for (size_t i = 0; i < inp_len; ++i) {
      auto wid = input[b][i];
      input_map[wid] = i + 1;
    }
  }

  const auto& weights = nn.getLinearWeights(0);
  auto wt_shape = weights.shape();
  BOOST_TEST(wt_shape[0] == out_len);
  BOOST_TEST(wt_shape[1] == inp_len + 1);

  for (size_t i = 0; i < wt_shape[0]; ++i) {
    for (size_t j = 0; j < wt_shape[1]; ++j) {
      input_map[weights[i][j]] = i * (wt_shape[1] + 1) + j + 1;
    }
  }

  auto res = circ.evaluate(input_map);
  BOOST_TEST(res.size() == out_len * batch_size);
  for (size_t b = 0; b < batch_size; ++b) {
    BOOST_TEST(res[2 * b] == 6);
    BOOST_TEST(res[2 * b + 1] == 13);
  }
}

BOOST_AUTO_TEST_CASE(relu) {
  const size_t len = 5;

  NeuralNetwork<Ring> nn;
  auto input = nn.newInput<2>({1, len});
  auto output = nn.relu(input);
  nn.setOutput(output);
  const auto& circ = nn.getCircuit();

  auto inp_shape = input.shape();

  std::unordered_map<wire_t, Ring> input_map;
  std::vector<Ring> exp_out(len);
  for (size_t i = 0; i < len; ++i) {
    auto wid = input[0][i];

    if (i % 2 == 0) {
      input_map[wid] = -1 - i;
      exp_out[i] = 0;
    } else {
      input_map[wid] = i + 1;
      exp_out[i] = i + 1;
    }
  }

  auto res = circ.evaluate(input_map);
  BOOST_TEST(res == exp_out);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(liquidity_matching)

BOOST_AUTO_TEST_CASE(netting) {
  SoDoGridLock<Ring> gl(3);
  std::unordered_map<wire_t, Ring> imap;

  wire_t wtxn_amt;
  wtxn_amt = gl.newTransaction(0, 1);
  imap[wtxn_amt] = 25;
  wtxn_amt = gl.newTransaction(1, 2);
  imap[wtxn_amt] = 30;
  wtxn_amt = gl.newTransaction(2, 0);
  imap[wtxn_amt] = 35;

  auto balances = gl.initBalances({20, 20, 20}, imap);
  auto selected = gl.initSelectedSet(imap);
  auto new_selected = gl.updateSelectedTransactions(balances, selected);

  auto circ = gl.getCircuit();
  for (auto w : new_selected) {
    circ.setAsOutput(w);
  }

  auto res = circ.evaluate(imap);
  BOOST_TEST(res[0] == 1);
  BOOST_TEST(res[1] == 0);
  for (size_t i = 2; i < res.size(); i++) {
    BOOST_TEST(res[i] == 1);
  }
}

BOOST_AUTO_TEST_CASE(deadlock) {
  SoDoGridLock<Ring> gl(3);
  std::unordered_map<wire_t, Ring> imap;

  wire_t wtxn_amt;
  wtxn_amt = gl.newTransaction(0, 1);
  imap[wtxn_amt] = 25;
  wtxn_amt = gl.newTransaction(1, 2);
  imap[wtxn_amt] = 30;
  wtxn_amt = gl.newTransaction(2, 0);
  imap[wtxn_amt] = 35;

  auto balances = gl.initBalances({20, 10, 0}, imap);
  auto selected = gl.initSelectedSet(imap);
  auto new_selected = gl.updateSelectedTransactions(balances, selected);

  auto circ = gl.getCircuit();
  for (auto w : new_selected) {
    circ.setAsOutput(w);
  }

  auto res = circ.evaluate(imap);
  BOOST_TEST(res[0] == 0);
}

BOOST_AUTO_TEST_SUITE_END()
