#include <io/netmp.h>
#include <quadsquad/offline_evaluator.h>
#include <quadsquad/online_evaluator.h>
#include <utils/circuit.h>
#include <utils/liquidity_matching.h>

#include <algorithm>
#include <boost/program_options.hpp>
#include <cmath>
#include <iostream>
#include <memory>
#include <random>
#include <unordered_map>

#include "utils.h"

using namespace quadsquad;
using json = nlohmann::json;
namespace bpo = boost::program_options;

std::unordered_map<utils::wire_t, Ring> generateRandomTxns(
    utils::SoDoGridLock<Ring>& gl, size_t num_banks, size_t num_txns,
    size_t seed) {
  std::unordered_map<utils::wire_t, Ring> imap;

  std::mt19937 gen(seed);
  std::uniform_int_distribution<> amt_dist(1, 1000);
  std::uniform_int_distribution<> bank_dist(0, static_cast<int>(num_banks) - 1);

  // Initialize transactions.
  for (size_t i = 0; i < num_txns; i++) {
    auto src = bank_dist(gen);
    auto dest = bank_dist(gen);
    auto amt = amt_dist(gen);
    auto wtxn_amt = gl.newTransaction(src, dest);
    imap[wtxn_amt] = static_cast<Ring>(amt);
  }

  // Initialize balances.
  std::vector<Ring> balances(num_banks);
  for (size_t i = 0; i < num_banks; i++) {
    balances[i] = amt_dist(gen);
  }
  auto wbalances = gl.initBalances(balances, imap);
  auto wselected = gl.initSelectedSet(imap);
  gl.updateSelectedTransactions(wbalances, wselected);

  return imap;
}

void benchmark(const bpo::variables_map& opts) {
  bool save_output = false;
  std::string save_file;
  if (opts.count("output") != 0) {
    save_output = true;
    save_file = opts["output"].as<std::string>();
  }

  auto pid = opts["pid"].as<size_t>();
  auto security_param = opts["security-param"].as<size_t>();
  auto threads = opts["threads"].as<size_t>();
  auto seed = opts["seed"].as<size_t>();
  auto repeat = opts["repeat"].as<size_t>();
  auto port = opts["port"].as<int>();
  auto num_banks = opts["num-banks"].as<size_t>();
  auto num_txns = opts["num-txns"].as<size_t>();

  std::shared_ptr<io::NetIOMP<4>> network = nullptr;
  if (opts["localhost"].as<bool>()) {
    network = std::make_shared<io::NetIOMP<4>>(pid, port, nullptr, true);
  } else {
    std::ifstream fnet(opts["net-config"].as<std::string>());
    if (!fnet.good()) {
      fnet.close();
      throw std::runtime_error("Could not open network config file");
    }
    json netdata;
    fnet >> netdata;
    fnet.close();

    std::vector<std::string> ipaddress(4);
    std::array<char*, 4> ip{};
    for (size_t i = 0; i < 4; ++i) {
      ipaddress[i] = netdata[i].get<std::string>();
      ip[i] = ipaddress[i].data();
    }

    network = std::make_shared<io::NetIOMP<4>>(pid, port, ip.data(), false);
  }

  json output_data;
  output_data["details"] = {{"pid", pid},
                            {"security_param", security_param},
                            {"threads", threads},
                            {"seed", seed},
                            {"num_banks", num_banks},
                            {"repeat", repeat},
                            {"num_txns", num_txns}};
  output_data["benchmarks"] = json::array();

  std::cout << "--- Details ---\n";
  for (const auto& [key, value] : output_data["details"].items()) {
    std::cout << key << ": " << value << "\n";
  }
  std::cout << std::endl;

  utils::SoDoGridLock<Ring> gl(num_banks);
  auto imap = generateRandomTxns(gl, num_banks, num_txns, seed);
  auto circ = gl.getCircuit().orderGatesByLevel();

  std::unordered_map<utils::wire_t, int> input_pid_map;
  for (const auto& it : imap) {
    input_pid_map[it.first] = 0;
  }

  std::cout << "--- Circuit ---\n";
  std::cout << circ << std::endl;

  emp::PRG prg(&emp::zero_block, seed);

  for (size_t r = 0; r < repeat; ++r) {
    auto preproc =
        OfflineEvaluator::dummy(circ, input_pid_map, security_param, pid, prg);

    OnlineEvaluator eval(pid, network, std::move(preproc), circ, security_param,
                         threads, seed);

    network->sync();

    eval.setInputs(imap);
    StatsPoint start(*network);
    for (size_t i = 0; i < circ.gates_by_level.size(); ++i) {
      eval.evaluateGatesAtDepth(i);
    }
    StatsPoint end(*network);

    auto rbench = end - start;
    output_data["benchmarks"].push_back(rbench);

    size_t bytes_sent = 0;
    for (const auto& val : rbench["communication"]) {
      bytes_sent += val.get<int64_t>();
    }

    std::cout << "--- Repetition " << r + 1 << " ---\n";
    std::cout << "time: " << rbench["time"] << " ms\n";
    std::cout << "sent: " << bytes_sent << " bytes\n";

    if (save_output) {
      saveJson(output_data, save_file);
    }
    std::cout << std::endl;
  }

  output_data["stats"] = {{"peak_virtual_memory", peakVirtualMemory()},
                          {"peak_resident_set_size", peakResidentSetSize()}};

  std::cout << "--- Statistics ---\n";
  for (const auto& [key, value] : output_data["stats"].items()) {
    std::cout << key << ": " << value << "\n";
  }
  std::cout << std::endl;

  if (save_output) {
    saveJson(output_data, save_file);
  }
}

// clang-format off
bpo::options_description programOptions() {
  bpo::options_description desc("Following options are supported by config file too.");
  desc.add_options()
    ("num-banks,b", bpo::value<size_t>()->required(), "Number of banks.")
    ("num-txns,x", bpo::value<size_t>()->required(), "Number of transactions.")
    ("pid,p", bpo::value<size_t>()->required(), "Party ID.")
    ("security-param", bpo::value<size_t>()->default_value(128), "Security parameter in bits.")
    ("threads,t", bpo::value<size_t>()->default_value(6), "Number of threads (recommended 6).")
    ("seed", bpo::value<size_t>()->default_value(200), "Value of the random seed.")
    ("net-config", bpo::value<std::string>(), "Path to JSON file containing network details of all parties.")
    ("localhost", bpo::bool_switch(), "All parties are on same machine.")
    ("port", bpo::value<int>()->default_value(10000), "Base port for networking.")
    ("output,o", bpo::value<std::string>(), "File to save benchmarks.")
    ("repeat,r", bpo::value<size_t>()->default_value(1), "Number of times to run benchmarks.");

  return desc;
}
// clang-format on

int main(int argc, char* argv[]) {
  auto prog_opts(programOptions());

  bpo::options_description cmdline(
      "Benchmark online phase for multiplication gates.");
  cmdline.add(prog_opts);
  cmdline.add_options()(
      "config,c", bpo::value<std::string>(),
      "configuration file for easy specification of cmd line arguments")(
      "help,h", "produce help message");

  bpo::variables_map opts;
  bpo::store(bpo::command_line_parser(argc, argv).options(cmdline).run(), opts);

  if (opts.count("help") != 0) {
    std::cout << cmdline << std::endl;
    return 0;
  }

  if (opts.count("config") > 0) {
    std::string cpath(opts["config"].as<std::string>());
    std::ifstream fin(cpath.c_str());

    if (fin.fail()) {
      std::cerr << "Could not open configuration file at " << cpath << "\n";
      return 1;
    }

    bpo::store(bpo::parse_config_file(fin, prog_opts), opts);
  }

  // Validate program options.
  try {
    bpo::notify(opts);

    // Check if output file already exists.
    if (opts.count("output") != 0) {
      std::ifstream ftemp(opts["output"].as<std::string>());
      if (ftemp.good()) {
        ftemp.close();
        throw std::runtime_error("Output file aready exists.");
      }
      ftemp.close();
    }

    if (!opts["localhost"].as<bool>() && (opts.count("net-config") == 0)) {
      throw std::runtime_error("Expected one of 'localhost' or 'net-config'");
    }
  } catch (const std::exception& ex) {
    std::cerr << ex.what() << std::endl;
    return 1;
  }

  try {
    benchmark(opts);
  } catch (const std::exception& ex) {
    std::cerr << ex.what() << "\nFatal error" << std::endl;
    return 1;
  }

  return 0;
}
