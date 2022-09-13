#include <io/netmp.h>
#include <quadsquad/offline_evaluator.h>
#include <utils/circuit.h>

#include <algorithm>
#include <boost/program_options.hpp>
#include <cmath>
#include <iostream>
#include <memory>

#include "utils.h"

using namespace quadsquad;
using json = nlohmann::json;
namespace bpo = boost::program_options;

utils::Circuit<Ring> generateCircuit(size_t num_mult_gates) {
  utils::Circuit<Ring> circ;

  std::vector<utils::wire_t> inputs(num_mult_gates);
  std::generate(inputs.begin(), inputs.end(),
                [&]() { return circ.newInputWire(); });

  std::vector<utils::wire_t> outputs(num_mult_gates);
  for (size_t i = 0; i < num_mult_gates - 1; ++i) {
    outputs[i] = circ.addGate(utils::GateType::kMul, inputs[i], inputs[i + 1]);
  }
  outputs[num_mult_gates - 1] = circ.addGate(
      utils::GateType::kMul, inputs[num_mult_gates - 1], inputs[0]);

  return circ;
}

void benchmark(const bpo::variables_map& opts) {
  bool save_output = false;
  std::string save_file;
  if (opts.count("output") != 0) {
    save_output = true;
    save_file = opts["output"].as<std::string>();
  }

  auto gates = opts["gates"].as<size_t>();
  auto pid = opts["pid"].as<size_t>();
  auto security_param = opts["security-param"].as<size_t>();
  auto cm_threads = opts["cm-threads"].as<size_t>();
  auto cp_threads = opts["cp-threads"].as<size_t>();
  auto seed = opts["seed"].as<size_t>();
  auto repeat = opts["repeat"].as<size_t>();
  auto port = opts["port"].as<int>();

  std::shared_ptr<io::NetIOMP<4>> network1 = nullptr;
  std::shared_ptr<io::NetIOMP<4>> network2 = nullptr;
  if (opts["localhost"].as<bool>()) {
    network1 = std::make_shared<io::NetIOMP<4>>(pid, port, nullptr, true);
    network2 = std::make_shared<io::NetIOMP<4>>(pid, port + 100, nullptr, true);
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

    network1 = std::make_shared<io::NetIOMP<4>>(pid, port, ip.data(), false);
    network2 =
        std::make_shared<io::NetIOMP<4>>(pid, port + 100, ip.data(), false);
  }

  json output_data;
  output_data["details"] = {{"gates", gates},
                            {"pid", pid},
                            {"security_param", security_param},
                            {"cm_threads", cm_threads},
                            {"cp_threads", cp_threads},
                            {"seed", seed},
                            {"repeat", repeat}};
  output_data["benchmarks"] = json::array();

  std::cout << "--- Details ---\n";
  for (const auto& [key, value] : output_data["details"].items()) {
    std::cout << key << ": " << value << "\n";
  }
  std::cout << std::endl;

  auto circ = generateCircuit(gates).orderGatesByLevel();

  std::unordered_map<utils::wire_t, int> input_pid_map;
  for (const auto& g : circ.gates_by_level[0]) {
    if (g->type == utils::GateType::kInp) {
      input_pid_map[g->out] = 0;
    }
  }

  initNTL(cp_threads);

  for (size_t r = 0; r < repeat; ++r) {
    OfflineEvaluator eval(pid, network1, network2, circ, security_param,
                          cm_threads, seed);

    network1->sync();
    network2->sync();

    network1->resetStats();
    network2->resetStats();
    CommPoint net1_st(*network1);
    CommPoint net2_st(*network2);
    TimePoint start;
    eval.run(input_pid_map);
    TimePoint end;
    CommPoint net1_ed(*network1);
    CommPoint net2_ed(*network2);

    auto time = end - start;
    auto res1 = net1_ed - net1_st;
    auto res2 = net2_ed - net2_st;
    for (size_t i = 0; i < res1.size(); ++i) {
      res1[i] += res2[i];
    }

    nlohmann::json rbench = {{"time", time}};
    rbench["communication"] = res1;

    size_t bytes_sent = 0;
    for (auto val : res1) {
      bytes_sent += val;
    }

    std::cout << "--- Repetition " << r + 1 << " ---\n";
    std::cout << "time: " << time << " ms\n";
    std::cout << "sent: " << bytes_sent << " bytes\n";
    std::cout << "throughput: " << (gates * 1e3) / time
              << " triples per second\n";

    output_data["benchmarks"].push_back(std::move(rbench));
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
    ("gates,g", bpo::value<size_t>()->required(), "Number of multiplication gates.")
    ("pid,p", bpo::value<size_t>()->required(), "Party ID.")
    ("security-param", bpo::value<size_t>()->default_value(128), "Security parameter in bits.")
    ("cm-threads", bpo::value<size_t>()->default_value(7), "Number of threads for communication (recommended value is at least 7).")
    ("cp-threads", bpo::value<size_t>()->default_value(1), "Number of threads for computation (recommended value close to number of cores).")
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
      "Benchmark triple generation in offline phase.");
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
