#pragma once

#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <emp-tool/emp-tool.h>

#include <algorithm>
#include <memory>
#include <string>
#include <unordered_map>

#include "../io/netmp.h"
#include "../utils/circuit.h"
#include "jump_provider.h"
#include "ot_provider.h"
#include "preproc.h"
#include "quadsquad/rand_gen_pool.h"
#include "sharing.h"
#include "types.h"

namespace quadsquad {
class OfflineEvaluator {
  int id_;
  int security_param_;
  RandGenPool rgen_;
  std::shared_ptr<io::NetIOMP<4>> network_;
  std::shared_ptr<io::NetIOMP<4>> network_ot_;
  utils::LevelOrderedCircuit circ_;
  std::shared_ptr<ThreadPool> tpool_;
  PreprocCircuit<Ring> preproc_;
  std::vector<std::unique_ptr<OTProvider>> ot_;
  JumpProvider jump_;
  NTL::ZZ_pContext ZZ_p_ctx_;
  NTL::ZZ_pEContext ZZ_pE_ctx_;

  // Data members used for book-keeping across methods.
  std::vector<utils::FIn2Gate> mult_gates_;
  std::array<std::vector<Ring>, 3> ab_terms_;
  std::array<std::vector<Ring>, 6> c_terms_;

  NTL::Vec<NTL::ZZ_pE> zk_prove_;
  std::array<NTL::Vec<NTL::ZZ_pE>, 3> zk_verify_;
  std::array<NTL::Vec<NTL::ZZ_pE>, 3> zk_verify_inputs_;
  std::vector<std::future<void>> zk_res_;
  std::array<NTL::ZZ_pE, 3> zk_exp_val_;
  std::array<std::vector<NTL::ZZ_pE>, 3> zk_check_;
  NTL::ZZ_pE x_interp_;
  std::array<NTL::ZZ_pE, 3> lagrange_coeff_d3_;
  NTL::Mat<NTL::ZZ_pE> opp_check_shares_;

  // Verifier computation for prover 'prover_id' for an iteration in the
  // recursive distributed ZKP.
  void zkVerifyRecursiveIter(int prover_id,
                             const NTL::Vec<NTL::ZZ_pE>& poly_p_shares);

  // Verifier computation for prover 'prover_id' for the base distributed ZKP.
  void zkVerifyBase(int prover_id, emp::block cc_key,
                    const NTL::Vec<NTL::ZZ_pE>& pi);

  // Used for running common coin protocol. Returns common random PRG key which
  // is then used to generate randomness for common coin output.
  emp::block commonCoinKey();

 public:
  // `network_1` and `network_2` are required to be distinct.
  // `network_2` is used for OT while `network_1` is used for all other tasks.
  OfflineEvaluator(int my_id, std::shared_ptr<io::NetIOMP<4>> network1,
                   std::shared_ptr<io::NetIOMP<4>> network2,
                   utils::LevelOrderedCircuit circ, int security_param,
                   int threads, int seed = 200);

  // Generate sharing of a random unknown value.
  static void randomShare(RandGenPool& rgen, ReplicatedShare<Ring>& share);
  // Generate sharing of a random value known to dealer (called by all parties
  // except the dealer).
  static void randomShareWithParty(int id, int dealer, RandGenPool& rgen,
                                   ReplicatedShare<Ring>& share);
  // Generate sharing of a random value known to party. Should be called by
  // dealer when other parties call other variant.
  static void randomShareWithParty(int id, RandGenPool& rgen,
                                   ReplicatedShare<Ring>& share, Ring& secret);

  // Following methods implement various preprocessing subprotocols.

  // Set masks for each wire. Should be called before running any of the other
  // subprotocols.
  void setWireMasks(
      const std::unordered_map<utils::wire_t, int>& input_pid_map);
  // Computes S_1 and S_2 summands.
  void computeABCrossTerms();
  // Computes S_0 summands by running instances of disMult.
  void computeCCrossTerms();
  // Combines all computed summands to create output shares. Should be called
  // after 'computeABCrossTerms' and 'computeCCrossTerms' terminate.
  void combineCrossTerms();
  // Runs distributed ZKP to verify behaviour in 'computeABCrossTerms'. Should
  // be called after 'computeABCrossTerms' terminates.
  void distributedZKP();
  // Compute output commitments. Should be called after 'combineCrossTerms'.
  void computeOutputCommitments();

  PreprocCircuit<Ring> getPreproc();

  // Efficiently runs above subprotocols.
  PreprocCircuit<Ring> run(
      const std::unordered_map<utils::wire_t, int>& input_pid_map);

  // Insecure preprocessing. All preprocessing data is generated in clear but
  // cast in a form that can be used in the online phase.
  static PreprocCircuit<Ring> dummy(
      const utils::LevelOrderedCircuit& circ,
      const std::unordered_map<utils::wire_t, int>& input_pid_map,
      size_t security_param, int pid, emp::PRG& prg);
};
};  // namespace quadsquad
