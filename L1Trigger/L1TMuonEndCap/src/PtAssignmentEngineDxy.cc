#include "L1Trigger/L1TMuonEndCap/interface/PtAssignmentEngineDxy.h"

#include <cassert>
#include <iostream>
#include <sstream>

#include "helper.h"  // assert_no_abort


PtAssignmentEngineDxy::PtAssignmentEngineDxy() {

  std::string cmssw_base_ = std::getenv("CMSSW_BASE");


  pbFileNameDxy_ = "/src/L1Trigger/L1TMuonEndCap/data/emtfpp_tf_graphs/model_graph.displ.3.pb";
  pbFileNameDxy_ = cmssw_base_ + pbFileNameDxy_;
  inputNameDxy_ = "batch_normalization_1_input";
  outputNamesDxy_ = {"dense_5/BiasAdd"};

  graphDefDxy_ = tensorflow::loadGraphDef(pbFileNameDxy_);
  assert(graphDefDxy_ != nullptr);
  sessionDxy_ = tensorflow::createSession(graphDefDxy_);
  assert(sessionDxy_ != nullptr);
}

PtAssignmentEngineDxy::~PtAssignmentEngineDxy() {

  tensorflow::closeSession(sessionDxy_);
  delete graphDefDxy_;
}

void PtAssignmentEngineDxy::configure(
    int verbose
) {
  verbose_          = verbose;

  // std::string cmssw_base_ = std::getenv("CMSSW_BASE");

  // graphDefDxy_      = graphDefDxy;
  // sessionDxy_       = sessionDxy;
  // pbFileNameDxy_    = pbFileNameDxy;
  // inputNameDxy_     = inputNameDxy;
  // outputNamesDxy_   = outputNamesDxy;

  // pbFileNameDxy_ = "/src/L1Trigger/L1TMuonEndCap/data/emtfpp_tf_graphs/model_graph.displ.3.pb";
  // pbFileNameDxy_ = cmssw_base_ + pbFileNameDxy_;
  // inputNameDxy_ = "batch_normalization_1_input";
  // outputNamesDxy_ = {"dense_5/BiasAdd"};

  // graphDefDxy_ = tensorflow::loadGraphDef(pbFileNameDxy_);
  // assert(graphDefDxy_ != nullptr);
  // sessionDxy_ = tensorflow::createSession(graphDefDxy_);
  // assert(sessionDxy_ != nullptr);
}

const PtAssignmentEngineDxyAux& PtAssignmentEngineDxy::aux() const {
  static const PtAssignmentEngineDxyAux instance;
  return instance;
}

void PtAssignmentEngineDxy::calculate_pt_dxy(const EMTFTrack& track, emtf::Feature& feature, emtf::Prediction& prediction) const {
  // This is called for each track instead of for entire track collection as was done in Phase-2 implementation 
  preprocessing_dxy(track, feature);
  call_tensorflow_dxy(feature, prediction);
  postprocessing_dxy(track, feature, prediction);
  
  return;
}

void PtAssignmentEngineDxy::preprocessing_dxy(const EMTFTrack& track, emtf::Feature& feature) const {
    static std::array<float, emtf::NLAYERS> x_phi;   // delta-phis = (raw phis - road_phi_median)
    static std::array<float, emtf::NLAYERS> x_theta; // raw thetas
    static std::array<float, emtf::NLAYERS> x_bend;
    static std::array<float, emtf::NLAYERS> x_qual;
    static std::array<float, emtf::NLAYERS> x_time;

    // Initialize to zeros
    x_phi.fill(0);
    x_theta.fill(0);
    x_bend.fill(0);
    x_qual.fill(0);
    x_time.fill(0);

    // Set the values

    // Determine hit_lay of EMTF++ from station, ring information
    for (const auto& hit : track.Hits()) {
      int32_t hit_lay = 99;
      if (hit.Station() == 1) {
        if (hit.Is_CSC()) {
          if (hit.Ring() == 1 || hit.Ring() == 4) {
            hit_lay = 0;            
          } else if (hit.Ring() == 2 || hit.Ring() == 3) {
            hit_lay = 1;
          }
        } else if (hit.Is_RPC()) {
            hit_lay = 5;
        }
      } else if (hit.Station() == 2) {
        if (hit.Is_CSC()) {
          hit_lay = 2;
        } else if (hit.Is_RPC()) {
          hit_lay = 6;
        }
      } else if (hit.Station() == 3) {
        if (hit.Is_CSC()) {
          hit_lay = 3;
        } else if (hit.Is_RPC()) {
          hit_lay = 7;
        }
      } else if (hit.Station() == 4) {
        if (hit.Is_CSC()) {
          hit_lay = 4;
        } else if (hit.Is_RPC()) {
          hit_lay = 8;
        }
      }

      // Drop iRPC
      if ((hit_lay == 7 || hit_lay == 8) && hit.Ring() == 1)
        continue;
      // Drop GE1/1, GE2/1, ME0, DT and skip failed layer mapping
      if (hit_lay >= 9)
        continue;

      assert(std::abs(x_phi.at(hit_lay)) < 1e-7);   // sanity check
      x_phi[hit_lay] = hit.Phi_fp();  // uses old_emtf_phi
      assert(std::abs(x_theta.at(hit_lay)) < 1e-7); // sanity check
      x_theta[hit_lay] = hit.Theta_fp();
      assert(std::abs(x_bend.at(hit_lay)) < 1e-7);  // sanity check
      x_bend[hit_lay] = hit.Bend();  // uses old_emtf_bend
      assert(std::abs(x_qual.at(hit_lay)) < 1e-7);  // sanity check
      x_qual[hit_lay] = hit.Quality();
      assert(std::abs(x_time.at(hit_lay)) < 1e-7);  // sanity check
      x_time[hit_lay] = hit.Time();
    }

    // Mimic Phase-1 EMTF input calculations
    // 6 delta Phis: S1-S2, S1-S3, S1-S4, S2-S3, S2-S4, S3-S4
    // 6 delta Thetas: S1-S2, S1-S3, S1-S4, S2-S3, S2-S4, S3-S4
    // 4 bends : set to zero if no CSC hit and thus RPC hit is used
    // 1 FR bit: for ME1 only
    // 1 Ring bit: for ME1 only
    // 1 track Theta taken from stub coordinate in ME2, ME3, ME4 (in this priority)
    // 4 RPC bits indicating if ME or RE hit was used in each station (S1, S2, S3, S4)
    // Total: 23 variables
    static std::array<float, 6> x_dphi;
    static std::array<float, 6> x_dtheta;

    static std::array<float, 4> x_phi_emtf;   // temporary
    static std::array<float, 4> x_theta_emtf; // temporary
    static std::array<float, 4> x_bend_emtf;
    static std::array<float, 1> x_fr_emtf;
    static std::array<float, 1> x_trk_theta;
    static std::array<float, 1> x_me11ring;
    static std::array<float, 4> x_rpcbit;

    // Initialize to zeros
    x_dphi.fill(0);
    x_dtheta.fill(0);
    //
    x_phi_emtf.fill(0);
    x_theta_emtf.fill(0);
    x_bend_emtf.fill(0);
    x_fr_emtf.fill(0);
    x_trk_theta.fill(0);
    x_me11ring.fill(0);
    x_rpcbit.fill(0);

    // Station 1
    if (x_theta[0] > 1e-7) {  // ME1/1
      x_phi_emtf[0]   = x_phi[0];
      x_theta_emtf[0] = x_theta[0];
      x_bend_emtf[0]  = x_bend[0];
      x_fr_emtf[0]    = (x_qual[0] > 0) ? 1 : -1;
    } else if (x_theta[1] > 1e-7) {  // ME1/2
      x_phi_emtf[0]   = x_phi[1];
      x_theta_emtf[0] = x_theta[1];
      x_bend_emtf[0]  = x_bend[1];
      x_fr_emtf[0]    = (x_qual[1] > 0) ? 1 : -1;
      x_me11ring[0]   = 1;
    } else if (x_theta[5] > 1e-7) {  // RE1
      x_phi_emtf[0]   = x_phi[5];
      x_theta_emtf[0] = x_theta[5];
      x_bend_emtf[0]  = 0;
      x_rpcbit[0]     = 1;
    }

    // Station 2
    if (x_theta[2] > 1e-7) {  // ME2
      x_phi_emtf[1]   = x_phi[2];
      x_theta_emtf[1] = x_theta[2];
      x_bend_emtf[1]  = x_bend[2];
    } else if (x_theta[6] > 1e-7) {  // RE2
      x_phi_emtf[1]   = x_phi[6];
      x_theta_emtf[1] = x_theta[6];
      x_bend_emtf[1]  = 0;
      x_rpcbit[1]     = 1;
    }

    // Station 3
    if (x_theta[3] > 1e-7) {  // ME3
      x_phi_emtf[2]   = x_phi[3];
      x_theta_emtf[2] = x_theta[3];
      x_bend_emtf[2]  = x_bend[3];
    } else if (x_theta[7] > 1e-7) {  // RE3
      x_phi_emtf[2]   = x_phi[7];
      x_theta_emtf[2] = x_theta[7];
      x_bend_emtf[2]  = 0;
      x_rpcbit[2]     = 1;
    }

    // Station 4
    if (x_theta[4] > 1e-7) {  // ME4
      x_phi_emtf[3]   = x_phi[4];
      x_theta_emtf[3] = x_theta[4];
      x_bend_emtf[3]  = x_bend[4];
    } else if (x_theta[8] > 1e-7) {  // RE4
      x_phi_emtf[3]   = x_phi[8];
      x_theta_emtf[3] = x_theta[8];
      x_bend_emtf[3]  = 0;
      x_rpcbit[3]     = 1;
    }

    // Set x_trk_theta
    if (x_theta[2] > 1e-7) {
      x_trk_theta[0] = x_theta[2];
    } else if (x_theta[3] > 1e-7) {
      x_trk_theta[0] = x_theta[3];
    } else if (x_theta[4] > 1e-7) {
      x_trk_theta[0] = x_theta[4];
    }

    // Set x_dphi, x_dtheta
    auto calc_delta = [](float a, float b) {
      if ((a > 1e-7) && (b > 1e-7)) {
        return static_cast<float>(a-b);
      }
      return 0.f;
    };

    x_dphi[0] = calc_delta(x_phi_emtf[0], x_phi_emtf[1]);
    x_dphi[1] = calc_delta(x_phi_emtf[0], x_phi_emtf[2]);
    x_dphi[2] = calc_delta(x_phi_emtf[0], x_phi_emtf[3]);
    x_dphi[3] = calc_delta(x_phi_emtf[1], x_phi_emtf[2]);
    x_dphi[4] = calc_delta(x_phi_emtf[1], x_phi_emtf[3]);
    x_dphi[5] = calc_delta(x_phi_emtf[2], x_phi_emtf[3]);

    x_dtheta[0] = calc_delta(x_theta_emtf[0], x_theta_emtf[1]);
    x_dtheta[1] = calc_delta(x_theta_emtf[0], x_theta_emtf[2]);
    x_dtheta[2] = calc_delta(x_theta_emtf[0], x_theta_emtf[3]);
    x_dtheta[3] = calc_delta(x_theta_emtf[1], x_theta_emtf[2]);
    x_dtheta[4] = calc_delta(x_theta_emtf[1], x_theta_emtf[3]);
    x_dtheta[5] = calc_delta(x_theta_emtf[2], x_theta_emtf[3]);

    // Pack the 23 features + 13 zeros used for padding
    feature = {{
      x_dphi     [0], x_dphi     [1], x_dphi     [2], x_dphi     [3], x_dphi     [4], x_dphi     [5],
      x_dtheta   [0], x_dtheta   [1], x_dtheta   [2], x_dtheta   [3], x_dtheta   [4], x_dtheta   [5],
      x_bend_emtf[0], x_bend_emtf[1], x_bend_emtf[2], x_bend_emtf[3], x_fr_emtf  [0], x_trk_theta[0],
      x_me11ring [0], x_rpcbit   [0], x_rpcbit   [1], x_rpcbit   [2], x_rpcbit   [3], 0             ,
      0             , 0             , 0             , 0             , 0             , 0             ,
      0             , 0             , 0             , 0             , 0             , 0
    }};
    return;
  }

  void PtAssignmentEngineDxy::call_tensorflow_dxy(const emtf::Feature& feature, emtf::Prediction& prediction) const {
    const int nfeatures_displ = 23;    // 23 features
    static tensorflow::Tensor input(tensorflow::DT_FLOAT, { 1, nfeatures_displ });
    static std::vector<tensorflow::Tensor> outputs;
    //assert(feature.size() == NFEATURES);

    float* d = input.flat<float>().data();
    //std::copy(feature.begin(), feature.end(), d);
    std::copy(feature.begin(), feature.begin() + nfeatures_displ, d);
    tensorflow::run(sessionDxy_, { { inputNameDxy_, input } }, outputNamesDxy_, &outputs);
    assert(outputs.size() == 1);
    //assert(prediction.size() == NPREDICTIONS);

    const float reg_pt_scale = 100.;  // a scale factor applied to regression during training
    const float reg_dxy_scale = 1.0;  // a scale factor applied to regression during training

    prediction.at(0) = outputs[0].matrix<float>()(0, 0);
    prediction.at(1) = outputs[0].matrix<float>()(0, 1);

    // Remove scale factor used during training
    prediction.at(0) /= reg_pt_scale;
    prediction.at(1) /= reg_dxy_scale;
    return;
  }

  void PtAssignmentEngineDxy::postprocessing_dxy(const EMTFTrack& track, const emtf::Feature& feature, emtf::Prediction& prediction) const {
    // Demote tracks with less than 3 stations (CSC/RPC)
    bool demote = false;
    int track_mode = 0;

    for (const auto& hit : track.Hits()) {
      int32_t hit_lay = 99;
      if (hit.Station() == 1) {
        if (hit.Is_CSC()) {
          if (hit.Ring() == 1 || hit.Ring() == 4) {
            hit_lay = 0;            
          } else if (hit.Ring() == 2 || hit.Ring() == 3) {
            hit_lay = 1;
          }
        } else if (hit.Is_RPC()) {
            hit_lay = 5;
        }
      } else if (hit.Station() == 2) {
        if (hit.Is_CSC()) {
          hit_lay = 2;
        } else if (hit.Is_RPC()) {
          hit_lay = 6;
        }
      } else if (hit.Station() == 3) {
        if (hit.Is_CSC()) {
          hit_lay = 3;
        } else if (hit.Is_RPC()) {
          hit_lay = 7;
        }
      } else if (hit.Station() == 4) {
        if (hit.Is_CSC()) {
          hit_lay = 4;
        } else if (hit.Is_RPC()) {
          hit_lay = 8;
        }
      } 

      // Drop iRPC
      if ((hit_lay == 7 || hit_lay == 8) && hit.Ring() == 1)
        continue;
      // Drop GE1/1, GE2/1, ME0, DT and skip failed layer mapping
      if (hit_lay >= 9)
        continue;


      if (hit_lay == 0 || hit_lay == 1 || hit_lay == 5) {
        track_mode |= (1 << 3);
      } else if (hit_lay == 2 || hit_lay == 6) {
        track_mode |= (1 << 2);
      } else if (hit_lay == 3 || hit_lay == 7) {
        track_mode |= (1 << 1);
      } else if (hit_lay == 4 || hit_lay == 8) {
        track_mode |= (1 << 0);
      }
    }

    if (!(track_mode == 11 || track_mode == 13 || track_mode == 14 || track_mode == 15)) {
      demote = true;
    }
    if (demote) {
      prediction.at(0) = 0.5;  // demote to 2 GeV
    }
    return;
  }
