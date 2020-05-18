#ifndef L1TMuonEndCap_Common_h
#define L1TMuonEndCap_Common_h

#include "DataFormats/L1TMuon/interface/CPPFDigi.h"
#include "DataFormats/L1TMuon/interface/EMTFHit.h"
#include "DataFormats/L1TMuon/interface/EMTFRoad.h"
#include "DataFormats/L1TMuon/interface/EMTFTrack.h"

#include "L1Trigger/L1TMuon/interface/GeometryTranslator.h"
#include "L1Trigger/L1TMuon/interface/MuonTriggerPrimitive.h"
#include "L1Trigger/L1TMuon/interface/MuonTriggerPrimitiveFwd.h"

#include "L1Trigger/L1TMuon/interface/TTGeometryTranslator.h"
#include "L1Trigger/L1TMuon/interface/TTMuonTriggerPrimitive.h"

#include "L1Trigger/L1TMuonEndCap/interface/EMTFSubsystemTag.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <array>

// Resolve namespaces

typedef l1t::EMTFHit EMTFHit;
typedef l1t::EMTFHitCollection EMTFHitCollection;
typedef l1t::EMTFRoad EMTFRoad;
typedef l1t::EMTFRoadCollection EMTFRoadCollection;
typedef l1t::EMTFTrack EMTFTrack;
typedef l1t::EMTFTrackCollection EMTFTrackCollection;
typedef l1t::EMTFPtLUT EMTFPtLUT;

typedef L1TMuon::GeometryTranslator GeometryTranslator;
typedef L1TMuon::TriggerPrimitive TriggerPrimitive;
typedef L1TMuon::TriggerPrimitiveCollection TriggerPrimitiveCollection;

typedef L1TMuon::TTGeometryTranslator TTGeometryTranslator;
typedef L1TMuon::TTTriggerPrimitive TTTriggerPrimitive;
typedef L1TMuon::TTTriggerPrimitiveCollection TTTriggerPrimitiveCollection;

typedef TriggerPrimitive::CSCData CSCData;
typedef TriggerPrimitive::RPCData RPCData;
typedef TriggerPrimitive::GEMData GEMData;
typedef TTTriggerPrimitive::TTData TTData;

typedef emtf::CSCTag CSCTag;
typedef emtf::RPCTag RPCTag;
typedef emtf::CPPFTag CPPFTag;
typedef emtf::GEMTag GEMTag;
typedef emtf::IRPCTag IRPCTag;
typedef emtf::ME0Tag ME0Tag;
typedef emtf::TTTag TTTag;

namespace emtf {

  // Constants

  // from DataFormats/MuonDetId/interface/CSCDetId.h
  constexpr int MIN_ENDCAP = 1;
  constexpr int MAX_ENDCAP = 2;

  // from DataFormats/MuonDetId/interface/CSCTriggerNumbering.h
  constexpr int MIN_TRIGSECTOR = 1;
  constexpr int MAX_TRIGSECTOR = 6;
  constexpr int NUM_SECTORS = 12;

  // Zones
  constexpr int NUM_ZONES = 4;
  constexpr int NUM_ZONE_HITS = 160;

  // Stations
  constexpr int NUM_STATIONS = 4;
  constexpr int NUM_STATION_PAIRS = 6;

  // Fixed-size arrays
  template <typename T>
  using sector_array = std::array<T, NUM_SECTORS>;
  template <typename T>
  using zone_array = std::array<T, NUM_ZONES>;

  // NN features and predictions
  constexpr int NLAYERS = 16;      // 5 (CSC) + 4 (RPC) + 3 (GEM) + 4 (DT)
  constexpr int NFEATURES = 36;    // NN features
  constexpr int NPREDICTIONS = 9;  // NN outputs

  using Feature = std::array<float, NFEATURES>;
  using Prediction = std::array<float, NPREDICTIONS>;

}  // namespace emtf

#endif
