#ifndef L1TMuonEndCap_PtAssignmentEngineDxyAux_h
#define L1TMuonEndCap_PtAssignmentEngineDxyAux_h

#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <array>

class PtAssignmentEngineDxyAux {
public:
  // Functions for GMT quantities
  int getGMTPtDxy(float pt) const;

  float getPtFromGMTPtDxy(int gmt_pt_dxy) const;

  int getGMTDxy(float dxy) const;

  float getTriggerPtDxy(float pt_pred) const;
};

#endif
