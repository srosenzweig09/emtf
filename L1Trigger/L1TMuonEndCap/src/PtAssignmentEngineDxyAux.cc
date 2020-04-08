#include "L1Trigger/L1TMuonEndCap/interface/PtAssignmentEngineDxyAux.h"


int PtAssignmentEngineDxyAux::getGMTPtDxy(float pt) const {
// compressed pt = pt*1 (scale) + 1 (pt = 0 is empty candidate)
  int gmt_pt_dxy = (pt * 1) + 1;
  gmt_pt_dxy = (gmt_pt_dxy > 255) ? 255 : gmt_pt_dxy;
  return gmt_pt_dxy;
}

float PtAssignmentEngineDxyAux::getPtFromGMTPtDxy(int gmt_pt_dxy) const {
  float pt = (gmt_pt_dxy <= 0) ?  0 : 1.0 * (gmt_pt_dxy-1);
  return pt;
}

int PtAssignmentEngineDxyAux::getGMTDxy(float dxy) const {
  int gmt_dxy = 0;
  if (std::abs(dxy) < 20.) {
    gmt_dxy = 0;
  } else if (std::abs(dxy) < 50.) {
    gmt_dxy = 1;
  } else if (std::abs(dxy) < 80.) {
    gmt_dxy = 2;
  } else {
    gmt_dxy = 3;
  }
  return gmt_dxy;
}

float PtAssignmentEngineDxyAux::getTriggerPtDxy(float pt_pred) const {
  // Used for get_trigger_pt()
  float s_min {0.};
  float s_max {0.};
  int   s_nbins {0};
  float s_step {0.};
  std::array<float, 120> s_lut {};
  std::array<float, 120> s_lut_wp50 {};

  s_min   = 0.;
  s_max   = 60.;
  s_nbins = 120;
  s_step  = (s_max - s_min)/float(s_nbins);
  s_lut   = {{  1.8219,  1.5725,  1.6237,  1.8672,  2.2298,  2.6692,  3.1724,  3.7226,
                4.3005,  4.8945,  5.4972,  6.1065,  6.7217,  7.3443,  7.9797,  8.6255,
                9.2779,  9.9339, 10.5868, 11.2340, 11.8745, 12.5056, 13.1352, 13.7707,
               14.3985, 15.0216, 15.6519, 16.2946, 16.9297, 17.5598, 18.1981, 18.8567,
               19.5345, 20.2317, 20.9849, 21.7932, 22.5650, 23.2764, 23.9233, 24.5326,
               25.1879, 25.9589, 26.8144, 27.6406, 28.3182, 28.9110, 29.4817, 30.0894,
               30.8001, 31.5674, 32.3055, 33.0457, 33.8479, 34.6975, 35.4941, 36.2179,
               36.9157, 37.6592, 38.5602, 39.6237, 40.7733, 41.9798, 43.2775, 44.6862,
               45.9872, 46.8917, 47.5905, 48.2057, 48.8099, 49.4649, 50.1705, 50.8610,
               51.5614, 52.2918, 53.0282, 53.7657, 54.5035, 55.2414, 55.9793, 56.7173,
               57.4553, 58.1933, 58.9314, 59.6694, 60.4074, 61.1454, 61.8834, 62.6214,
               63.3594, 64.0973, 64.8353, 65.5733, 66.3113, 67.0493, 67.7873, 68.5253,
               69.2633, 70.0012, 70.7392, 71.4772, 72.2152, 72.9532, 73.6912, 74.4292,
               75.1671, 75.9051, 76.6431, 77.3811, 78.1191, 78.8571, 79.5950, 80.3330,
               81.0710, 81.8090, 82.5470, 83.2849, 84.0229, 84.7609, 85.4989, 86.2369}};
  s_lut_wp50 = {{1.8124,  1.6471,  1.6755,  1.8367,  2.0737,  2.3461,  2.6370,  2.9560,
                3.3248,  3.7365,  4.1760,  4.6336,  5.1040,  5.5834,  6.0695,  6.5604,
                7.0554,  7.5540,  8.0544,  8.5545,  9.0529,  9.5514, 10.0488, 10.5407,
               11.0263, 11.5075, 11.9870, 12.4668, 12.9474, 13.4297, 13.9161, 14.4090,
               14.9068, 15.4037, 15.8966, 16.3903, 16.8852, 17.3796, 17.8709, 18.3599,
               18.8473, 19.3375, 19.8375, 20.3540, 20.8927, 21.4490, 21.9967, 22.5160,
               23.0021, 23.4527, 23.8652, 24.2528, 24.6402, 25.0503, 25.4903, 25.9606,
               26.4660, 27.0031, 27.5589, 28.1126, 28.6454, 29.1493, 29.6322, 30.1029,
               30.5670, 31.0276, 31.4843, 31.9233, 32.3456, 32.7724, 33.2167, 33.6778,
               34.1510, 34.6287, 35.1127, 35.6217, 36.1572, 36.7039, 37.2606, 37.8230,
               38.3763, 38.9074, 39.4167, 39.9213, 40.4378, 40.9845, 41.5990, 42.2614,
               42.9157, 43.5348, 44.1085, 44.6446, 45.1498, 45.6289, 46.0819, 46.5207,
               46.9573, 47.3828, 47.7878, 48.1767, 48.5567, 48.9351, 49.3208, 49.7180,
               50.1278, 50.5593, 51.0135, 51.4887, 51.9777, 52.4705, 52.9646, 53.4593,
               53.9542, 54.4493, 54.9446, 55.4399, 55.9353, 56.4307, 56.9261, 57.4215}};
  assert(s_lut.size() == (size_t) s_nbins);
  assert(s_lut_wp50.size() == (size_t) s_nbins);

  float xml_pt = std::abs(1.0/pt_pred);
  if (xml_pt <= 2.) {  // do not use the LUT if below 2 GeV
    return xml_pt;
  }


  // xml_pt = std::clamp(xml_pt, s_min, static_cast<float>(s_max - 1e-5));
  if (xml_pt > static_cast<float>(s_max - 1e-5)) xml_pt = static_cast<float>(s_max - 1e-5);
  else if (xml_pt < s_min) xml_pt = s_min;

  xml_pt = (xml_pt - s_min) / (s_max - s_min) * float(s_nbins);  // convert to bin number
  int binx = static_cast<int>(xml_pt);
  binx = (binx == s_nbins-1) ? (binx-1) : binx;  // avoid boundary

  // int binx = digitize(xml_pt);
  float x0 = float(binx) * s_step;
  float x1 = float(binx+1) * s_step;
  float y0 = s_lut.at(binx);
  float y1 = s_lut.at(binx+1);
  // float trg_pt = interpolate(xml_pt, x0, x1, y0, y1);
  float trg_pt = (xml_pt - x0) / (x1 - x0) * (y1 - y0) + y0;
  return trg_pt;
}