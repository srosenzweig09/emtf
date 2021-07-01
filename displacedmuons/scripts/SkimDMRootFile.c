/* Project: EMTF / Displaced Muons
 * Author: Suzanne Rosenzweig
 * This script skims a large collection of displaced muon gun samples and saves only the events with an emtf muon that can be matched to a gmt/gen muon.
 *
 */

#include <iostream>
// #include "TH1F.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TChain.h"
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <vector>
#include "TString.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
// #include "TCanvas.h"
// #include "TH1F.h"
// #include "TH2F.h"
// #include "analysis.h"
// #include "CutReader.h"

float DPhi(double phi1,double phi2){
  float temp=phi1-phi2;
  if (temp>3.14) temp=temp-6.28;
  if (temp<-3.14) temp=temp+6.28;
  return temp;
}

int SkimDMRootFile(int argc, char* argv[]){

  bool v8 = false;
  int min = atoi(argv[2]);
  int max = atoi(argv[3]);

  TString file1, file2;

  if (v8) {
    file1 = "samples/v8/DisplacedMuGun_flatPt2to1000_negEndcap_flatXYZEtaPhi_11_3_0_pre5_NNv8_1M.root";
    file2 = "samples/v8/DisplacedMuGun_flatPt2to1000_posEndcap_flatXYZEtaPhi_11_3_0_pre5_NNv8_1M.root";

    TFile *fout =new TFile("samples/v8/DisplacedMuGun_flatPt2to1000_flatXYZEtaPhi_11_3_0_pre5_NNv8_skimmed.root","RECREATE");
  }
  else {
    file1 = "samples/v6/DisplacedMuGun_flatPt2to1000_negEndcap_flatXYZEtaPhi_11_3_0_pre5_NNv6_2M.root";
    file2 = "samples/v6/DisplacedMuGun_flatPt2to1000_posEndcap_flatXYZEtaPhi_11_3_0_pre5_NNv6_2M.root";

    // TFile *fout =new TFile("samples/v6/DisplacedMuGun_flatPt2to1000_flatXYZEtaPhi_11_3_0_pre5_NNv6_skimmed_1.root","RECREATE");
    TFile *fout =new TFile(argv[1],"RECREATE");
  }

  // load trees
  TString tree = "EMTFNtuple/tree";
  
  TTree * t1 =new TTree("tree","tree");

  TChain * cc=new TChain(tree);

  cc->Add(file1);
  cc->Add(file2);

  TTreeReader reader(cc);

  TTreeReaderValue<int32_t> emtfTrackSize(reader,"emtfTrack_size");
  TTreeReaderArray<float  > emtfTrackPt(reader,"emtfTrack_pt");
  TTreeReaderArray<float  > emtfTrackPtDxy(reader,"emtfTrack_pt_dxy");
  TTreeReaderArray<float  > emtfTrackDxy(reader,"emtfTrack_dxy");
  TTreeReaderArray<float  > emtfTrackPhi(reader,"emtfTrack_phi");
  TTreeReaderArray<float  > emtfTrackTheta(reader,"emtfTrack_theta");
  TTreeReaderArray<float  > emtfTrackEta(reader,"emtfTrack_eta");
  TTreeReaderArray<int    > emtfTrackGMTPhi(reader,"emtfTrack_GMT_phi");
  TTreeReaderArray<int    > emtfTrackGMTEta(reader,"emtfTrack_GMT_eta");
  TTreeReaderArray<short  > emtfTrackMode(reader,"emtfTrack_mode");
  TTreeReaderArray<short  > emtfTrackEndcap(reader,"emtfTrack_endcap");
  TTreeReaderArray<short  > emtfTrackSector(reader,"emtfTrack_sector");
  TTreeReaderArray<short  > emtfTrackBX(reader,"emtfTrack_bx");

  TTreeReaderValue<int32_t> gmtMuonSize(reader,"gmtMuon_size");
  TTreeReaderArray<float  > gmtMuonPhi(reader,"gmtMuon_phi");
  TTreeReaderArray<float  > gmtMuonEta(reader,"gmtMuon_eta");

  TTreeReaderArray<short  > genPartCharge(reader,"genPart_q");
  TTreeReaderArray<float  > genPartPt(reader,"genPart_pt");
  TTreeReaderArray<float  > genPartDxy(reader,"genPart_dxy");
  TTreeReaderArray<float  > genPartEta(reader,"genPart_eta");
  TTreeReaderArray<float  > genPartPhi(reader,"genPart_phi");
  TTreeReaderArray<float  > genPartVx(reader,"genPart_vx");
  TTreeReaderArray<float  > genPartVy(reader,"genPart_vy");
  TTreeReaderArray<float  > genPartVz(reader,"genPart_vz");

  std::vector<float> emtfTrack_pt, emtfTrack_pt_dxy, emtfTrack_dxy, emtfTrack_phi, emtfTrack_eta, emtfTrack_dR, gen_pt, gen_d0, gen_vz, gen_eta, gen_phi, gen_etaStar, gen_phiStar, gen_dR, matched_gen_pt, matched_gen_d0, matched_gen_vz, matched_gen_eta, matched_gen_phi, matched_gen_etaStar, matched_gen_phiStar, matched_gen_dR;

  std::vector<int>   gen_q, matched_gen_q;

  t1->Branch("emtfTrack_pt",&emtfTrack_pt);
  t1->Branch("emtfTrack_pt_dxy",&emtfTrack_pt_dxy);
  t1->Branch("emtfTrack_dxy",&emtfTrack_dxy);
  t1->Branch("emtfTrack_phi",&emtfTrack_phi);
  t1->Branch("emtfTrack_eta",&emtfTrack_eta);
  t1->Branch("emtfTrack_dR",&emtfTrack_dR);

  t1->Branch("gen_q",&gen_q);
  t1->Branch("gen_pt",&gen_pt);
  t1->Branch("gen_d0",&gen_d0);
  t1->Branch("gen_vz",&gen_vz);
  t1->Branch("gen_eta",&gen_eta);
  t1->Branch("gen_phi",&gen_phi);
  t1->Branch("gen_etaStar",&gen_etaStar);
  t1->Branch("gen_phiStar",&gen_phiStar);
  t1->Branch("gen_dR",&gen_dR);

  t1->Branch("matched_gen_q",&matched_gen_q);
  t1->Branch("matched_gen_pt",&matched_gen_pt);
  t1->Branch("matched_gen_d0",&matched_gen_d0);
  t1->Branch("matched_gen_vz",&matched_gen_vz);
  t1->Branch("matched_gen_eta",&matched_gen_eta);
  t1->Branch("matched_gen_phi",&matched_gen_phi);
  t1->Branch("matched_gen_etaStar",&matched_gen_etaStar);
  t1->Branch("matched_gen_phiStar",&matched_gen_phiStar);
  t1->Branch("matched_gen_dR",&matched_gen_dR);

  int matchedMu = 0;
  float z_ME2 = 808.0; // ME2 z coordinate in [cm]
  float r = 0.;
  float rStar = 0.;
  float etaStar_gen = 0.;
  float phiStar_gen = 0.;

  int eventCount = 0;
  int saveCount = 0;

  while(reader.Next()){
    if (eventCount % 10000 == 0) {
      std::cout << eventCount << " events read!" << std::endl;
    } 
    eventCount++;
    if (eventCount > max) {break;}
    if (eventCount < min) {continue;}
    saveCount++;
    
    int i = 0;

    // if (genPartVz[i] > 100) {
    //     continue;
    // }

    // std::cout << eventCount << " events read!" << std::endl;

    // clear stuff
    gen_pt.clear();
    gen_eta.clear();
    gen_phi.clear();
    gen_etaStar.clear();
    gen_phiStar.clear();
    gen_d0.clear();
    gen_dR.clear();
    gen_q.clear();
    gen_vz.clear();

    matched_gen_pt.clear();
    matched_gen_eta.clear();
    matched_gen_phi.clear();
    matched_gen_etaStar.clear();
    matched_gen_phiStar.clear();
    matched_gen_d0.clear();
    matched_gen_dR.clear();
    matched_gen_q.clear();
    matched_gen_vz.clear();

    emtfTrack_pt.clear();
    emtfTrack_pt_dxy.clear();
    emtfTrack_eta.clear();
    emtfTrack_phi.clear();
    emtfTrack_dxy.clear();
    emtfTrack_dR.clear();


    // convert displaced eta and phi to prompt-like and calculate dR between gen and l1 muon
    if (genPartEta[i] > 0) r = abs(z_ME2 - genPartVz[i])/abs(TMath::SinH(genPartEta[i]));
    else r = abs(-z_ME2 - genPartVz[i])/abs(TMath::SinH(genPartEta[i]));     
    float xStar = genPartVx[i] + r * TMath::Cos(genPartPhi[i]);
    float yStar = genPartVy[i] + r * TMath::Sin(genPartPhi[i]);
    rStar = TMath::Sqrt(xStar * xStar + yStar * yStar);
    
    etaStar_gen = TMath::ASinH(z_ME2/rStar) * (genPartEta[i]/abs(genPartEta[i]));

    // if (abs(etaStar_gen) < 1.2 || abs(etaStar_gen) > 2.5 ) continue;

    if (xStar >= 0) phiStar_gen = TMath::ATan(yStar/xStar); 
    else if (yStar >= 0 && xStar < 0) phiStar_gen = TMath::Pi() + TMath::ATan(yStar/xStar); 
    else if (yStar <= 0 && xStar < 0) phiStar_gen = TMath::ATan(yStar/xStar) - TMath::Pi(); 

    float d0 = -999.0;
    float invPt = genPartCharge[i]/genPartPt[i];
    invPt = (std::abs(invPt) < 1./10000) ? (invPt < 0 ? -1./10000 : +1./10000) : invPt;
    double R = -1.0 / (0.003 * 3.811 * invPt);
    float xc = genPartVx[i] - (R * std::sin(genPartPhi[i]));
    float yc = genPartVy[i] + (R * std::cos(genPartPhi[i]));
    d0 = R - ((R < 0 ? -1. : +1.) * TMath::Sqrt(xc*xc + yc*yc));

    gen_pt.push_back(genPartPt[i]);
    gen_vz.push_back(genPartVz[i]);
    gen_eta.push_back(genPartEta[i]);
    gen_phi.push_back(genPartPhi[i]);
    gen_etaStar.push_back(etaStar_gen);
    gen_phiStar.push_back(phiStar_gen);
    gen_d0.push_back(d0);
    gen_q.push_back(genPartCharge[i]);

    bool flag;
    if (*emtfTrackSize == 2) {

        int j = i + 1;

        if (emtfTrackMode[j] > emtfTrackMode[i]) {
            i = j;
        }

    } else {

      for (int j = 0; j<*emtfTrackSize; j++) {
        if (emtfTrackMode[j] > emtfTrackMode[i]) {
          i = j;
        }
      }

    }

    if (*emtfTrackSize < 1) {
      gen_dR.push_back(-1);
      t1->Fill();
      continue;
    }

    if (*gmtMuonSize < 1) {
      gen_dR.push_back(-1);
      t1->Fill();
      continue;
    }
    
    // // Quality Cuts!
    // if (emtfTrackMode[i] < 11) {
    //   t1->Fill();
    //   continue;
    //   }
    // if (emtfTrackMode[i] == 12) {
    //   t1->Fill();
    //   continue;
    //   }

    float GMTEta = emtfTrackGMTEta[i] * 0.010875;

    int globPhi = (emtfTrackSector[i] - 1) * 96 + emtfTrackGMTPhi[i];

    globPhi = (globPhi + 600) % 576;

    float GMTPhi = globPhi * 0.010908;

    float dR_gmt_emtf = TMath::Sqrt((gmtMuonEta[i]-GMTEta)*(gmtMuonEta[i]-GMTEta)+DPhi(gmtMuonPhi[i],GMTPhi)*DPhi(gmtMuonPhi[i],GMTPhi));

    float dR_gmt_gen  = TMath::Sqrt((gmtMuonEta[i]-etaStar_gen)*(gmtMuonEta[i]-etaStar_gen)+DPhi(gmtMuonPhi[i],phiStar_gen)*DPhi(gmtMuonPhi[i],phiStar_gen));


    // // LOW PT MUONS SNEAKING IN? ADJUST DR CUT AND SEE IF THIS GOES AWAY.
    // if (dR_gmt_gen > 0.6) {
    //   t1->Fill();
    //   continue;
    //   }
    // if (dR_gmt_emtf > 5) {
    //   t1->Fill();
    //   continue;
    // }

    
    // std::cout << gmtMuonEta[i] << " gmtMuonEta[i]" << std::endl;

    matched_gen_pt.push_back(genPartPt[i]);
    matched_gen_vz.push_back(genPartVz[i]);
    matched_gen_eta.push_back(genPartEta[i]);
    matched_gen_phi.push_back(genPartPhi[i]);
    matched_gen_etaStar.push_back(etaStar_gen);
    matched_gen_phiStar.push_back(phiStar_gen);
    matched_gen_d0.push_back(d0);
    matched_gen_dR.push_back(dR_gmt_gen);
    matched_gen_q.push_back(genPartCharge[i]);

    gen_dR.push_back(dR_gmt_gen);

    emtfTrack_pt.push_back(emtfTrackPt[i]);
    emtfTrack_pt_dxy.push_back(emtfTrackPtDxy[i]);
    emtfTrack_eta.push_back(gmtMuonEta[i]);
    emtfTrack_phi.push_back(gmtMuonPhi[i]);
    emtfTrack_dxy.push_back(emtfTrackDxy[i]);
    emtfTrack_dR.push_back(dR_gmt_emtf);

    matchedMu++;



    t1->Fill();

  } // end event loop

  std::cout << "Events read:   " << eventCount << std::endl;
  std::cout << "Events saved:  " << saveCount << std::endl;
  std::cout << "Muons matched: " << matchedMu << std::endl;

  t1->Write();

  return 0;
} // end function


int main(int argc, char* argv[]){
  SkimDMRootFile(argc, argv);
  return 0;
}