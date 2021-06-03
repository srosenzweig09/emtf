/* L1T analysis package 
 * module: main module of Orthogonal dataset method 
 * Author: George Karathanasis, georgios.karathanasis@cern.ch
 *
 *
 */

#include <iostream>
#include "TH1F.h"
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
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
// #include "analysis.h"
// #include "CutReader.h"

float DPhi(double phi1,double phi2){
  float temp=phi1-phi2;
  if (temp>3.14) temp=temp-6.28;
  if (temp<-3.14) temp=temp+6.28;
  return temp;
}

int skim_dm_root(){

  TString file1 = "samples/v6/DisplacedMuGun_flatPt2to1000_negEndcap_flatXYZEtaPhi_11_3_0_pre5_NNv6_2M.root";
  TString file2 = "samples/v6/DisplacedMuGun_flatPt2to1000_posEndcap_flatXYZEtaPhi_11_3_0_pre5_NNv6_2M.root";

  // load trees
  TString tree = "EMTFNtuple/tree";

  TFile *fout =new TFile("samples/v6/DisplacedMuGun_flatPt2to1000_flatXYZEtaPhi_11_3_0_pre5_NNv6_skimmed.root","RECREATE");
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

  std::vector<float> emtfTrack_pt, emtfTrack_pt_dxy, emtfTrack_dxy, emtfTrack_phi, emtfTrack_eta, gen_pt, gen_d0, gen_eta, gen_phi, gen_etaStar, gen_phiStar;

  std::vector<int>   gen_q;

  t1->Branch("emtfTrack_pt",&emtfTrack_pt);
  t1->Branch("emtfTrack_pt_dxy",&emtfTrack_pt_dxy);
  t1->Branch("emtfTrack_dxy",&emtfTrack_dxy);
  t1->Branch("emtfTrack_phi",&emtfTrack_phi);
  t1->Branch("emtfTrack_eta",&emtfTrack_eta);

  t1->Branch("gen_q",&gen_q);
  t1->Branch("gen_pt",&gen_pt);
  t1->Branch("gen_d0",&gen_d0);
  t1->Branch("gen_eta",&gen_eta);
  t1->Branch("gen_phi",&gen_phi);
  t1->Branch("gen_etaStar",&gen_etaStar);
  t1->Branch("gen_phiStar",&gen_phiStar);

  int eventCount = 0;
  float z_ME2 = 808.0; // ME2 z coordinate in [cm]
  float r = 0.;
  float rStar = 0.;
  float etaStar_gen = 0.;
  float phiStar_gen = 0.;

  while(reader.Next()){
    if (eventCount % 10000 == 0) {
      std::cout << eventCount << " events read!" << std::endl;
    } 
    // if (eventCount > 100) {
    //     continue;
    // }
    eventCount++;

    if (*emtfTrackSize < 1) {
        continue;
    }

    if (*gmtMuonSize < 1) {
        continue;
    }

    
    int i = 0;
    if (genPartVz[i] > 100) {
        continue;
    }

    // std::cout << eventCount << " events read!" << std::endl;

    // clear stuff
    gen_pt.clear();
    gen_eta.clear();
    gen_phi.clear();
    gen_etaStar.clear();
    gen_phiStar.clear();
    gen_d0.clear();
    gen_q.clear();

    emtfTrack_pt.clear();
    emtfTrack_pt_dxy.clear();
    emtfTrack_eta.clear();
    emtfTrack_phi.clear();
    emtfTrack_dxy.clear();

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

    bool flag;
    if (*emtfTrackSize == 2) {

        int j = i + 1;

        if (emtfTrackMode[j] > emtfTrackMode[i]) {
            i = j;
        }

    } 

    if (flag) {std::cout << i << " = i" << std::endl;}

    gen_pt.push_back(genPartPt[i]);
    gen_eta.push_back(genPartEta[i]);
    gen_phi.push_back(genPartPhi[i]);
    gen_etaStar.push_back(etaStar_gen);
    gen_phiStar.push_back(phiStar_gen);
    gen_d0.push_back(d0);
    gen_q.push_back(genPartCharge[i]);

    
    // std::cout << gmtMuonEta[i] << " gmtMuonEta[i]" << std::endl;

    emtfTrack_pt.push_back(emtfTrackPt[i]);
    emtfTrack_pt_dxy.push_back(emtfTrackPtDxy[i]);
    emtfTrack_eta.push_back(gmtMuonEta[i]);
    emtfTrack_phi.push_back(gmtMuonPhi[i]);
    emtfTrack_dxy.push_back(emtfTrackDxy[i]);

    t1->Fill();

  } // end event loop

  t1->Write();

  return 0;
} // end function
