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

int skim_large_root(){

  TString file1 = "samples/DisplacedMuGun_flatPt2to1000_posEndcap_flatD0Z0EtaPhi_11_3_0_pre5_2M5.root";
  TString file2 = "samples/DisplacedMuGun_flatPt2to1000_negEndcap_flatD0Z0EtaPhi_11_3_0_pre5.root";

  // load trees
  TString tree = "EMTFNtuple/tree";

  TFile *fout =new TFile("DisplacedMuGun_flatPt2to1000_flatD0Z0EtaPhi_11_3_0_pre5_2M5_skimmed.root","RECREATE");
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
  TTreeReaderArray<float  > gmtMuonPt(reader,"gmtMuon_pt");
  TTreeReaderArray<float  > gmtMuonPtDxy(reader,"gmtMuon_pt_dxy");
  TTreeReaderArray<short  > gmtMuonDxy(reader,"gmtMuon_dxy");
  TTreeReaderArray<float  > gmtMuonPhi(reader,"gmtMuon_phi");
  TTreeReaderArray<float  > gmtMuonEta(reader,"gmtMuon_eta");
  TTreeReaderArray<short  > gmtMuonCharge(reader,"gmtMuon_q");
  TTreeReaderArray<short  > gmtMuonQual(reader,"gmtMuon_qual");

  TTreeReaderArray<short  > genPartCharge(reader,"genPart_q");
  TTreeReaderArray<float  > genPartPt(reader,"genPart_pt");
  TTreeReaderArray<float  > genPartDxy(reader,"genPart_dxy");
  TTreeReaderArray<float  > genPartEta(reader,"genPart_eta");
  TTreeReaderArray<float  > genPartPhi(reader,"genPart_phi");
  TTreeReaderArray<short  > genPartID(reader,"genPart_ID");
  TTreeReaderArray<int32_t> genPartParentID(reader,"genPart_parentID");
  TTreeReaderArray<float  > genPartVx(reader,"genPart_vx");
  TTreeReaderArray<float  > genPartVy(reader,"genPart_vy");
  TTreeReaderArray<float  > genPartVz(reader,"genPart_vz");

  std::vector<float> emtfTrack_pt, emtfTrack_pt_dxy, emtfTrack_dxy, emtfTrack_phi, emtfTrack_theta, emtfTrack_eta, gmtMuon_pt, gmtMuon_pt_dxy, gmtMuon_phi, gmtMuon_eta, gen_pt, gen_dxy, gen_d0, gen_eta, gen_phi, gen_etaStar, gen_phiStar, gen_vx, gen_vy, gen_vz, gen_dR, dR_gmt_emtf;

  std::vector<int>   emtfTrack_size, emtfTrack_GMT_phi, emtfTrack_GMT_eta, emtfTrack_mode, emtfTrack_endcap, emtfTrack_sector, emtfTrack_bx, gmtMuon_size, gmtMuon_dxy, gmtMuon_q, gmtMuon_qual, gen_q, gen_ID, gen_parentID, gen_gmt_idx, gmt_emtf_idx;

  t1->Branch("dR_gmt_emtf",&dR_gmt_emtf);
  t1->Branch("gmt_emtf_idx",&gmt_emtf_idx);

  t1->Branch("emtfTrack_size",&emtfTrack_size);
  t1->Branch("emtfTrack_pt",&emtfTrack_pt);
  t1->Branch("emtfTrack_pt_dxy",&emtfTrack_pt_dxy);
  t1->Branch("emtfTrack_dxy",&emtfTrack_dxy);
  t1->Branch("emtfTrack_phi",&emtfTrack_phi);
  t1->Branch("emtfTrack_theta",&emtfTrack_theta);
  t1->Branch("emtfTrack_eta",&emtfTrack_eta);
  t1->Branch("emtfTrack_GMT_phi",&emtfTrack_GMT_phi);
  t1->Branch("emtfTrack_GMT_eta",&emtfTrack_GMT_eta);
  t1->Branch("emtfTrack_mode",&emtfTrack_mode);
  t1->Branch("emtfTrack_endcap",&emtfTrack_endcap);
  t1->Branch("emtfTrack_sector",&emtfTrack_sector);
  t1->Branch("emtfTrack_bx",&emtfTrack_bx);

  t1->Branch("gmtMuon_size",&gmtMuon_size);
  t1->Branch("gmtMuon_pt",&gmtMuon_pt);
  t1->Branch("gmtMuon_pt_dxy",&gmtMuon_pt_dxy);
  t1->Branch("gmtMuon_dxy",&gmtMuon_dxy);
  t1->Branch("gmtMuon_phi",&gmtMuon_phi);
  t1->Branch("gmtMuon_eta",&gmtMuon_eta);
  t1->Branch("gmtMuon_q",&gmtMuon_q);
  t1->Branch("gmtMuon_qual",&gmtMuon_qual);

  t1->Branch("gen_q",&gen_q);
  t1->Branch("gen_pt",&gen_pt);
  t1->Branch("gen_dxy",&gen_dxy);
  t1->Branch("gen_d0",&gen_d0);
  t1->Branch("gen_eta",&gen_eta);
  t1->Branch("gen_phi",&gen_phi);
  t1->Branch("gen_eta",&gen_etaStar);
  t1->Branch("gen_phi",&gen_phiStar);
  t1->Branch("gen_ID",&gen_ID);
  t1->Branch("gen_parentID",&gen_parentID);
  t1->Branch("gen_vx",&gen_vx);
  t1->Branch("gen_vy",&gen_vy);
  t1->Branch("gen_vz",&gen_vz);
  t1->Branch("gen_dR",&gen_dR);

//   float dR, l1mu_idx, dR_EMTF, idx_EMTF;
//   std::vector<float> dR_gmt;
//   std::vector<int> gmt_idx;

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

    // clear stuff
    emtfTrack_pt.clear();
    emtfTrack_pt_dxy.clear();
    emtfTrack_dxy.clear();
    emtfTrack_phi.clear();
    emtfTrack_theta.clear();
    emtfTrack_eta.clear();
    emtfTrack_size.clear();
    emtfTrack_GMT_phi.clear();
    emtfTrack_GMT_eta.clear();
    emtfTrack_mode.clear();
    emtfTrack_endcap.clear();
    emtfTrack_sector.clear();
    emtfTrack_bx.clear();

    gmtMuon_pt.clear();
    gmtMuon_pt_dxy.clear();
    gmtMuon_dxy.clear();
    gmtMuon_phi.clear();
    gmtMuon_eta.clear();
    gmtMuon_size.clear();
    gmtMuon_q.clear();
    gmtMuon_qual.clear();
    
    gen_pt.clear();
    gen_dxy.clear();
    gen_d0.clear();
    gen_eta.clear();
    gen_phi.clear();
    gen_etaStar.clear();
    gen_phiStar.clear();
    gen_vx.clear();
    gen_vy.clear();
    gen_vz.clear();
    gen_dR.clear();
    gen_q.clear();
    gen_ID.clear();
    gen_parentID.clear();
    gen_gmt_idx.clear();
    gmt_emtf_idx.clear();

    gen_gmt_idx.clear();
    dR_gmt_emtf.clear();
    
    

    

    // loop over GEN muons
    int ngen = 0;
    bool skip_flag = false;
    
    int i = 0;
    //   std::cout << "genPartParentID[i] = " << genPartParentID[i] << std::endl;
    if(abs(genPartID[i]) != 13 ) continue;
    //   if(abs(genPartParentID[i]) != 6000113 ) continue;
    ngen++;

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

//   float Lxy = TMath::Sqrt(genPartVx[i]*genPartVx[i] + genPartVy[i]*genPartVy[i]);

    

    // GEN muon GMT muon matching and GMT muon EMTF track matching
    float dR = 5.0;
    int idx = -99;
    float dR_EMTF = 5.0;
    int idx_EMTF = -99;
    // Loop over gmt muons and find matching gen muon
    for (unsigned int j=0; j<*gmtMuonSize; j++){
    float dR_new = TMath::Sqrt((gmtMuonEta[j]-etaStar_gen)*(gmtMuonEta[j]-etaStar_gen)+DPhi(gmtMuonPhi[j],phiStar_gen)*DPhi(gmtMuonPhi[j],phiStar_gen));
    if(dR_new > dR){
        continue;
    }
    else{
        dR = dR_new;
        idx = j;
    }
    }
    if (dR > 0.6) {
        continue;
    }
    
    // Find EMTF tracks matched to GMT muon and fill L1 muon info 
    for (int i=0; i<*gmtMuonSize; i++){
      if (idx != i){
        continue;
      } 

      
    //   GMT muon EMTF track matching
      for (int j=0; j<*emtfTrackSize; j++){
        // if(emtfTrackBX[j] != 0) continue;

        float GMTEta = emtfTrackGMTEta[j] * 0.010875;

        int globPhi = (emtfTrackSector[j] - 1) * 96 + emtfTrackGMTPhi[j];
        globPhi = (globPhi + 600) % 576;
        float GMTPhi = globPhi * 0.010908;

        float dR_new_EMTF = TMath::Sqrt((gmtMuonEta[i]-GMTEta)*(gmtMuonEta[i]-GMTEta)+DPhi(gmtMuonPhi[i],GMTPhi)*DPhi(gmtMuonPhi[i],GMTPhi));

        if(dR_new_EMTF > dR_EMTF){
          continue;
        }
        else{
          idx_EMTF = j;
          dR_EMTF = dR_new_EMTF;
        }

      }

    //   if (dR_EMTF > 1.0) {
    //       continue;
    //   }

    if (idx_EMTF < 0) {
        continue;
    }
    int GMT_dxy = -1;

    if (abs(emtfTrackDxy[idx_EMTF]) < 25){
        GMT_dxy = 0;
    }
    else if (abs(emtfTrackDxy[idx_EMTF]) < 50){
        GMT_dxy = 1;
    }
    else if (abs(emtfTrackDxy[idx_EMTF]) < 75){
        GMT_dxy = 2;
    }
    else {
        GMT_dxy = 3;
    }
        

    gen_dR.push_back(dR);
    gen_gmt_idx.push_back(idx);
    
    gen_pt.push_back(genPartPt[i]);
//   std::cout << "gen_pt["<<i<<"] = " << genPartPt[i] << std::endl;
    gen_eta.push_back(genPartEta[i]);
    gen_phi.push_back(genPartPhi[i]);
//   gen_Lxy.push_back(Lxy);
    gen_etaStar.push_back(etaStar_gen);
    gen_phiStar.push_back(phiStar_gen);
    gen_vx.push_back(genPartVx[i]);
    gen_vy.push_back(genPartVy[i]);
    gen_vz.push_back(genPartVz[i]);
    gen_parentID.push_back(genPartParentID[i]);
    gen_ID.push_back(i);
    gen_d0.push_back(d0);
    gen_q.push_back(genPartCharge[i]);

    gmt_emtf_idx.push_back(idx_EMTF);
    dR_gmt_emtf.push_back(dR_EMTF);

    emtfTrack_pt.push_back(emtfTrackPt[idx_EMTF]);
    emtfTrack_pt_dxy.push_back(emtfTrackPtDxy[idx_EMTF]);
    emtfTrack_eta.push_back(gmtMuonEta[i]);
    emtfTrack_phi.push_back(gmtMuonPhi[i]);
    // emtfTrack_qual.push_back(gmtMuonQual[i]);
    emtfTrack_dxy.push_back(emtfTrackDxy[idx_EMTF]);

        
    }

    t1->Fill();

  } // end event loop

  t1->Write();

  return 0;
} // end function
