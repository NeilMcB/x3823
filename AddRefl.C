#include <string>
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TChain.h"
#include "TMath.h"
#include "TLorentzVector.h"

float mPion =  139.57061;
float mJpsi = 3096.870;
float mKaon =  493.667;

TLorentzVector toFourVector(float px, float py, float pz, float m) 
{
  const TVector3 vec = TVector3(px, py, pz);
  return TLorentzVector(vec, TMath::Sqrt(m*m + vec.Mag2()));  
}

std::string AddRefl(std::string fileName = "data.root", std::string treeName = "DecayTree" , std::string trailer = "_reflections.root" )
{
  // get the input
  TChain* decaytree = new TChain(treeName.c_str());
  decaytree->Add(fileName.c_str());
	
  // make the output
  std::string outputName = fileName.substr(0,fileName.size() - 5);
  outputName += trailer;
  TFile* outFile  =new TFile(outputName.c_str(),"RECREATE");
   
  std::cout << "Reading: " << treeName << " from " << fileName  << " to " << outputName << std::endl;

  // clone the tree..
  TTree*  newtree = decaytree->CloneTree(-1);

  // Access required branches
  double k_px;  newtree->SetBranchAddress("bplus_PX_K",&k_px);
  double k_py;  newtree->SetBranchAddress("bplus_PY_K",&k_py);
  double k_pz;  newtree->SetBranchAddress("bplus_PZ_K",&k_pz);

  double jpsi_px;  newtree->SetBranchAddress("bplus_PX_jpsi",&jpsi_px);
  double jpsi_py;  newtree->SetBranchAddress("bplus_PY_jpsi",&jpsi_py);
  double jpsi_pz;  newtree->SetBranchAddress("bplus_PZ_jpsi",&jpsi_pz);
  
  double pion_px1;  newtree->SetBranchAddress("bplus_PX_pi1",&pion_px1);
  double pion_py1;  newtree->SetBranchAddress("bplus_PY_pi1",&pion_py1);
  double pion_pz1;  newtree->SetBranchAddress("bplus_PZ_pi1",&pion_pz1);
  double pion_px2;  newtree->SetBranchAddress("bplus_PX_pi2",&pion_px2);
  double pion_py2;  newtree->SetBranchAddress("bplus_PY_pi2",&pion_py2);
  double pion_pz2;  newtree->SetBranchAddress("bplus_PZ_pi2",&pion_pz2);

  float smass; newtree->SetBranchAddress("scaledmass", &smass);

  int k_ID ; newtree->SetBranchAddress("kaon_ID" ,&k_ID); 
  int p1_ID; newtree->SetBranchAddress("pion1_ID",&p1_ID);
  int p2_ID; newtree->SetBranchAddress("pion2_ID",&p2_ID);

  //add the new branches 
  float r_k2p; TBranch* branch_r_k2p = newtree->Branch("reflection_k2p", &r_k2p, "reflection_k2p/F");
  float r_jkp; TBranch* branch_r_jkp = newtree->Branch("reflection_jkp", &r_jkp, "reflection_jkp/F");  
 
  // loop and make the conversion
  int num_entries  = newtree->GetEntries();
  for (int i = 0; i <num_entries; ++i) {
    // Access current entry
    newtree->GetEntry(i);

    // *** Plot kaon->pion mis-ID *** //
    // Sub in pion mass for Kaon
    TLorentzVector k2pvec = toFourVector(k_px    , k_py    , k_pz    , mPion);
    TLorentzVector pivec1 = toFourVector(pion_px1, pion_py1, pion_pz1, mPion);
    TLorentzVector pivec2 = toFourVector(pion_px2, pion_py2, pion_pz2, mPion);
    TLorentzVector jvec   = toFourVector(jpsi_px , jpsi_py , jpsi_pz , mJpsi);
    // Calculate B 4-vector
    TLorentzVector B_k2p = jvec + k2pvec + pivec1 + pivec2 ; 
    r_k2p = B_k2p.M();
    // Fill branch
    branch_r_k2p->Fill();

    // *** Plot jkp system for lower sideband *** //
    if ((smass > 5202) && (smass < 5250)) {
      TLorentzVector kvec = toFourVector(k_px, k_py, k_pz, mKaon);
      // Pair positive kaon with negative pion
      TLorentzVector pvec;
      if (k_ID ==  321) {
        pvec = (p1_ID == -211) ? toFourVector(pion_px1, pion_py1, pion_pz1, mPion) : toFourVector(pion_px2, pion_py2, pion_pz2, mPion); 
      }
      if (k_ID == -321) {
        pvec = (p1_ID ==  211) ? toFourVector(pion_px1, pion_py1, pion_pz1, mPion) : toFourVector(pion_px2, pion_py2, pion_pz2, mPion);
      }
      // Combine to find mass of system
      TLorentzVector jkp = jvec + kvec + pvec;
      r_jkp = jkp.M();
    }
    else {
      // Dummy value to maintain branch size
      r_jkp = -0.1;
    }
    // Fill branch
    branch_r_jkp->Fill(); 
  }

  // Write to new file and close
  newtree->Write();
  outFile->Close();

  return outputName;
}
