#include <string>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TCut.h"
#include "TChain.h"
#include "TMath.h"

using namespace std;

// helper for making the output name
std::string createOutputName(std::string& name, std::string& trailer){
  std::string outputName = name.substr(0, name.size() - 5);
  outputName += trailer;
  return outputName;
}

// helper for opening the file
TFile* openOutput(std::string& tree, std::string& input, std::string& output) {
  TFile* outFile  =new TFile(output.c_str(),"RECREATE");
  std::cout << "Reading: " << tree << " from " << input  << " to " << output << std::endl;
  return outFile;
}

// helper for finalizing
void finalize(TTree* cutTree, TFile* output){
  TTree*  newtree = cutTree->CloneTree(-1);
  newtree->Write();
  output->Close();
}

/* Now the main business */
string preprocess(string fileName, bool sideband = false, string trailer = "_cut.root", string treeName = "DecayTree"){

  // get the input
  TChain* decaytree = new TChain(treeName.c_str());
  decaytree->Add(fileName.c_str());
	
  // make the output file name 
  string outputName = createOutputName(fileName, trailer);
 
  // make the output file
  TFile* outFile = openOutput(treeName,fileName,outputName); 

  
  /*** Cuts to the data ***/

  // prepare range of cuts to the data
  string cutString1 = "(mjpipi - 3096.90 - mpipi < 500)";                               // Q factor
  string cutString2 = "(!(((mppp > 5264) && (mppp < 5294)) && (kaon_ProbNNpi > 0.4)))"; // Kaon -> Pion misidentification in B mass region  
  string cutString3 = "(!((scaledmass > 5400) && ((mjpk > 5264) && (mjpk < 5294))))";   // Remove B -> J/psi K* signal from upper sideband  
  string cutString4 = sideband ? "" : "(scaledmass < 5400)"                             // Remove sideband
  string cutString  = cutString1 + " && " + cutString2 + " && " + cutString3 + " && " + cutString4; // Combine
  // apply to data
  TCut cut = cutString.c_str();
  TTree* newtree = decaytree->CopyTree(cut);
  std::cout << " Selected " <<  decaytree->GetEntries() << " " << newtree->GetEntries() << std::endl;
  


  /*** New Branches ***/

  // access required branches
  double pi1_ip; newtree->SetBranchAddress("pion1_IPCHI2_OWNPV"  , &pi1_ip);
  double pi2_ip; newtree->SetBranchAddress("pion2_IPCHI2_OWNPV"  , &pi2_ip);
  double k_ip;   newtree->SetBranchAddress("kaon_IPCHI2_OWNPV"   , &k_ip);
  double mup_ip; newtree->SetBranchAddress("muplus_IPCHI2_OWNPV" , &mup_ip);
  double mum_ip; newtree->SetBranchAddress("muminus_IPCHI2_OWNPV", &mum_ip);

  double bp_ip; newtree->SetBranchAddress("bplus_IPCHI2_OWNPV", &bp_ip);
  double bp_fd; newtree->SetBranchAddress("bplus_FDCHI2_OWNPV", &bp_fd);

  double mup_pt; newtree->SetBranchAddress("muplus_PT" , &mup_pt);
  double mum_pt; newtree->SetBranchAddress("muminus_PT", &mum_pt);

  // add new branches
  double ip1, log_ip1; TBranch* branch_log_ip1 = newtree->Branch("log_1_IPCHI2_OWNPV", &log_ip1, "log_1_IPCHI2_OWNPV/D");
  double ip2, log_ip2; TBranch* branch_log_ip2 = newtree->Branch("log_2_IPCHI2_OWNPV", &log_ip2, "log_2_IPCHI2_OWNPV/D");
  double ip3, log_ip3; TBranch* branch_log_ip3 = newtree->Branch("log_3_IPCHI2_OWNPV", &log_ip3, "log_3_IPCHI2_OWNPV/D");
  double ip4, log_ip4; TBranch* branch_log_ip4 = newtree->Branch("log_4_IPCHI2_OWNPV", &log_ip4, "log_4_IPCHI2_OWNPV/D");
  double ip5, log_ip5; TBranch* branch_log_ip5 = newtree->Branch("log_5_IPCHI2_OWNPV", &log_ip5, "log_5_IPCHI2_OWNPV/D");

  double log_ipb; TBranch* branch_log_ipb = newtree->Branch("log_bplus_IPCHI2_OWNPV", &log_ipb, "log_bplus_IPCHI2_OWNPV/D");
  double log_fdb; TBranch* branch_log_fdb = newtree->Branch("log_bplus_FDCHI2_OWNPV", &log_fdb, "log_bplus_FDCHI2_OWNPV/D");

  double mu_pt_max; TBranch* branch_mu_pt_max = newtree->Branch("mu_PT_max", &mu_pt_max, "mu_PT_max/D");
  double mu_pt_min; TBranch* branch_mu_pt_min = newtree->Branch("mu_PT_min", &mu_pt_min, "mu_PT_min/D");

  // loop through entries in tree
  int num_entries = newtree->GetEntries();
  for (int i = 0; i < num_entries; ++i) {
    // access current entry
    newtree->GetEntry(i);

    // store IPCHI2 values in vector to sort
    vector<double> ipchi2s;
    ipchi2s.push_back(pi1_ip);
    ipchi2s.push_back(pi2_ip);
    ipchi2s.push_back(k_ip);
    ipchi2s.push_back(mup_ip);
    ipchi2s.push_back(mum_ip);
    sort(ipchi2s.begin(), ipchi2s.end());

    // convert to log values and store
    ip1 = ipchi2s[0]; log_ip1 = TMath::Log(ip1); branch_log_ip1->Fill();
    ip2 = ipchi2s[1]; log_ip2 = TMath::Log(ip2); branch_log_ip2->Fill();
    ip3 = ipchi2s[2]; log_ip3 = TMath::Log(ip3); branch_log_ip3->Fill();
    ip4 = ipchi2s[3]; log_ip4 = TMath::Log(ip4); branch_log_ip4->Fill();
    ip5 = ipchi2s[4]; log_ip5 = TMath::Log(ip5); branch_log_ip5->Fill();

    // convert bplus vars to log and store
    log_ipb = TMath::Log(bp_ip); branch_log_ipb->Fill();
    log_fdb = TMath::Log(bp_fd); branch_log_fdb->Fill();

    // store muon pt as max and min
    mu_pt_max = max(mup_pt, mum_pt); branch_mu_pt_max->Fill();
    mu_pt_min = min(mup_pt, mum_pt); branch_mu_pt_min->Fill();
  }

  finalize(newtree,outFile);  

  return outputName;
}

