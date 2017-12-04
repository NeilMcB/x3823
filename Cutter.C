#include <string>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TCut.h"
#include "TChain.h"
#include "TMath.h"
/* 
*  Example of how to write a reduced ntuple
*  Two selections, Matts (add on to roadmap):
*  Use: .L ApplySelection.C+
*        ApplySelection("theFile.root")
* Giving extra arguments allows to change ntuple and output name
*        ApplySelection("EtaDown.root", "rhoGammaDown/rhoGammaDown", "_Selected.root");
*/ 


/** First some helpers to avoid duplicating code */

std::string createOutputName(std::string& name, std::string& trailer){
  // helper for making the output name
  std::string outputName = name.substr(0,name.size() - 5);
  outputName += trailer;
  return outputName;
}

TFile* openOutput(std::string& tree, std::string& input, std::string& output) {
  // helper for opening the file
  TFile* outFile  =new TFile(output.c_str(),"RECREATE");
  std::cout << "Reading: " << tree << " from " << input  << " to " << output << std::endl;
  return outFile;
}

void finalize(TTree* cutTree, TFile* output){
  // helper for finalizing
  TTree*  newtree = cutTree->CloneTree(-1);
  newtree->Write();
  output->Close();
}

/* Now the main business */

std::string Cutter(std::string fileName = "mcall.root", std::string treeName = "psiCand", std::string trailer = "_cut.root", std::string cutString  = ""){

  // Matts selection, generally applied on top of the road map

  // get the input
  TChain* decaytree = new TChain(treeName.c_str());
  decaytree->Add(fileName.c_str());
	
  // make the output file name 
  std::string outputName = createOutputName(fileName, trailer);
 
  // make the output file
  TFile* outFile = openOutput(treeName,fileName,outputName); 
  TCut cut = cutString.c_str();

  TTree* smalltree = decaytree->CopyTree(cut);
  std::cout << " Selected " <<  decaytree->GetEntries() << " " << smalltree->GetEntries() << std::endl;
  
  finalize(smalltree,outFile);  

  
  return outputName;
}

