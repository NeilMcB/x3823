#include "TROOT.h"
#include "TMath.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooTruthModel.h"
#include "RooConstVar.h"
#include "RooExponential.h"
#include "RooResolutionModel.h"
#include "RooGaussModel.h"
#include "RooCBShape.h"
#include "RooAddModel.h"
#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooHistPdf.h"
#include "RooProdPdf.h"
#include "RooVoigtian.h"
#include "RooDecay.h"
#include "RooStats/ModelConfig.h"
#include "TCanvas.h"
#include "TTree.h"
#include "RooWorkspace.h"
#include "TCut.h"
#include "TCanvas.h"
#include "RooChebychev.h"
#include "RooDstD0BG.h"
#include "TFile.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "TAxis.h"
#include "RooDataHist.h"
//#include "ksfun.h"
#include "TSystem.h"
#include "RooFitResult.h"
#include "TLatex.h"
#include <string>
#include <sstream>
#include <utility>
//#include "RooStudentT.h"
#include "fitModels.C"
#include "TH1F.h"
#include "RooBackPdf.h"
#include "RooDataHist.h"
#include "TH1F.h"
#include "TLine.h"
#include "TGraphErrors.h"
#include "RooFFTConvPdf.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "X3872.h"
using namespace RooFit;
using namespace RooStats;

 // D0 & D*0
  const double m_D0     = 1864.83  * 1e-3;
  const double m_D0star = 2006.85  * 1e-3 ;




void plotHan(std::string canname  = "can.eps", std::string blurb = "" ,double gval = 0.3, double Egval = -11e-3, double rhoVal = 0.007, double omegaVal = 0.036, double g0val = 0, double minmass = 3850e-3, double maxmass =3890e-3 ){


 RooRealVar m("m","m",minmass, maxmass);
  
 RooRealVar* mx =new RooRealVar("mx", "mx", Egval + m_D0  + m_D0star ,  0.,1 );
 RooRealVar* g0 =new RooRealVar("g0", "g0", g0val,  0.,1 ); g0->setConstant(true);
 RooRealVar* g =new RooRealVar("g", "g", gval,  0., 1); g->setConstant(true);
 RooRealVar* f_rho =new RooRealVar("f_rho", "f_rho", rhoVal,  0., 5);f_rho->setConstant(true);
 RooRealVar* f_omega =new RooRealVar("f_omega", "f_omega", omegaVal,  0., 5); f_omega->setConstant(true);

 X3872::Hanhard* h1 = hanhard(m,minmass, maxmass, mx,g,f_rho,f_omega, g0, "hanhard1");
 
 std::cout << h1->R_rho() << std::endl;
 std::cout << h1->R_omega() << std::endl;
  
  TCanvas* can = new TCanvas("can","can", 800, 600);
  RooPlot* mframe = m.frame();

  h1->plotOn(mframe, LineColor(2));

   mframe->SetTitle("");
   TAxis* xachse = mframe->GetXaxis(); TAxis* yachse = mframe->GetYaxis();
   xachse->SetTitleFont (132);
   yachse->SetTitleFont (132);
   xachse->SetLabelFont (132);
   yachse->SetLabelFont (132); 
   xachse->SetTitle("m(J/#psi #pi^{+} #pi^{#font[122]{ -} }) [GeV/c^{2}]"); 
   std::string aystring = massCandidates(minmass,maxmass,50); 
   yachse->SetTitle(aystring.c_str()); 
   yachse->SetTitleOffset(0.8); 
   xachse->SetTitleOffset(0.9);
   xachse->SetTitleSize(0.05); 
   yachse->SetTitleSize(0.055);

   mframe->Draw();

 
   TLatex *   tex2 = new TLatex(3.88, 0.04 ,blurb.c_str());
   tex2->SetTextFont(132);
   tex2->SetLineWidth(2);
   tex2->Draw();

   can->Print(canname.c_str());
   std::string epstopdf = "epstopdf " + canname ;
   gSystem->Exec(epstopdf.c_str());

}

void plotWithScale(std::string canname = "can.eps", std::string blurb = "" ,double lambda = 1, double gval = 0.3, double Egval = -11e-3, double rhoVal = 0.007, double omegaVal = 0.036, double g0val = 0){

  plotHan(canname, blurb, lambda*gval, lambda*Egval, lambda*rhoVal, lambda*omegaVal, lambda*g0val);

  
}

void makeplots(){

  // scan through and make plots
  plotWithScale("lambda_2.eps","#lambda = 0.2",0.2);
  plotWithScale("lambda_4.eps","#lambda = 0.4",0.4); 
  plotWithScale("lambda_6.eps","#lambda = 0.6",0.6); 
  plotWithScale("lambda_8.eps","#lambda = 0.8",0.8); 
  plotWithScale("lambda_10.eps","#lambda = 1",1); 
  plotWithScale("lambda_12.eps","#lambda = 1.2",1.2);

  // Belle set table 1 0907.4901
  plotHan("belle1.eps","Belle 1", 0.3,-12.8e-3 , 0.0077, 0.0407,1.1e-3);
  plotHan("belle2.eps","Belle 2", 0.137,-12.3e-3 ,0.000471, 0.00271,1e-3);
  plotHan("belle3.eps","Belle 3", 0.091,-7.8e-3 , 0.00009, 0.00523,2e-3);

  // Babar set table 2
  plotHan("babar4.eps","Babar 4", 0.225,-9.7e-3 , 0.0065, 0.0360,1e-3);
  plotHan("babar5.eps","Babar 5", 0.145,-6.0e-3 , 0.0040, 0.0230,2e-3);
  plotHan("babar6.eps","Babar 6", 0.08,-8.4e-3 , 0.0002, 0.0010,1e-3);
  plotHan("babar7.eps","Babar 7", 0.09,-9.0e-3 , 0.0005, 0.0029,2e-3);

  
}


