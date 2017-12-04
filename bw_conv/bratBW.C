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
#include "X3872_constants.h"
//#include "RooBackPdf.h"
#include "RooDataHist.h"
#include "TH1F.h"
#include "TLine.h"
#include "TGraphErrors.h"
#include "RooFFTConvPdf.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "Braatan.h"
using namespace RooFit;
using namespace RooStats;

 // D0 & D*0
//  const double m_D0     = 1864.83  * 1e-3;
//  const double m_D0star = 2006.85  * 1e-3 ;

using namespace X3872;



void bratBW(std::string canname  = "can.eps", std::string blurb = "" , double minmass = 3840e-3, double maxmass =3891e-3 ){


 RooRealVar m("m","m",minmass, maxmass);


 RooRealVar* mx =new RooRealVar("mx", "mx", 3.8717 ,  0.,5 );
 RooRealVar* width =new RooRealVar("width", "width", 1.4e-3, 0, 1);

 RooRelBreitWigner* bWigner = bwg2(m,minmass,maxmass,mx,width,"bw", 0.717, true);

 RooDataSet* dset = bWigner->generate(m,1e5);
 
 RooRealVar* gRe =new RooRealVar("gRe", "gRe", 38.4e-3 ,  -0.01,1 ); gRe->setConstant(true);
 RooRealVar* gIm =new RooRealVar("gIm", "gIm", 3e-3, 0, 1);

 RooRealVar* thres = new RooRealVar("thres","thres", m_D0  + m_D0star, 0, 5); //thres->setConstant(true);
 X3872::Braatan* brat = new X3872::Braatan("brat","brat" , m, *gRe, *gIm , *thres); 
 
 brat->fitTo(*dset); 
 
  TCanvas* can = new TCanvas("can","can", 800, 600);
  RooPlot* mframe = m.frame();

  dset->plotOn(mframe);
  brat->plotOn(mframe, LineColor(2));

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

   // can->Print(canname.c_str());
   // std::string epstopdf = "epstopdf " + canname ;
   // gSystem->Exec(epstopdf.c_str());

   //std::cout << brat->EX() << " " << brat->GammaX() << std::endl;
   //std::cout << brat->fwhm() << std::endl;
}


