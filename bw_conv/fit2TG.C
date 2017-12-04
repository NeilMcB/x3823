#include "TROOT.h"
#include "TMath.h"
#include "RooWindowPdf.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooRelBreitWigner.h"
#include "TSystem.h"
#include "TFile.h"
#include "RooProdPdf.h"
#include "RooFFTConvPdf.h"
#include "RooAddPdf.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TPad.h"
#include "TAxis.h"
#include "TLine.h"
#include "RooExponential.h"
#include "RooHist.h"
#include "RooFitResult.h"
#include "RooNumConvPdf.h"
#include "RooConstVar.h"
#include "RooDCBShape.h"

using namespace RooFit;

void fitMCpsi2TG(std::string filename="../data/mcpsi_Qcut_vetos.root", 
                   std::string treename="DecayTree", std::string mbranch="mjpipi",
                   double minmass = 3666., double maxmass = 3706.,
                   double fval  = 0.7    , bool fixf  = true,   
                   int nbins = 10000, double buff = 0.5){

    // Load missing packages
    gSystem->Load("./RooRelBreitWigner_cc.so");
    gSystem->Load("./RooWindowPdf_cxx.so");
    gSystem->Load("./RooDCBShape_cxx.so");

    // Load data
    std::cout << "Loading tree " << treename << " from " << filename << std::endl;
    TFile* theFile = new TFile(filename.c_str()); 
    TTree* tree = (TTree*)theFile->Get(treename.c_str());

    // Prepare variable for fitting to
    std::cout << "Filling dataset." << std::endl;
    RooRealVar* m = new RooRealVar(mbranch.c_str(), mbranch.c_str(), minmass, maxmass);
    // Initialise and fill dataset
    RooDataSet* data = new RooDataSet("data", "data", RooArgSet(*m), Import(*tree)); 

    // Initialise Breit-Wigner variables
    double wval   = 0.304;
    double mval   = 3686.11;
	double avalue = 518.;
	double blow   = 3685.25;
	double bhi    = 3690.65;
    //std::string par1 = "2.48";
    //std::string par2 = "10.5";
    //std::string par3 = "0.595*1000";

    // *** Initialise RooFit Variables ***
    std::cout << "Initialising variables." << std::endl;
    // Min and Max RooWindowPDF
    RooRealVar* hi  = new RooRealVar("hi" , "hi" , bhi);  
    RooRealVar* low = new RooRealVar("low", "low", blow);
    // Mass and Width to fit  
    RooRealVar* mx    = new RooRealVar("mx"   , "mx"   , mval, mval-20., mval+20.);
    RooRealVar* width = new RooRealVar("width", "width", wval,    0.01, 5*wval); width->setConstant(true);
    // Variables for DCB shape
    RooRealVar* sig = new RooRealVar("sig", "sig", 2., 1., 5.);
    RooRealVar* nL  = new RooRealVar("nL", "nL", 1.); nL->setConstant(true); // We expect this value from QED
    RooRealVar* nR  = new RooRealVar("nR", "nR", 2.5, .5, 20.);
    RooRealVar* aL  = new RooRealVar("aL", "aL", 2., .5, 10.);
    RooRealVar* aR  = new RooRealVar("aR", "aR", 3., .5, 10.);
    // Vars for Breit-Wigner
    RooRealVar* ma = new RooRealVar("ma", "ma", 3096.916);
    RooRealVar* mb = new RooRealVar("mb", "mb", 507.);
    RooRealVar* radius = new RooRealVar("radius", "radius", 3e-3);
    RooRealVar* spin = new RooRealVar("spin", "spin", 0);  
    
    // *** Initialise RooFit Models *** //
    // BreitWigner model
    RooRelBreitWigner* bw = new RooRelBreitWigner("bw", "bw", *m, *mx, *width, *radius, *ma, *mb, *spin);
    // Window due to limited MC width
    RooWindowPdf* windowpdf = new RooWindowPdf("window","window", *mx, *hi, *low);
    // Windowed Breit-Wigner
    RooProdPdf* wbW1 = new RooProdPdf("wbW", "wbW", RooArgList(*bw, *windowpdf));
    // Initialise models for fitting
    RooDCBShape* dcb = new RooDCBShape("dcb", "dcb", *m, RooConst(0.), *sig, *aL, *aR, *nL, *nR);
    
    // *** Convolve bw resonance with resolution function *** //
    std::cout << "Convolving pdfs." << std::endl;
    // Fast fourier transform convolution
    m->setBins(nbins,"cache");  // what is cache?
    // Convolve breit-wigner with models
    RooFFTConvPdf* sigmodel = new RooFFTConvPdf("sigmodel","sigmodel", *m, *bw, *dcb);
    sigmodel->setBufferStrategy(RooFFTConvPdf::Extend);
    sigmodel->setBufferFraction(buff); 

    // *** Fit and plot ***
    RooFitResult* fresult = sigmodel->fitTo(*data, Save(), RooFit::Optimize(kFALSE)); 
    // Prepare canvas
    RooPlot* mframe = m->frame();
    TCanvas* can1 = new TCanvas("can1", "can1", 400, 500);
    TPad* pad1 = new TPad("pad1","This is pad1",.01,.01,.95,.77);
    TPad* pad2 = new TPad("pad2","This is pad2",.01,.76,.95,.97);
    pad1->Draw();
    pad2->Draw();
    pad1->cd();
    // Plot data and model
    data->plotOn(mframe, Binning(50));
    sigmodel->plotOn(mframe, LineColor(2)); 
    // Format plot
    mframe->SetTitle("");
    TAxis* xachse = mframe->GetXaxis(); TAxis* yachse = mframe->GetYaxis();
    xachse->SetTitleFont (132);
    yachse->SetTitleFont (132);
    xachse->SetLabelFont (132);
    yachse->SetLabelFont (132); 
    xachse->SetTitle("m(J/#psi #pi^{+} #pi^{-}) [GeV/c^{2}]"); 
    yachse->SetTitle("candidates"); 
    yachse->SetTitleOffset(1.05); 
    xachse->SetTitleOffset(0.95);
    xachse->SetTitleSize(0.05); 
    yachse->SetTitleSize(0.05);
    // Draw
    mframe->Draw();
    // Plot pulls
    pad2->cd();
    RooPlot* frame2 =  m->frame();
    RooHist* phist = mframe->pullHist( );
    frame2->addPlotable(phist,"P");
    frame2->SetMinimum(-6.5); 
    frame2->SetMaximum(6.5);
    // Make pretty
    xachse = frame2->GetXaxis(); yachse = frame2->GetYaxis();
    xachse->SetTitleFont (132);
    yachse->SetTitleFont (132);
    xachse->SetLabelFont (132);
    yachse->SetLabelFont (132); 
    yachse->SetTitleOffset(0.15); 
    yachse->SetTitleSize(0.2);
    yachse->SetLabelSize(0.08);
    xachse->SetLabelOffset(0.2);
    yachse->SetNdivisions(110);
    frame2->SetXTitle("");
    frame2->SetYTitle("pull");
    frame2->SetTitle("");
    frame2->Draw() ; 
    TLine* line = new TLine(minmass,0, maxmass, 0);
    // Draw
    line->Draw();
    can1->SaveAs("resolution_fit_psi_2TG.pdf");
}

void fitMCx38722TG(std::string filename="../data/mcx_Qcut_vetos.root", 
                   std::string treename="DecayTree", std::string mbranch="mjpipi",
                   double minmass = 3852., double maxmass = 3892.,
                   double fval  = 0.7    , bool fixf  = true,   
                   int nbins = 10000, double buff = 0.5){

    // Load missing packages
    gSystem->Load("./RooRelBreitWigner_cc.so");
    gSystem->Load("./RooWindowPdf_cxx.so");
    gSystem->Load("./RooDCBShape_cxx.so");

    // Load data
    std::cout << "Loading tree " << treename << " from " << filename << std::endl;
    TFile* theFile = new TFile(filename.c_str()); 
    TTree* tree = (TTree*)theFile->Get(treename.c_str());

    // Prepare variable for fitting to
    std::cout << "Filling dataset." << std::endl;
    RooRealVar* m = new RooRealVar(mbranch.c_str(), mbranch.c_str(), minmass, maxmass);
    // Initialise and fill dataset
    RooDataSet* data = new RooDataSet("data", "data", RooArgSet(*m), Import(*tree)); 

    // Initialise Breit-Wigner variables
    double wval   = 0.317;
    double mval   = 3871.5;
    double avalue = 717.;
    double blow   = 3869.3;
    double bhi    = 3977.;

    // *** Initialise RooFit Variables ***
    std::cout << "Initialising variables." << std::endl;
    // Min and Max RooWindowPDF
    RooRealVar* hi  = new RooRealVar("hi" , "hi" , bhi);  
    RooRealVar* low = new RooRealVar("low", "low", blow);
    // Mass and Width to fit  
    RooRealVar* mx    = new RooRealVar("mx"   , "mx"   , mval, mval-20., mval+20.);
    RooRealVar* width = new RooRealVar("width", "width", wval,    0.01, 5*wval); width->setConstant(true);
    // Variables for DCB shape
    RooRealVar* sig = new RooRealVar("sig", "sig", 2., 1., 5.);
    RooRealVar* nL  = new RooRealVar("nL", "nL", 1.); nL->setConstant(true); // We expect this value from QED
    RooRealVar* nR  = new RooRealVar("nR", "nR", 2.5, .5, 20.);
    RooRealVar* aL  = new RooRealVar("aL", "aL", 2., .5, 10.); 
    RooRealVar* aR  = new RooRealVar("aR", "aR", 3., .5, 10.); 
    // Vars for Breit-Wigner
    RooRealVar* ma = new RooRealVar("ma", "ma", 3096.916);
    RooRealVar* mb = new RooRealVar("mb", "mb", 507.);
    RooRealVar* radius = new RooRealVar("radius", "radius", 3e-3);
    RooRealVar* spin = new RooRealVar("spin", "spin", 0);  
    
    // *** Initialise RooFit Models *** //
    // BreitWigner model
    RooRelBreitWigner* bw = new RooRelBreitWigner("bw", "bw", *m, *mx, *width, *radius, *ma, *mb, *spin);
    // Window due to limited MC width
    RooWindowPdf* windowpdf = new RooWindowPdf("window","window", *mx, *hi, *low);
    // Windowed Breit-Wigner
    RooProdPdf* wbW1 = new RooProdPdf("wbW", "wbW", RooArgList(*bw, *windowpdf));
    // Initialise models for fitting
    RooDCBShape* dcb = new RooDCBShape("dcb", "dcb", *m, RooConst(0.), *sig, *aL, *aR, *nL, *nR);
    
    // *** Convolve bw resonance with resolution function *** //
    std::cout << "Convolving pdfs." << std::endl;
    // Fast fourier transform convolution
    m->setBins(nbins,"cache");  // what is cache?
    // Convolve breit-wigner with models
    RooFFTConvPdf* sigmodel = new RooFFTConvPdf("sigmodel","sigmodel", *m, *bw, *dcb);
    sigmodel->setBufferStrategy(RooFFTConvPdf::Extend);
    sigmodel->setBufferFraction(buff); 

    // *** Fit and plot ***
    RooFitResult* fresult = sigmodel->fitTo(*data, Save(), RooFit::Optimize(kFALSE)); 
    // Prepare canvas
    RooPlot* mframe = m->frame();
    TCanvas* can1 = new TCanvas("can1", "can1", 400, 500);
    TPad* pad1 = new TPad("pad1","This is pad1",.01,.01,.95,.77);
    TPad* pad2 = new TPad("pad2","This is pad2",.01,.76,.95,.97);
    pad1->Draw();
    pad2->Draw();
    pad1->cd();
    // Plot data and model
    data->plotOn(mframe, Binning(50));
    sigmodel->plotOn(mframe, LineColor(2)); 
    // Format plot
    mframe->SetTitle("");
    TAxis* xachse = mframe->GetXaxis(); TAxis* yachse = mframe->GetYaxis();
    xachse->SetTitleFont (132);
    yachse->SetTitleFont (132);
    xachse->SetLabelFont (132);
    yachse->SetLabelFont (132); 
    xachse->SetTitle("m(J/#psi #pi^{+} #pi^{-}) [GeV/c^{2}]"); 
    yachse->SetTitle("candidates"); 
    yachse->SetTitleOffset(1.05); 
    xachse->SetTitleOffset(0.95);
    xachse->SetTitleSize(0.05); 
    yachse->SetTitleSize(0.05);
    // Draw
    mframe->Draw();
    // Plot pulls
    pad2->cd();
    RooPlot* frame2 =  m->frame();
    RooHist* phist = mframe->pullHist( );
    frame2->addPlotable(phist,"P");
    frame2->SetMinimum(-6.5); 
    frame2->SetMaximum(6.5);
    // Make pretty
    xachse = frame2->GetXaxis(); yachse = frame2->GetYaxis();
    xachse->SetTitleFont (132);
    yachse->SetTitleFont (132);
    xachse->SetLabelFont (132);
    yachse->SetLabelFont (132); 
    yachse->SetTitleOffset(0.15); 
    yachse->SetTitleSize(0.2);
    yachse->SetLabelSize(0.08);
    xachse->SetLabelOffset(0.2);
    yachse->SetNdivisions(110);
    frame2->SetXTitle("");
    frame2->SetYTitle("pull");
    frame2->SetTitle("");
    frame2->Draw() ; 
    TLine* line = new TLine(minmass,0, maxmass, 0);
    // Draw
    line->Draw();
    can1->SaveAs("resolution_fit_x_2TG.pdf");
}

void fitData2CB(std::string filename="../data/data_Qcut_vetos.root", 
                std::string treename="DecayTree", std::string mbranch="mjpipi",
                double minmass = 3666., double maxmass = 3706.,
                double fval  = 0.7    , bool fixf  = true,   
                int nbins = 10000, double buff = 0.5){

    // Load missing packages
    gSystem->Load("./RooRelBreitWigner_cc.so");
    gSystem->Load("./RooWindowPdf_cxx.so");
    gSystem->Load("./RooIpatia2_cxx.so");

    // Load data
    std::cout << "Loading tree " << treename << " from " << filename << std::endl;
    TFile* theFile = new TFile(filename.c_str()); 
    TTree* tree = (TTree*)theFile->Get(treename.c_str());

    // Prepare variable for fitting to
    std::cout << "Filling dataset." << std::endl;
    RooRealVar* m = new RooRealVar(mbranch.c_str(), mbranch.c_str(), minmass, maxmass);
    // Initialise and fill dataset
    RooDataSet* data = new RooDataSet("data", "data", RooArgSet(*m), Import(*tree)); 

    // Initialise Breit-Wigner variables
    double wval   = 0.304;
    double mval   = 3686.11;
    double avalue = 518.;
    double blow   = 3685.25;
    double bhi    = 3690.65;
    //std::string par1 = "2.48";
    //std::string par2 = "10.5";
    //std::string par3 = "0.595*1000";

    // *** Initialise RooFit Variables ***
    std::cout << "Initialising variables." << std::endl;
    RooRealVar* hi  = new RooRealVar("hi" , "hi" , bhi);  
    RooRealVar* low = new RooRealVar("low", "low", blow);
    // Mass and Width to fit  
    RooRealVar* mx    = new RooRealVar("mx"   , "mx"   , mval, mval-20., mval+20.);
    RooRealVar* width = new RooRealVar("width", "width", wval,    0.01, 5*wval); width->setConstant(true);
    // Variables for DCB shape
    RooRealVar* sig = new RooRealVar("sig", "sig", 2., 1., 5.);
    RooRealVar* nL  = new RooRealVar("nL", "nL", 1.); nL->setConstant(true); // We expect this value from QED
    RooRealVar* nR  = new RooRealVar("nR", "nR", 2.00000e+01); nR->setConstant(true);
    RooRealVar* aL  = new RooRealVar("aL", "aL", 2.36917e+00); aL->setConstant(true);
    RooRealVar* aR  = new RooRealVar("aR", "aR", 1.74857e+00); aR->setConstant(true);
    // Vars for Breit-Wigner
    RooRealVar* ma = new RooRealVar("ma", "ma", 3096.916);
    RooRealVar* mb = new RooRealVar("mb", "mb", 507.);
    RooRealVar* radius = new RooRealVar("radius", "radius", 3e-3);
    RooRealVar* spin = new RooRealVar("spin", "spin", 0);  
    // Variables for exp background
    RooRealVar* c = new RooRealVar("a", "a", 0.0001, -1, 1.);

    
    // *** Initialise RooFit Models *** //
    // BreitWigner model
    RooRelBreitWigner* bw = new RooRelBreitWigner("bw", "bw", *m, *mx, *width, *radius, *ma, *mb, *spin);
    // Window due to limited MC width
    RooWindowPdf* windowpdf = new RooWindowPdf("window","window", *mx, *hi, *low);
    // Windowed Breit-Wigner
    RooProdPdf* wbW = new RooProdPdf("wbW", "wbW", RooArgList(*bw, *windowpdf));
    // Initialise models for fitting
    RooDCBShape* dcb = new RooDCBShape("dcb", "dcb", *m, RooConst(0.), *sig, *aL, *aR, *nL, *nR);
    RooExponential* bgr = new RooExponential("bgr", "bgr", *m, *c);
    RooRealVar* sig_yield = new RooRealVar("sig_yield", "sig_yield", 55000., 20000., 70000.);
    RooRealVar* bgr_yield = new RooRealVar("bgr_yield", "bgr_yield", 20000.,     0., 50000.); 
    
    // *** Convolve bw resonance with resolution function *** //
    std::cout << "Convolving pdfs." << std::endl;
    // Fast fourier transform convolution
    m->setBins(nbins,"cache");  // what is cache?
    // Convolve breit-wigner with models
    RooFFTConvPdf* sigmodel = new RooFFTConvPdf("sigmodel","sigmodel", *m, *bw, *dcb);
    sigmodel->setBufferStrategy(RooFFTConvPdf::Extend);
    sigmodel->setBufferFraction(buff); 

    // *** Add signal and background
    RooAddPdf* fitmodel = new RooAddPdf("fitmodel", "fitmodel", RooArgList(*sigmodel, *bgr), RooArgList(*sig_yield, *bgr_yield));

    // *** Fit and plot ***
    RooFitResult* fresult = fitmodel->fitTo(*data, Save(), RooFit::Optimize(kFALSE)); 
    // Prepare canvas
    RooPlot* mframe = m->frame();
    TCanvas* can1 = new TCanvas("can1", "can1", 400, 500);
    TPad* pad1 = new TPad("pad1","This is pad1",.01,.01,.95,.77);
    TPad* pad2 = new TPad("pad2","This is pad2",.01,.76,.95,.97);
    pad1->Draw();
    pad2->Draw();
    pad1->cd();
    // Plot data and model
    data->plotOn(mframe, Binning(50));
    fitmodel->plotOn(mframe, LineColor(2));
    fitmodel->plotOn(mframe, Components("bgr"), LineColor(4), LineStyle(2)); 
    fitmodel->plotOn(mframe, LineColor(2));
    // Format plot
    mframe->SetTitle("");
    TAxis* xachse = mframe->GetXaxis(); TAxis* yachse = mframe->GetYaxis();
    xachse->SetTitleFont (132);
    yachse->SetTitleFont (132);
    xachse->SetLabelFont (132);
    yachse->SetLabelFont (132); 
    xachse->SetTitle("m(J/#psi #pi^{+} #pi^{-}) [GeV/c^{2}]"); 
    yachse->SetTitle("candidates"); 
    yachse->SetTitleOffset(1.05); 
    xachse->SetTitleOffset(0.95);
    xachse->SetTitleSize(0.05); 
    yachse->SetTitleSize(0.05);
    // Draw
    mframe->Draw();
    // Plot pulls
    pad2->cd();
    RooPlot* frame2 =  m->frame();
    RooHist* phist = mframe->pullHist( );
    frame2->addPlotable(phist,"P");
    frame2->SetMinimum(-6.5); 
    frame2->SetMaximum(6.5);
    // Make pretty
    xachse = frame2->GetXaxis(); yachse = frame2->GetYaxis();
    xachse->SetTitleFont (132);
    yachse->SetTitleFont (132);
    xachse->SetLabelFont (132);
    yachse->SetLabelFont (132); 
    yachse->SetTitleOffset(0.15); 
    yachse->SetTitleSize(0.2);
    yachse->SetLabelSize(0.08);
    xachse->SetLabelOffset(0.2);
    yachse->SetNdivisions(110);
    frame2->SetXTitle("");
    frame2->SetYTitle("pull");
    frame2->SetTitle("");
    frame2->Draw() ; 
    TLine* line = new TLine(minmass,0, maxmass, 0);
    // Draw
    line->Draw();
    can1->SaveAs("resolution_fit_data_2TG.pdf");
}