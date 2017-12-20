#include "TROOT.h"
#include "TMath.h"
#include "RooWindowPdf.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
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

using namespace RooFit;

void fitMCpsi2CB(std::string filename="../../data/mcpsi_cut.root", 
                 std::string treename="DecayTree", std::string mbranch="mjpipi",
                 double minmass = 3666., double maxmass = 3706.,
                 double fval  = 0.7    , bool fixf  = false,   
                 int nbins = 10000     , double buff = 0.5){

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

    // *** Initialise RooFit Variables ***
    std::cout << "Initialising variables." << std::endl;
    // Min and Max RooWindowPDF
    RooRealVar* hi  = new RooRealVar("hi" , "hi" , bhi);  
    RooRealVar* low = new RooRealVar("low", "low", blow);
    // Mass and Width to fit  
    RooRealVar* mx    = new RooRealVar("mx"   , "mx"   , mval, mval-20., mval+20.);
    RooRealVar* width = new RooRealVar("width", "width", wval,    0.01, 5*wval); width->setConstant(true);
    // Variables for DCB shape
    RooRealVar* f  = new RooRealVar("f" , "f" , fval, 0., 1.); f->setConstant(fixf);
    RooRealVar* s1 = new RooRealVar("s1", "s1", 2., 1., 5.);
    RooRealVar* s2 = new RooRealVar("s2", "s2", 3., 1., 5.);
    RooRealVar* n1 = new RooRealVar("n1", "n1", 1.); n1->setConstant(true); // We expect this value from QED
    RooRealVar* n2 = new RooRealVar("n2", "n2", 2.5, .5, 10.);
    RooFormulaVar* a1  = new RooFormulaVar("a1", "a1", "2.48*pow(59.5 * s1, 10.5)/(1 + pow(59.5 * s1 , 10.5))", RooArgSet(*s1));
    RooFormulaVar* a2  = new RooFormulaVar("a2", "a2", "2.48*pow(59.5 * s2, 10.5)/(1 + pow(59.5 * s2 , 10.5))", RooArgSet(*s2));
    // Vars for Breit-Wigner
    RooRealVar* ma = new RooRealVar("ma", "ma", 3096.916);
    RooRealVar* mb = new RooRealVar("mb", "mb", 507.);
    RooRealVar* radius = new RooRealVar("radius", "radius", 3e-3);
    RooRealVar* spin = new RooRealVar("spin", "spin", 0);  
    
    // *** Initialise RooFit Models *** //
    // BreitWigner model
    RooRelBreitWigner* bw1 = new RooRelBreitWigner("bw1", "bw1", *m, *mx, *width, *radius, *ma, *mb, *spin);
    RooRelBreitWigner* bw2 = new RooRelBreitWigner("bw2", "bw2", *m, *mx, *width, *radius, *ma, *mb, *spin);
    // Window due to limited MC width
    RooWindowPdf* windowpdf1 = new RooWindowPdf("window1","window1", *mx, *hi, *low);
    RooWindowPdf* windowpdf2 = new RooWindowPdf("window2","window2", *mx, *hi, *low);
    // Windowed Breit-Wigner
    RooProdPdf* wbW1 = new RooProdPdf("wbW1", "wbW1", RooArgList(*bw1, *windowpdf1));
    RooProdPdf* wbW2 = new RooProdPdf("wbW2", "wbW2", RooArgList(*bw2, *windowpdf2));
    // Initialise models for fitting
    RooCBShape* cb1 = new RooCBShape("cb1", "cb1", *m, RooConst(0.), *s1, *a1, *n1);
    RooCBShape* cb2 = new RooCBShape("cb2", "cb2", *m, RooConst(0.), *s2, *a2, *n2);
    
    // *** Convolve bw resonance with resolution function *** //
    std::cout << "Convolving pdfs." << std::endl;
    // Fast fourier transform convolution
    m->setBins(nbins,"cache");  // what is cache?
    // Convolve breit-wigner with models
    RooFFTConvPdf* fun1 = new RooFFTConvPdf("fun1","fun1", *m, *bw1, *cb1);
    RooFFTConvPdf* fun2 = new RooFFTConvPdf("fun2","fun2", *m, *bw2, *cb2);
    // Specify behaviour at boundaries
    fun1->setBufferStrategy(RooFFTConvPdf::Extend);
    fun1->setBufferFraction(buff); 
    fun2->setBufferStrategy(RooFFTConvPdf::Extend);
    fun2->setBufferFraction(buff);

    // *** Add two CB Shapes together ***
    RooAbsPdf* sigmodel = new RooAddPdf("sigmodel", "sigmodel", RooArgList(*fun1, *fun2), *f);

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
    sigmodel->plotOn(mframe, Components("cb1") , LineColor(4), LineStyle(2));
 	sigmodel->plotOn(mframe, Components("cb2") , LineColor(3), LineStyle(3));
    sigmodel->plotOn(mframe, LineColor(2)); 
    // Format plot
    mframe->SetTitle("");
    TAxis* xachse = mframe->GetXaxis(); TAxis* yachse = mframe->GetYaxis();
    xachse->SetTitleFont (132);
    yachse->SetTitleFont (132);
    xachse->SetLabelFont (132);
    yachse->SetLabelFont (132); 
    xachse->SetTitle("m(J/#psi #pi^{+} #pi^{-}) [MeV/c^{2}]"); 
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
    can1->SaveAs("resolution_fit_psi_2CB.pdf");
}


void fitMCx38722CB(std::string filename="../../data/mcx_cut.root", 
                   std::string treename="DecayTree", std::string mbranch="mjpipi",
                   double minmass = 3852., double maxmass = 3892.,
                   double fval  = 0.7    , bool fixf  = false,   
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
    RooRealVar* f  = new RooRealVar("f" , "f" , fval, 0., 1.); f->setConstant(fixf);
    RooRealVar* s1 = new RooRealVar("s1", "s1", 2., 1., 5.);
    RooRealVar* s2 = new RooRealVar("s2", "s2", 3., 1., 5.);
    RooRealVar* n1 = new RooRealVar("n1", "n1", 1.); n1->setConstant(true); // We expect this value from QED
    RooRealVar* n2 = new RooRealVar("n2", "n2", 1.15, .5, 10.); n2->setConstant(true); // Fixing from psi(2S)
    RooFormulaVar* a1  = new RooFormulaVar("a1", "a1", "2.48*pow(59.5 * s1, 10.5)/(1 + pow(59.5 * s1 , 10.5))", RooArgSet(*s1));
    RooFormulaVar* a2  = new RooFormulaVar("a2", "a2", "2.48*pow(59.5 * s2, 10.5)/(1 + pow(59.5 * s2 , 10.5))", RooArgSet(*s2));
    // Vars for Breit-Wigner
    RooRealVar* ma = new RooRealVar("ma", "ma", 3096.916);
    RooRealVar* mb = new RooRealVar("mb", "mb", 507.);
    RooRealVar* radius = new RooRealVar("radius", "radius", 3e-3);
    RooRealVar* spin = new RooRealVar("spin", "spin", 0);  
    
    // *** Initialise RooFit Models *** //
    // BreitWigner model
    RooRelBreitWigner* bw1 = new RooRelBreitWigner("bw1", "bw1", *m, *mx, *width, *radius, *ma, *mb, *spin);
    RooRelBreitWigner* bw2 = new RooRelBreitWigner("bw2", "bw2", *m, *mx, *width, *radius, *ma, *mb, *spin);
    // Window due to limited MC width
    RooWindowPdf* windowpdf1 = new RooWindowPdf("window1","window1", *mx, *hi, *low);
    RooWindowPdf* windowpdf2 = new RooWindowPdf("window2","window2", *mx, *hi, *low);
    // Windowed Breit-Wigner
    RooProdPdf* wbW1 = new RooProdPdf("wbW1", "wbW1", RooArgList(*bw1, *windowpdf1));
    RooProdPdf* wbW2 = new RooProdPdf("wbW2", "wbW2", RooArgList(*bw2, *windowpdf2));
    // Initialise models for fitting
    RooCBShape* cb1 = new RooCBShape("cb1", "cb1", *m, RooConst(0.), *s1, *a1, *n1);
    RooCBShape* cb2 = new RooCBShape("cb2", "cb2", *m, RooConst(0.), *s2, *a2, *n2);
    
    // *** Convolve bw resonance with resolution function *** //
    std::cout << "Convolving pdfs." << std::endl;
    // Fast fourier transform convolution
    m->setBins(nbins,"cache");  // what is cache?
    // Convolve breit-wigner with models
    RooFFTConvPdf* fun1 = new RooFFTConvPdf("fun1","fun1", *m, *bw1, *cb1);
    RooFFTConvPdf* fun2 = new RooFFTConvPdf("fun2","fun2", *m, *bw2, *cb2);
    // Specify behaviour at boundaries
    fun1->setBufferStrategy(RooFFTConvPdf::Extend);
    fun1->setBufferFraction(buff); 
    fun2->setBufferStrategy(RooFFTConvPdf::Extend);
    fun2->setBufferFraction(buff);

    // *** Add two Gaussians together ***
    RooAbsPdf* sigmodel = new RooAddPdf("sigmodel", "sigmodel", RooArgList(*fun1, *fun2), *f);

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
    sigmodel->plotOn(mframe, Components("cb1") , LineColor(4), LineStyle(2));
    sigmodel->plotOn(mframe, Components("cb2") , LineColor(3), LineStyle(3));
    sigmodel->plotOn(mframe, LineColor(2)); 
    // Format plot
    mframe->SetTitle("");
    TAxis* xachse = mframe->GetXaxis(); TAxis* yachse = mframe->GetYaxis();
    xachse->SetTitleFont (132);
    yachse->SetTitleFont (132);
    xachse->SetLabelFont (132);
    yachse->SetLabelFont (132); 
    xachse->SetTitle("m(J/#psi #pi^{+} #pi^{-}) [MeV/c^{2}]"); 
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
    can1->SaveAs("resolution_fit_x_2CB.pdf");
}

void fitData2CB(std::string filename="../../data/data_cut.root", 
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
    // Min and Max RooWindowPDF
    RooRealVar* hi  = new RooRealVar("hi" , "hi" , bhi);  
    RooRealVar* low = new RooRealVar("low", "low", blow);
    // Mass and Width to fit  
    RooRealVar* mx    = new RooRealVar("mx"   , "mx"   , mval, mval-20., mval+20.);
    RooRealVar* width = new RooRealVar("width", "width", wval,    0.01, 5*wval); width->setConstant(true);
    // Variables for DCB shape
    RooRealVar* f  = new RooRealVar("f" , "f" , fval, 0., 1.); f->setConstant(fixf);
    RooRealVar* s1 = new RooRealVar("s1", "s1", 2., 1., 5.);  s1->setConstant(true); // Take from data
    RooRealVar* s2 = new RooRealVar("s2", "s2", 3., 1., 5.);  s2->setConstant(true); // Take from data
    RooRealVar* sc = new RooRealVar("sc", "sc", 1., 0., 2.);  // Scaling factor
    RooFormulaVar* sc_s1 = new RooFormulaVar("sc_s1", "sc_s1", "sc*s1", RooArgSet(*s1, *sc));
    RooFormulaVar* sc_s2 = new RooFormulaVar("sc_s2", "sc_s2", "sc*s2", RooArgSet(*s2, *sc));
    RooRealVar* n1 = new RooRealVar("n1", "n1", 1.);          n1->setConstant(true); // We expect this value from QED
    RooRealVar* n2 = new RooRealVar("n2", "n2", 1.15459e+00); n2->setConstant(true);
    RooFormulaVar* a1  = new RooFormulaVar("a1", "a1", "2.48*pow(59.5 * sc_s1, 10.5)/(1 + pow(59.5 * sc_s1 , 10.5))", RooArgSet(*sc_s1));
    RooFormulaVar* a2  = new RooFormulaVar("a2", "a2", "2.48*pow(59.5 * sc_s2, 10.5)/(1 + pow(59.5 * sc_s2 , 10.5))", RooArgSet(*sc_s2));
    //RooRealVar* a1 = new RooRealVar("a1", "a1", 2., .5, 10.);
    //RooRealVar* a2 = new RooRealVar("a2", "a2", 3., .5, 10.);
    // Variables for exp background
    RooRealVar* c = new RooRealVar("a", "a", 0.0001, -1, 1.);
    // Vars for Breit-Wigner
    RooRealVar* ma = new RooRealVar("ma", "ma", 3096.916);
    RooRealVar* mb = new RooRealVar("mb", "mb", 507.);
    RooRealVar* radius = new RooRealVar("radius", "radius", 3e-3);
    RooRealVar* spin = new RooRealVar("spin", "spin", 0);  
    
    // *** Initialise RooFit Models *** //
    // BreitWigner model
    RooRelBreitWigner* bw1 = new RooRelBreitWigner("bw1", "bw1", *m, *mx, *width, *radius, *ma, *mb, *spin);
    RooRelBreitWigner* bw2 = new RooRelBreitWigner("bw2", "bw2", *m, *mx, *width, *radius, *ma, *mb, *spin);
    // Window due to limited MC width
    RooWindowPdf* windowpdf1 = new RooWindowPdf("window1","window1", *mx, *hi, *low);
    RooWindowPdf* windowpdf2 = new RooWindowPdf("window2","window2", *mx, *hi, *low);
    // Windowed Breit-Wigner
    RooProdPdf* wbW1 = new RooProdPdf("wbW1", "wbW1", RooArgList(*bw1, *windowpdf1));
    RooProdPdf* wbW2 = new RooProdPdf("wbW2", "wbW2", RooArgList(*bw2, *windowpdf2));
    // Initialise models for fitting
    RooCBShape* cb1 = new RooCBShape("cb1", "cb1", *m, RooConst(0.), *sc_s1, *a1, *n1);
    RooCBShape* cb2 = new RooCBShape("cb2", "cb2", *m, RooConst(0.), *sc_s2, *a2, *n2);
    RooExponential* bgr = new RooExponential("bgr", "bgr", *m, *c);
    RooRealVar* sig_yield = new RooRealVar("sig_yield", "sig_yield", 55000., 20000., 70000.);
    RooRealVar* bgr_yield = new RooRealVar("bgr_yield", "bgr_yield", 20000.,     0., 50000.); 
    
    // *** Convolve bw resonance with resolution function *** //
    std::cout << "Convolving pdfs." << std::endl;
    // Fast fourier transform convolution
    m->setBins(nbins,"cache");  // what is cache?
    // Convolve breit-wigner with models
    RooFFTConvPdf* fun1 = new RooFFTConvPdf("fun1","fun1", *m, *bw1, *cb1);
    RooFFTConvPdf* fun2 = new RooFFTConvPdf("fun2","fun2", *m, *bw2, *cb2);
    fun1->setBufferStrategy(RooFFTConvPdf::Extend);
    fun1->setBufferFraction(buff); 
    fun2->setBufferStrategy(RooFFTConvPdf::Extend);
    fun2->setBufferFraction(buff);

    // *** Add two Gaussians together ***
    RooAbsPdf* sigmodel = new RooAddPdf("sigmodel", "sigmodel", RooArgList(*fun1, *fun2), *f);
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
    fitmodel->plotOn(mframe, Components("bgr"), LineColor(5), LineStyle(4)); 
    sigmodel->plotOn(mframe, Components("cb1"), LineColor(4), LineStyle(2));
    sigmodel->plotOn(mframe, Components("cb2"), LineColor(3), LineStyle(3));
    fitmodel->plotOn(mframe, LineColor(2));
    // Format plot
    mframe->SetTitle("");
    TAxis* xachse = mframe->GetXaxis(); TAxis* yachse = mframe->GetYaxis();
    xachse->SetTitleFont (132);
    yachse->SetTitleFont (132);
    xachse->SetLabelFont (132);
    yachse->SetLabelFont (132); 
    xachse->SetTitle("m(J/#psi #pi^{+} #pi^{-}) [MeV/c^{2}]"); 
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
    can1->SaveAs("resolution_fit_data_2CB.pdf");
}


void fitX38232CB(std::string filename="../../data/data_cut.root", 
                 std::string treename="DecayTree", std::string mbranch="mjpipi",
                 double minmass = 3666., double maxmass = 3706.,   // ??? find defined mass range
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
    double wval   = 0.;
    double mval   = 3823.;    // ??? confirm from PDG
    double avalue = 518.;
    double blow   = 3685.25;  // ???
    double bhi    = 3690.65;  // ???

    // *** Initialise RooFit Variables ***
    std::cout << "Initialising variables." << std::endl;
    // Min and Max RooWindowPDF
    RooRealVar* hi  = new RooRealVar("hi" , "hi" , bhi);  
    RooRealVar* low = new RooRealVar("low", "low", blow);
    // Mass and Width to fit  
    RooRealVar* mx    = new RooRealVar("mx"   , "mx"   , mval, mval-20., mval+20.); mx->setConstant(false);
    RooRealVar* width = new RooRealVar("width", "width", wval,    0., 5*wval); width->setConstant(false);
    // Variables for DCB shape - take from MC and interpolation
    RooRealVar* f  = new RooRealVar("f" , "f" , 0.7, 0., 1.); f ->setConstant(true);
    RooRealVar* sc_s1 = new RooRealVar("sc_s1", "sc_s1", 2., 1., 5.);  sc_s1->setConstant(true); // Take from data
    RooRealVar* sc_s2 = new RooRealVar("sc_s2", "sc_s2", 3., 1., 5.);  sc_s2->setConstant(true); // Take from data
    RooRealVar* n1 = new RooRealVar("n1", "n1", 1.);          n1->setConstant(true); // We expect this value from QED
    RooRealVar* n2 = new RooRealVar("n2", "n2", 1.15459e+00); n2->setConstant(true);
    RooRealVar* a1 = new RooRealVar("a1", "a1", 2., .5, 10.);
    RooRealVar* a2 = new RooRealVar("a2", "a2", 3., .5, 10.);
    // Variables for exp background
    RooRealVar* c = new RooRealVar("a", "a", 0.0001, -1, 1.);
    // Vars for Breit-Wigner
    RooRealVar* ma = new RooRealVar("ma", "ma", 3096.916);
    RooRealVar* mb = new RooRealVar("mb", "mb", 507.);
    RooRealVar* radius = new RooRealVar("radius", "radius", 3e-3);
    RooRealVar* spin = new RooRealVar("spin", "spin", 0);  
    
    // *** Initialise RooFit Models *** //
    // BreitWigner model
    RooRelBreitWigner* bw1 = new RooRelBreitWigner("bw1", "bw1", *m, *mx, *width, *radius, *ma, *mb, *spin);
    RooRelBreitWigner* bw2 = new RooRelBreitWigner("bw2", "bw2", *m, *mx, *width, *radius, *ma, *mb, *spin);
    // Window due to limited MC width
    RooWindowPdf* windowpdf1 = new RooWindowPdf("window1","window1", *mx, *hi, *low);
    RooWindowPdf* windowpdf2 = new RooWindowPdf("window2","window2", *mx, *hi, *low);
    // Windowed Breit-Wigner
    RooProdPdf* wbW1 = new RooProdPdf("wbW1", "wbW1", RooArgList(*bw1, *windowpdf1));
    RooProdPdf* wbW2 = new RooProdPdf("wbW2", "wbW2", RooArgList(*bw2, *windowpdf2));
    // Initialise models for fitting
    RooCBShape* cb1 = new RooCBShape("cb1", "cb1", *m, RooConst(0.), *sc_s1, *a1, *n1);
    RooCBShape* cb2 = new RooCBShape("cb2", "cb2", *m, RooConst(0.), *sc_s2, *a2, *n2);
    RooExponential* bgrmodel = new RooExponential("bgrmodel", "bgrmodel", *m, *c);
    RooRealVar* sig_yield = new RooRealVar("sig_yield", "sig_yield", 55000., 20000., 70000.); // ??? look at data to find sensible values
    RooRealVar* bgr_yield = new RooRealVar("bgr_yield", "bgr_yield", 20000.,     0., 50000.); 
    
    // *** Convolve bw resonance with resolution function *** //
    std::cout << "Convolving pdfs." << std::endl;
    // Fast fourier transform convolution
    m->setBins(nbins,"cache");  // what is cache?
    // Convolve breit-wigner with models
    RooFFTConvPdf* fun1 = new RooFFTConvPdf("fun1","fun1", *m, *bw1, *cb1);
    RooFFTConvPdf* fun2 = new RooFFTConvPdf("fun2","fun2", *m, *bw2, *cb2);
    fun1->setBufferStrategy(RooFFTConvPdf::Extend);
    fun1->setBufferFraction(buff); 
    fun2->setBufferStrategy(RooFFTConvPdf::Extend);
    fun2->setBufferFraction(buff);

    // *** Add two Gaussians together ***
    RooAbsPdf* sigmodel = new RooAddPdf("sigmodel", "sigmodel", RooArgList(*fun1, *fun2), *f);
    // *** Add signal and background
    RooAddPdf* fitmodel = new RooAddPdf("fitmodel", "fitmodel", RooArgList(*sigmodel, *bgrmodel), RooArgList(*sig_yield, *bgr_yield));

    // *** Fit and plot ***
    RooFitResult* sig_result = fitmodel->fitTo(*data, Save(), RooFit::Optimize(kFALSE)); 
    RooFitResult* bgr_result = bgrmodel->fitTo(*data, Save(), RooFit::Optimize(kFALSE));
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
    bgrmodel->plotOn(mframe, LineColor(3));
    fitmodel->plotOn(mframe, LineColor(2));
    // Format plot
    mframe->SetTitle("");
    TAxis* xachse = mframe->GetXaxis(); TAxis* yachse = mframe->GetYaxis();
    xachse->SetTitleFont (132);
    yachse->SetTitleFont (132);
    xachse->SetLabelFont (132);
    yachse->SetLabelFont (132); 
    xachse->SetTitle("m(J/#psi #pi^{+} #pi^{-}) [MeV/c^{2}]"); 
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
    can1->SaveAs("resolution_fit_X3823_2CB.pdf");

    // ??? Determine significance from Greig's formula
}