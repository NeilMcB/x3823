#include "TROOT.h"
#include "TMath.h"
#include "RooWindowPdf.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooRelBreitWigner.h"
#include "TSystem.h"

using namespace RooFit;
using namespace RooStats;

// Main Function
void fitMCres(std::string fittype="2CB", std::string resname="psi2S", 
                 std::string filename="../data/mcpsi_Qcut_vetos.root", 
                 std::string treename="DecayTree", std::string mbranch="mjpipi",
                 double minmass = 3666., double maxmass = 3706.,
                 double fval  = 0.9    , bool fixf  = false,  
                 double n1val = 2.5    , bool fixN1 = false,
                 double n2val = 1.     , bool fixN2 = false, 
                 bool fixWidth = false,
                 int nbins = 10000, double buff = 0.5){
    
    // Load missing packages
    gSystem->Load("./RooRelBreitWigner_cc.so");
    gSystem->Load("./RooWindowPdf_cxx.so");

    // Load data
    std::cout << "Loading tree " << treename << " from " << filename << std::endl;
    TFile* theFile = new TFile(filename.c_str()); 
    TTree* tree = (TTree*)theFile->Get(treename.c_str());

    // Prepare variable for fitting to
    std::cout << "Filling dataset." << std::endl;
    RooRealVar* m = new RooRealVar(mbranch.c_str(), mbranch.c_str(), minmass, maxmass);
    // Initialise and fill dataset
    RooDataSet* data = new RooDataSet("data", "data", RooArgSet(*m), Import(*tree)); 
    

    // Initialise variables for X(3872) or psi(2S)
    std::cout << "Initialising variables." << std::endl;
    double mval ; double wval; double avalue; double par1, par2,par3;
    double bhi, blow;
    if (resname=="X"){
        wval = 0.317;
        mval = 3871.5;
        avalue = 717.; 
        bhi = 3877.;
        blow = 3869.3;
        par1 = 2.48;
        par2 = 8.1;
        par3 = 0.474*1000;
    }
    else {
        wval = 0.304;
        mval = 3686.11;
        avalue = 518.;
        blow = 3685.25;
        bhi = 3690.65;
        par1 = 2.48;
        par2 = 10.5;
        par3 = 0.595*1000;
    }
    // *** Initialise RooFit Variables ***
    // Min and Max RooWindowPDF
    RooRealVar* hi  = new RooRealVar("hi" , "hi" , bhi);  
    RooRealVar* low = new RooRealVar("low", "low", blow);
    // Mass and Width to fit  
    RooRealVar* mx    = new RooRealVar("mx"   , "mx"   , mval, mval-2., mval+2.);
    RooRealVar* width = new RooRealVar("width", "width", wval,    0.01, 5*wval); width->setConstant(true);
    
    // Sigma for Gauss shape
    RooRealVar* sg = new RooRealVar("sig_g", "sig_g", 0.1, 0.001, 3.);
    // Sigmas for 2CB Shape
    RooRealVar* s1  = new RooRealVar("sig_1", "sig_1", 0.003, 0.001,  3.);
    //RooRealVar* r12 = new RooRealVar("r_s12", "r_s12", 1.   , 0.001, 10.); 
    RooRealVar* s2 = new RooRealVar("sig_2", "sig_2", 4.); s2->setConstant(true);
    // ns for CB    
    RooRealVar* n1 = new RooRealVar("n1","n1",n1val, 0.5, 500); n1->setConstant(fixN1);
    RooRealVar* n2 = new RooRealVar("n2","n2",n2val, 0.5, 200); n2->setConstant(fixN2); 
    // alphas for CB
    RooRealVar* ac  = new RooRealVar("ac" , "ac" ,  5.,    .1, 10. );
    RooRealVar* ac2 = new RooRealVar("ac2", "ac2", -5., -10. ,  -.1);
    // Vars for Breit-Wigner
    RooRealVar* ma = new RooRealVar("ma", "ma", 3096.916);
    RooRealVar* mb = new RooRealVar("mb", "mb", 507.);
    RooRealVar* radius = new RooRealVar("radius", "radius", 3e-3);
    RooRealVar* spin = new RooRealVar("spin", "spin", 0);  
    // BreitWigner models
    RooRelBreitWigner* bw1 = new RooRelBreitWigner("bw1", "bw1", *m, *mx, *width, *radius, *ma, *mb, *spin);
    RooRelBreitWigner* bw2 = new RooRelBreitWigner("bw2", "bw2", *m, *mx, *width, *radius, *ma, *mb, *spin);
    // Window due to limited MC width
    RooWindowPdf* windowpdf = new RooWindowPdf("window","window", *mx, *hi, *low);
    // Windowed Breit-Wigner
    std::cout << "My code reaches here." << std::endl;
    RooProdPdf* wbWigner1 = new RooProdPdf("wb1", "wb1", RooArgList(*bw1, *windowpdf));
    std::cout << "But not here :(" << std::endl;
    RooProdPdf* wbWigner2 = new RooProdPdf("wb2", "wb2", RooArgList(*bw2, *windowpdf));
    // Fraction for combining convolved PDFs
    RooRealVar* f = new RooRealVar("f","f",fval, 0.0,1.0); f->setConstant(fixf);
    // Resolution value for fitting
    RooRealVar* meanres = new RooRealVar("meanres","meanres", 0);
    // Initialise models for fitting
    RooCBShape* cb1 = new RooCBShape("crystalball1", "cb(x,mean,sigma)", *m, *meanres, *s1, *ac , *n1);
    RooCBShape* cb2 = new RooCBShape("crystalball2", "cb(x,mean,sigma)", *m, *meanres, *s2, *ac2, *n2); // should this be ac2?
    RooGaussian* gauss2 = new RooGaussian("gauss2", "gauss", *m, *meanres, *sg);

    std::cout << "Convolving pdfs." << std::endl;

    // Fast fourier transform convolution
    m->setBins(nbins,"cache");  // what is cache?
    
    // Convolve breit-wigner with models
    RooFFTConvPdf* fun1 = new RooFFTConvPdf("convPdf1","convPdf1", *m, *wbWigner1, *cb1);
    RooFFTConvPdf* fun2 = new RooFFTConvPdf("convPdf2","convPdf2", *m, *wbWigner2, *cb2); 
    RooFFTConvPdf* fun3 = new RooFFTConvPdf("convPdf3","convPdf3", *m, *wbWigner2, *gauss2); 
    // Set the range in which the convolution is performed
    fun1->setBufferStrategy(RooFFTConvPdf::Extend);
    fun1->setBufferFraction(buff); 
    fun2->setBufferStrategy(RooFFTConvPdf::Extend);
    fun2->setBufferFraction(buff);
    fun3->setBufferStrategy(RooFFTConvPdf::Extend);
    fun3->setBufferFraction(buff);

    std::cout << "Comining pdfs." << std::endl;
    // Combine the pdfs
    RooAddPdf* sigmodel;
    if (fittype =="2CB") {
        sigmodel = new RooAddPdf("sigmodel", "sigmodel", RooArgList(*fun1, *fun2), *f);
    }
    else {
        sigmodel = new RooAddPdf("sigmodel", "sigmodel", RooArgList(*fun1, *fun3), *f);
    }

    // Fit model to data
    RooFitResult* fresult = sigmodel->fitTo(*data, Save(), RooFit::Optimize(kFALSE)); 
    std::cout << "Fit complete." << std::endl;
    

    // Make a pretty plot
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
    sigmodel->plotOn(mframe, Components("convPdf1"), LineColor(4), LineStyle(2));
    if (fittype =="2CB") {
      sigmodel->plotOn(mframe, Components("convPdf2"), LineColor(3), LineStyle(3));
    }
    else {
      sigmodel->plotOn(mframe, Components("convPdf3"), LineColor(3), LineStyle(3));
    }
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
    can1->SaveAs("resolution_fit.pdf");

    /*
    ROOT.gPad.SetTopMargin(0)
    ROOT.gPad.SetLeftMargin(0.15)
    ROOT.gPad.SetRightMargin(0.035)
    ROOT.gPad.SetPad(.01,.01,.95,.77)
    frame.SetTitle("Fitted Data Bmass")
    frame.SetMaximum(frame.GetMaximum()*1.1)
    frame.GetYaxis().SetTitleOffset(1.6)
    frame.Draw()
    # Plot pulls
    c.cd(1)
    ROOT.gPad.SetTopMargin(0)
    ROOT.gPad.SetLeftMargin(0.15)
    ROOT.gPad.SetRightMargin(0.035)
    ROOT.gPad.SetPad(.01,.76,.95,.97)
    # Determine pulls and format
    h_pull = frame.pullHist()
    h_pull.SetFillColor(15)
    h_pull.SetFillStyle(3144)
    # Add pulls to frame
    frame_pull.addPlotable(h_pull,'L3')
    frame_pull.GetYaxis().SetNdivisions(505)
    frame_pull.GetYaxis().SetLabelSize(0.20)
    frame_pull.SetTitle("")
    frame_pull.Draw()
    */

    // Print resolution
    if (fittype =="2CB") {
        std::cout << "CB sigma1: " << s1->getValV() << " CB sigma2: " << s2->getValV() << std::endl;
    }
    else {
        std::cout << "Gauss sigma: " << sg->getValV() << std::endl;
        std::cout << "Error: " << sg->getError() << std::endl;
    }
}




