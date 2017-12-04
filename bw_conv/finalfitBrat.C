#include "TROOT.h"
#include "TMath.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooAbsPdf.h"
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
#include "RooGenericPdf.h"
#include "RooProdPdf.h"
#include "RooVoigtian.h"
#include "RooDecay.h"
#include "RooStats/ModelConfig.h"
#include "TCanvas.h"
#include "TTree.h"
#include "RooWorkspace.h"
#include <map>
#include "TCut.h"
#include "TCanvas.h"
#include "RooChebychev.h"
#include "RooDstD0BG.h"
#include "TFile.h"
#include "TF1.h"
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
//#include "RooStudentT.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "RooStats/SPlot.h"
#include "Braatan.h"
#include "X3872_constants.h"
using namespace RooFit;
using namespace RooStats;

const int numberofbins = 80;

typedef struct {

  double value;
  double error;
 

} ValueWithError;

typedef struct {
  
  double minll;
  double chisqdof ;
  double probchi2;

  int status;
  int quality;
  ValueWithError gRe;
  ValueWithError gIm;
  ValueWithError sigma; 
  ValueWithError nsig;
  ValueWithError nback;
  ValueWithError backpar1;
  ValueWithError backpar2;

} Result;

void fillValueWithError(ValueWithError* val,RooRealVar* var){
  val->value = var->getVal();
  val->error = var->getError();
}

void addBlurb(TCanvas* can, std::string blurb){

  std::cout << "Adding Blurb " << can->GetName() << std::endl;
  TLatex *myLatex = new TLatex(0.5,0.5,"");
  myLatex->SetTextFont(132); 
  myLatex->SetTextColor(1); 
  myLatex->SetTextSize(0.06); 
  myLatex->SetNDC(kTRUE); 
  myLatex->SetTextAlign(11);  
  myLatex->SetTextSize(0.06); 
  myLatex->DrawLatex(0.55, 0.75,blurb.c_str()); 
   
}

//const std::string bkgfile11 =  "newmva11eta/allbackground_0.4_Clone_f.root";
//const std::string bkgfile12 =  "newmva12eta/allbackground_0.35_Clone_f.root";
const unsigned int nbin = 6;
std::string l0names[nbin] = {"newmva11eta/allsignal_0.4_Clone_f_angle_tune_cut_low_l0_small.root","newmva12eta/allsignal_0.35_Clone_f_angle_tuned_cut_low_l0_small.root",
			    "newmva11eta/allsignal_0.4_Clone_f_angle_tune_cut_mid_l0_small.root","newmva12eta/allsignal_0.35_Clone_f_angle_tuned_cut_mid_l0_small.root",
                            "newmva11eta/allsignal_0.4_Clone_f_angle_tune_cut_high_l0_small.root","newmva12eta/allsignal_0.35_Clone_f_angle_tuned_cut_high_l0_small.root"};

std::string fnames[nbin] = {"newmva11eta/allsignal_0.4_Clone_f_angle_tune_cut_low_small.root","newmva12eta/allsignal_0.35_Clone_f_angle_tuned_cut_low_small.root",
			    "newmva11eta/allsignal_0.4_Clone_f_angle_tune_cut_mid_small.root","newmva12eta/allsignal_0.35_Clone_f_angle_tuned_cut_mid_small.root",
                            "newmva11eta/allsignal_0.4_Clone_f_angle_tune_cut_high_small.root","newmva12eta/allsignal_0.35_Clone_f_angle_tuned_cut_high_small.root"};
/*
std::string fnames[nbin] = {"newmva11eta/allsignal_0.4_Clone_f_angle_tune_cut_low.root","newmva12eta/allsignal_0.35_Clone_f_angle_tuned_cut_low.root",
			    "newmva11eta/allsignal_0.4_Clone_f_angle_tune_cut_mid.root","newmva12eta/allsignal_0.35_Clone_f_angle_tuned_cut_mid.root",
                            "newmva11eta/allsignal_0.4_Clone_f_angle_tune_cut_high.root","newmva12eta/allsignal_0.35_Clone_f_angle_tuned_cut_high.root"};
*/
const std::string mbranch = "mgev";

/*const double sval1[nbin] = {2.65e-3, 2.68e-3};
const double sval2[nbin] = {3.98e-3, 4.01e-3};
*/
const double sval1[nbin] = {2.43e-3, 2.43e-3, 2.67e-3,2.68e-3, 3.03e-3, 3.05e-3};
const double sval2[nbin] = {3.52e-3, 3.51e-3, 3.91e-3, 3.91e-3, 3.98e-3, 4.01e-3};

const double psi_pdg =  3686.097e-3;

const double fval[nbin] = {0.7,0.7, 0.7, 0.7, 0.7, 0.7};


//checked
const double dcb_cbnval1[nbin]  = {1.9, 1.8, 1.9, 1.8, 1.9, 1.8 };
const double cbg_cbnval1[nbin]  = {1.41, 1.49, 1.23, 1.18, 0.99, 0.98 };
//const double cbg_cbnval1[nbin]  = {1.41, 1.47, 1.22, 1.18, 1.03, 0.99 };
const double dcb_cbenval1[nbin]  = {0.1,0.1,0.1, 0.1,0.1,0.1};
const double cbg_cbenval1[nbin]  = { 0.08,0.08,0.08, 0.08,0.09,0.08};

const double cbnval2[nbin]  = {1, 1, 1, 1,1,1};

//checked
const double dcb_scaleval[nbin] = {1.037, 1.0184, 1.032, 1.013, 1.0456, 1.0238};
const double dcb_escaleval[nbin] = {0.007, 0.007, 0.007, 0.007, 0.007, 0.006};
const double cbg_scaleval[nbin] = {1.027, 1.007, 1.027, 1.008, 1.043, 1.02 };
//const double cbg_scaleval[nbin] = {1.027, 1.007, 1.027, 1.005, 1.042, 1.02 };
const double cbg_escaleval[nbin] = {0.008, 0.007 , 0.007, 0.006, 0.008, 0.006};

// checked
const double psi_dcb_scaleval[nbin] = {1.010,1.016, 1.046, 1.0369, 1.0286 , 1.006 };
const double psi_dcb_escaleval[nbin] = {0.007, 0.006, 0.0082, 0.006, 0.009, 0.007};
const double psi_cbg_scaleval[nbin] = {1.033, 1.032, 1.054,  1.047, 1.037, 1.011};
const double psi_cbg_escaleval[nbin] = {0.008, 0.006, 0.009, 0.006, 0.009, 0.006};

const double psi_mcscale_cbg[nbin]  = {0.989, 1.002, 0.981, 0.981, 1.003,1.007};

const double psi_emcscale_cbg[nbin]  = {0.011, 0.009, 0.011, 0.007 ,0.012, 0.009 };


const double psi_mcscale_dcb[nbin]  = {1.014, 1.02, 0.993, 0.996, 1.014, 1.011 };
const double psi_emcscale_dcb[nbin]  = {0.008, 0.007, 0.01, 0.0080 , 0.012, 0.009 };

const double c_psx[6] = {3.89816e-01 , 4.61556e-01, 1.2, 1.2, 99, 100};

double tot_dcb_scaleval[nbin];
double tot_dcb_escaleval[nbin];
double tot_cbg_scaleval[nbin] ;
double tot_cbg_escaleval[nbin];

std::string blurbs[nbin] = { "#splitline{LHCb Preliminary}{#scale[1]{#sqrt{s} = 7 TeV, L =1 fb^{ #kern[0.3]{#font[122]{-}1}}}}" ,  "#splitline{LHCb Preliminary}{#scale[1]{#sqrt{s} = 8 TeV,  L =2 fb^{ #kern[0.3]{#font[122]{-}1}}}}",
			     "#splitline{LHCb Preliminary}{#scale[1]{#sqrt{s} = 7 TeV, L =1 fb^{ #kern[0.3]{#font[122]{-}1}}}}" ,  "#splitline{LHCb Preliminary}{#scale[1]{#sqrt{s} = 8 TeV, L =2 fb^{ #kern[0.3]{#font[122]{-}1}}}}",
			     "#splitline{LHCb Preliminary}{#scale[1]{#sqrt{s} = 7 TeV, L =1 fb^{ #kern[0.3]{#font[122]{-}1}}}}" ,  "#splitline{LHCb Preliminary}{#scale[1]{#sqrt{s} = 8 TeV, L =2 fb^{ #kern[0.3]{#font[122]{-}1}}}}"
};

// checked
const double escale_extra[nbin] = {0.4e-2,0.4e-2,0.7e-2, 0.5e-2,1e-2,1e-2};

const std::string labels[nbin] = {"11_low","12_low",
				  "11_mid","12_mid",
				  "11_high", "12_high"};


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


void totScale(){

  for (int i = 0; i < nbin; ++i){
    tot_dcb_scaleval[i] = psi_dcb_scaleval[i]* dcb_scaleval[i]    ;
    tot_dcb_escaleval[i] =  tot_dcb_scaleval[i]* sqrt(std::pow( psi_dcb_escaleval[i]/ psi_dcb_scaleval[i],2) + std::pow(dcb_escaleval[i]/dcb_scaleval[i],2) + std::pow(escale_extra[i],2) +std::pow(psi_emcscale_dcb[i]/psi_mcscale_dcb[i], 2.));
    //std::cout << "***" << labels[i] << "***" << std::endl;
    //std::cout << " scale " << labels[i] << " "  <<  tot_dcb_scaleval[i] <<  " " <<  tot_dcb_escaleval[i] << std::endl;
    //std::cout << "*** "<< tot_dcb_scaleval[i] << " " << dcb_escaleval[i] << std::endl;
  }


  for (int i = 0; i < nbin; ++i){
    tot_cbg_scaleval[i] = psi_cbg_scaleval[i]* cbg_scaleval[i]  ;
    tot_cbg_escaleval[i] =  tot_cbg_scaleval[i]* sqrt(std::pow( psi_cbg_escaleval[i]/ psi_cbg_scaleval[i],2) + std::pow(cbg_escaleval[i]/cbg_scaleval[i],2) + std::pow(escale_extra[i],2) + std::pow(psi_emcscale_cbg[i]/psi_mcscale_cbg[i], 2.));
    //std::cout << "***" << labels[i] << "***" << std::endl;
    //std::cout << " scale cbg " << labels[i] << "  "<<  tot_cbg_scaleval[i] <<  " " <<  tot_cbg_escaleval[i] << std::endl;
    //std::cout << tot_cbg_scaleval[i] << " " << cbg_escaleval[i] << std::endl;
  }
  
}

ValueWithError rmsDCB(int ibin){

  ValueWithError val;
  val.value =  tot_dcb_scaleval[ibin]*sqrt(0.7*sval1[ibin]*sval1[ibin] + 0.3*sval2[ibin]*sval2[ibin])  ;
  val.error = val.value* tot_dcb_escaleval[ibin]/tot_dcb_scaleval[ibin];
  return val;
}

ValueWithError rmsCBG(int ibin){

  ValueWithError val;
  val.value =  tot_cbg_scaleval[ibin]*sqrt(0.7*sval1[ibin]*sval1[ibin] + 0.3*sval2[ibin]*sval2[ibin])  ;
  val.error = val.value* tot_cbg_escaleval[ibin]/tot_cbg_scaleval[ibin];
  return val;
}


void printscale(){

  totScale(); 
 for (int i = 0; i < nbin; ++i){
   std::cout << "***" << labels[i] << "***" << std::endl;
    ValueWithError val =  rmsDCB(i);
   std::cout << " scale " << labels[i] << " "  <<  tot_dcb_scaleval[i] <<  " " <<  tot_dcb_escaleval[i] << std::endl;
   std::cout << "*** "<< tot_dcb_scaleval[i] << " " << dcb_escaleval[i] << " rms " << val.value<< std::endl;
  }


  for (int i = 0; i < nbin; ++i){
    ValueWithError val =  rmsCBG(i);
    std::cout << "***" << labels[i] << "***" << std::endl;
    std::cout << " scale cbg " << labels[i] << "  "<<  tot_cbg_scaleval[i] <<  " " <<  tot_cbg_escaleval[i] << std::endl;
    std::cout << tot_cbg_scaleval[i] << " " << cbg_escaleval[i]  << " rms " << val.value<< std::endl;
  }
 
  
}

const double polx2_a[nbin] = {3.2, 4.2, 3.9, 3.8, 4.3, 4.1 };
const double polx2_c1[nbin] = {0.13, 0.10, 0.12, 0.12, 0.12, 0.12};

const double polx3_a[nbin] = {2.8, 3,3.5 ,3.5 ,4,4};


const double polx_c1[nbin] = { 0.13, 0.13, 0.16, 0.16, 0.21, 0.2  };
const double cheb_c1[nbin] = {-0.02,-0.01, 0.05, -0.01, 0.09, 0.08 };
const double expexp_c[nbin] = { -0.5, -0.15, 1.2, 1.2, 2.4,2.1  };

const double c0[nbin] = {-0.02, -0.006, 0.05, 0.05, 0.05, 0.08};
const double c1[nbin] = {-0.02, -0.03, -0.03, 0.001, -0.03, 0.008};

const double t_a[nbin] = {6.98109e+01, 7.27265e+01, 6.07986e-01};
const double t_b[nbin] = {4.10447e-01,5.22739e-01,  3.57502e+00};

double nval[nbin] = {1780, 3425 , 1794, 3698, 1680, 3261};
double bval[nbin] = {19928, 38581, 19426,37180, 16899, 32958};

//double bval[nbin] = {0.9*19928, 0.9*38581, 0.9*19426,0.9*37180, 0.9*16899, 0.9*32958};

const double psimass_cbg[nbin] = {3.68597, 3.68601,3.68598, 3.68602, 3.6861 ,3.68609 };
const double epsimass_cbg[nbin] =  {1.593e-5, 1.2e-5, 1.91e-5, 1.4e-5, 2.5e-5, 1.72e-5};

const double psimass_dcb[nbin] = {3.68597, 3.6860,3.68598,3.68601 , 3.68609, 3.68609 };
const double epsimass_dcb[nbin] =  {1.59124e-5,  1.2e-5, 1.9e-5 , 1.4e-5, 2.5e-5, 1.73e-5};






const double mval = 3.8715;
const double avalue = 0.717;
const double frac = 0.08;
const std::string mstring = "1.63 + 209*";
const std::string par1= "2.28";
const std::string par2= "7.9";
const std::string par3= "0.464*1000";

double nentry[nbin];


void drawSlice(std::string type, std::string btype, const std::string& blurb, const std::string& label, RooCategory& sample, RooRealVar& m, RooDataSet& data, RooSimultaneous& simpdf ,double minmass, double maxmass){

  std::string canname = std::string("slice")  + label;
  TCanvas *can1 = new TCanvas(canname.c_str(), canname.c_str(),10,44,600,400);
  
  RooPlot* frame = m.frame(Bins(50));
  std::string slicename = std::string("sample==sample::") + label;
  data.plotOn(frame,Cut(slicename.c_str()));

  std::string signame = "signal_"+label;
  std::string bgname = "bg_"+label;
  
  //std::string signame2 = "gaussx1" + label;
  simpdf.plotOn(frame,Slice(sample,label.c_str()),ProjWData(sample,data), LineColor(4),  LineStyle(3), Components(bgname.c_str())) ;
  simpdf.plotOn(frame,Slice(sample,label.c_str()),ProjWData(sample,data), LineColor(1),  LineStyle(2), Components(signame.c_str())) ;
  simpdf.plotOn(frame,Slice(sample,label.c_str()),ProjWData(sample,data),  LineStyle(1),  LineColor(2)) ;

  frame->SetTitle("");
  TAxis* xachse = frame->GetXaxis(); TAxis* yachse = frame->GetYaxis();
  xachse->SetTitleFont (132);
  yachse->SetTitleFont (132);
  xachse->SetLabelFont (132);
  yachse->SetLabelFont (132); 
  xachse->SetTitle("m(J/#psi #pi^{+} #pi^{#font[122]{ -} }) [GeV/c^{2}]"); 
  std::string aystring = massCandidates(minmass,maxmass,numberofbins); 
  yachse->SetTitle(aystring.c_str()); 
  yachse->SetTitleOffset(0.85); 
  xachse->SetTitleOffset(0.9);
  xachse->SetTitleSize(0.05); 
  yachse->SetTitleSize(0.055);

  frame->Draw();
  
  addBlurb(can1,blurb);
  std::string epsname = "sim_" + type + "_" + btype +"_" + label + ".eps";

  can1->Print(epsname.c_str());
  std::string pdfname = "epstopdf " + epsname;
  gSystem->Exec(pdfname.c_str());
  
 
}

void filltree(std::string files[6], int ibin, RooCategory& cat, RooRealVar& m, RooDataSet& data, std::string mbranch , double minmass, double maxmass){

 std::string filename = files[ibin];
 TFile* theFile = new TFile(filename.c_str());
 TTree* tree = (TTree*)theFile->Get("psiCand");
 nentry[ibin] = tree->GetEntries();
 // nval[ibin] = tree->GetEntries()*0.1;
 //bval[ibin] = tree->GetEntries()*0.9;
 std::string label = labels[ibin];

 //std::stringstream stream;
 // streams << mbranch << " >" << minmass << "&&" << mbranch << "< " << maxmass;;
 //TCut cut = stream.str().c_str();
 

 unsigned int i= 0;
 double mass; tree->SetBranchAddress(mbranch.c_str(),&mass);
 unsigned int nentries = tree->GetEntries();
 for (int i = 0; i < nentries; ++i){
   tree->GetEntry(i);
   if (mass > minmass && mass < maxmass){
     m.setVal(mass);
     cat.setLabel(label.c_str());
     data.add(RooArgSet(m,cat));
   }
}

 theFile->Close();
 
}

RooAbsPdf* makepdf(int ibin, RooRealVar& m, double minmass, double maxmass, RooRealVar* gRe,
		   RooRealVar* gIm, bool fixRe = false , bool fixIm = false, bool fixScale =false, std::string type = "DCB", std::string btype= "cheb2",
		   RooRealVar* psivar =0, RooRealVar* scale=0, RooRealVar* n1 = 0){

  std::cout << "fix scale " << fixScale << std::endl;

  
    double thresval;
    if (type == "DCB"){
      thresval = psimass_dcb[ibin] - psi_pdg + X3872::m_D0  + X3872::m_D0star ;
    }
    else {
      thresval = psimass_cbg[ibin] - psi_pdg + X3872::m_D0  + X3872::m_D0star;
    }

    //     thresval =  X3872::m_D0  + X3872::m_D0star;
    // std::cout << "Threshold " << thresval << std::endl;
   std::string scalename = scale->GetName();

   
   std::string bwname = "brat_" + labels[ibin];
   RooRealVar* thres = new RooRealVar("thres","thres", thresval);
   std::cout << " *************** " <<std::endl;
   gRe->Print();
   X3872::Braatan* brat = new X3872::Braatan(bwname.c_str(), bwname.c_str() , m, *gRe, *gIm , *thres); 
  
   std::string sigma1name = "sigma1var" + labels[ibin] ;
   std::stringstream tval1; tval1 << sval1[ibin];
   RooFormulaVar* s1 = new RooFormulaVar(sigma1name.c_str(), sigma1name.c_str(), scaledSigma(tval1.str(),scalename).c_str(),*scale);
   
   std::string sigma2name = "sigma2var" + labels[ibin] ;
   std::stringstream tval2; tval2 << sval2[ibin];
   RooFormulaVar* s2 = new RooFormulaVar(sigma2name.c_str(), sigma2name.c_str(), scaledSigma(tval2.str(),scalename).c_str(),*scale);

   std::string fname = "f_" + labels[ibin] ;
   RooRealVar* f = new RooRealVar("f","f",fval[ibin], 0,1); f->setConstant(true);

   std::string atname1 = "atail" + labels[ibin];
   RooFormulaVar* ac1 = tailModelScale(*scale,tval1.str(),atname1.c_str(),par1, par2, par3);

   std::string n2name = "n2_" + labels[ibin];
   RooRealVar* n2 = new RooRealVar("n2","n2",cbnval2[ibin], 0,1); n2->setConstant(true);
   std::string atname2 = "atail2" + labels[ibin];
   RooFormulaVar* ac2 = tailModelScale(*scale,tval2.str(),atname2.c_str(),par1, par2, par3);

   RooRealVar* meanres = new RooRealVar("meanres","meanres", 0);

   std::string cbname1 = "cb1_" + labels[ibin];
   RooCBShape* cb1 = new RooCBShape(cbname1.c_str(),cbname1.c_str(),m, *meanres,*s1,*ac1,*n1);
   std::string cbname2 = "cb2_" + labels[ibin];
   RooCBShape* cb2 = new RooCBShape(cbname2.c_str(),cbname2.c_str(),m, *meanres,*s2,*ac2,*n2);

   std::string sbname = "gauss" + labels[ibin];
   RooGaussian* gauss = new RooGaussian(sbname.c_str(), sbname.c_str(), m, *meanres ,*s2);

   std::string csname1 = "convpdf1_" + labels[ibin] ;
   RooFFTConvPdf* tfun1 = fftpdf(csname1, m,  brat,cb1);
 
   std::string csname2 = "convpdf2_" + labels[ibin] ;
   RooFFTConvPdf* tfun2 = fftpdf(csname2, m,  brat,cb2);

   std::string csname3 = "convpdf3_" + labels[ibin] ;
   RooFFTConvPdf* tfun3 = fftpdf(csname3, m,  brat,gauss);

   // make the signal pdfs
   std::string signame = "signal_"+labels[ibin];  RooAddPdf* sigpdf ;
   if (type == "DCB"){
     sigpdf = new RooAddPdf(signame.c_str(), signame.c_str(),RooArgList(*tfun1,*tfun2),*f);
   }
   else {
     sigpdf = new RooAddPdf(signame.c_str(), signame.c_str(),RooArgList(*tfun1,*tfun3),*f);
   }
   

   // and the background pdfs
   std::string bgname = "bg_"+labels[ibin];
   RooAbsPdf* bg;
   if (btype == "cheb2"){
     bg = cheb2(m,bgname, false, false, c0[ibin],c1[ibin], -0.75, 0.75);
   }
   else if (btype == "cheb1") {
     bg = cheb1(m, bgname, false, cheb_c1[ibin], -0.2, 0.5);
   }
   else if (btype == "exp"){
     std::string bgname = "exp_"+labels[ibin];
     std::string bgname_var = "c0_"+labels[ibin];
     RooRealVar* tauback = new RooRealVar(bgname_var.c_str(),bgname_var.c_str(), expexp_c[ibin], -5, 5.);
     RooExponential* bexp = new RooExponential(bgname.c_str(),bgname.c_str(),m, *tauback);
     bg = bexp;
   }
   else if (btype == "polx"){
     RooRealVar* a = new RooRealVar("a","a",3.6); //a.setConstant(true); 
     //RooRealVar* b = new RooRealVar("b","b", -0.1, -1,0);// b.setConstant(true);
     std::string c0name = "c0_"+labels[ibin];
     RooRealVar* c0 = new RooRealVar(c0name.c_str(),c0name.c_str(), polx_c1[ibin], -0.05, 0.45);
     RooRealVar* mr = new RooRealVar("mr","mr", 3376.e-3, 2000e-3,4000e-3); mr->setConstant(true);
     RooGenericPdf* tbg = new RooGenericPdf(bgname.c_str(), bgname.c_str(), "TMath::Power((@0-@3),@1) *TMath::Exp(-@0/@2)" , RooArgList(m,*a,*c0,*mr));
     bg = tbg; 
  }
   else if  (btype == "polx2"){
     std::string aname = "a_"+labels[ibin];
     RooRealVar* a = new RooRealVar(aname.c_str(),aname.c_str(),polx2_a[ibin], 1, 6); //a.setConstant(true); 
     //RooRealVar* b = new RooRealVar("b","b", -0.1, -1,0);// b.setConstant(true);
     std::string c0name = "c0_"+labels[ibin];
     RooRealVar* c0 = new RooRealVar(c0name.c_str(),c0name.c_str(), 0.16, -0.05, 0.45); c0->setConstant(true);
     RooRealVar* mr = new RooRealVar("mr","mr", 3376.e-3, 2000e-3,4000e-3); mr->setConstant(true);
     RooGenericPdf* tbg = new RooGenericPdf(bgname.c_str(), bgname.c_str(), "TMath::Power((@0-@3),@1) *TMath::Exp(-@0/@2)" , RooArgList(m,*a,*c0,*mr));
     bg = tbg; 
  }
   else if (btype == "polx3"){
     std::string aname = "a_"+labels[ibin];
     RooRealVar* a = new RooRealVar(aname.c_str(),aname.c_str(),polx3_a[ibin], 1, 7); //a.setConstant(true);
     //RooRealVar* b = new RooRealVar("b","b", -0.1, -1,0);// b.setConstant(true);
     std::string c0name = "c0_"+labels[ibin];
     RooRealVar* c0 = new RooRealVar(c0name.c_str(),c0name.c_str(), 0.17, 0., 0.45); c0->setConstant(true);
     RooRealVar* mr = new RooRealVar("mr","mr", 3376.e-3, 2000e-3,4000e-3); mr->setConstant(true);
       RooRealVar* mr2 = new RooRealVar("mr2","m2", 5.3, 2000e-3,6000e-3); mr2->setConstant(true);
     
       RooGenericPdf* tbg = new RooGenericPdf(bgname.c_str(), bgname.c_str(), "TMath::Exp(-@0/@2)*TMath::Power((@0*@0-@3*@3)*(@4*@4 - @0*@0), @1)" , RooArgList(m,*a,*c0,*mr,*mr2));
     bg = tbg; 
   }
   
   else if (btype == "pol"){
     //RooRealVar* a = new RooRealVar("a","a",3.6, 1,6); //a.setConstant(true); 
     //RooRealVar* b = new RooRealVar("b","b", -0.1, -1,0);// b.setConstant(true);
     std::string c0name = "c0_"+labels[ibin];
     RooRealVar* c0 = new RooRealVar(c0name.c_str(),c0name.c_str(),1, -0.9, 2);
     RooRealVar* mr = new RooRealVar("mr","mr", 3376.e-3, 2000e-3,4000e-3); mr->setConstant(true);
     RooGenericPdf* tbg = new RooGenericPdf(bgname.c_str(), bgname.c_str(), "TMath::Power((@0-@2),@1)" , RooArgList(m,*c0,*mr));
     bg = tbg;
   }

   else if (btype == "psx"){
      std::string c0name = "c0_"+labels[ibin];
     RooRealVar* c0 = new RooRealVar(c0name.c_str(),c0name.c_str(), c_psx[ibin], 0.001, 300);
     RooRealVar* mr = new RooRealVar("mr","mr", 3376.e-3, 2000e-3,4000e-3); mr->setConstant(true);
     RooGenericPdf* tbg = new RooGenericPdf(bgname.c_str(), bgname.c_str(), "TMath::Sqrt(TMath::Power(@0-@2,2))*TMath::Exp(-@0/@1)" , RooArgList(m,*c0,*mr));
     bg = tbg;
   }
   else if (btype == "pol2"){
     std::string c0name = "c0_"+labels[ibin];
     RooRealVar* c0 = new RooRealVar(c0name.c_str(),c0name.c_str(), -0.84, -2 , 1);
     RooRealVar* mr = new RooRealVar("mr","mr", 3376.e-3, 2000e-3,4000e-3); mr->setConstant(true);
     RooGenericPdf* tbg = new RooGenericPdf(bgname.c_str(), bgname.c_str(), "@1*TMath::Power(@0-@2,2) + @0-@2" , RooArgList(m,*c0,*mr));
     bg = tbg;
   }
   else if (btype == "cheb2f"){ 
     RooChebychev* bgchebf = cheb2(m,bgname, false, true, 0.3, -0.11, -2,2, -2,2);
     bg = bgchebf;
   }
    
   // full pdf
   std::string pdfname = "model" + labels[ibin];
   std::string nsigname = "nsig_" + labels[ibin];
   std::string nbgname = "nbg_" + labels[ibin];
   RooRealVar* nsig = new RooRealVar(nsigname.c_str(),nsigname.c_str(), nval[ibin], 0 , nentry[ibin]);
   RooRealVar* nback = new RooRealVar(nbgname.c_str(),nbgname.c_str(), bval[ibin],  0 ,  nentry[ibin]);
   RooAddPdf* model = new RooAddPdf(pdfname.c_str(),pdfname.c_str(), RooArgList(*sigpdf,*bg), RooArgList(*nsig,*nback));

   //   std::string fpdfname = "fmodel" + labels[ibin];
   //RooAbsPdf* fmodel = new RooProdPdf(fpdfname.c_str(),"model", RooArgSet(*model, *gscale, *gpsi));
   //model->Print();
   //nsig->Print();
   //nback->Print();
   //   return fixScale ==false ? fmodel : model ;
   return model;
}



Result* datafitbwcon(std::string files[6], bool makeplots = true,std::string type = "DCB", std::string btype = "polx", double re = 10e-3, double im = 12e-3, bool fixRe = false , bool fixIm = false,  bool fixScale = true, double mval = 0.186, bool closefiles = false,  double minmass = 3.832, double maxmass = 3.912){

  totScale();

 Result* res = new Result();

 RooCategory sample("sample","sample") ;
 int ibin;
 for (ibin =0 ; ibin < nbin; ++ibin){
   sample.defineType(labels[ibin].c_str()) ;
 }
 RooRealVar m("m","m",minmass, maxmass);
 m.setBins(50000,"cache"); // crucial line !!!
 
 // dataset
 RooDataSet data("DS", "DS", RooArgSet(m,sample));

 for (ibin = 0; ibin < nbin; ++ibin){
   filltree(files,ibin,sample, m, data, mbranch, minmass, maxmass );
 }

 RooArgSet cons = RooArgSet("cons");
 
 // breit wigner stuff
 RooRealVar* gRe =new RooRealVar("gRe", "gRe",re, -0.1, 1); gRe->setConstant(fixRe);
 RooRealVar* gIm =new RooRealVar("gIm", "gIm", im, -0.1,1); gIm->setConstant(fixIm);
 
 
 RooSimultaneous simPdf("simPdf","simultaneous pdf",sample);

  // add them all together
  std::vector<RooAbsPdf*> pdfs;
  for (ibin =0 ; ibin < nbin; ++ibin ){


      // switch for models
    double tscalevar;  double etscalevar; double tnvar; double etnvar; double  psimassvar; double epsimassvar;
    if (type == "DCB"){
      tscalevar = tot_dcb_scaleval[ibin];
      etscalevar =  tot_dcb_escaleval[ibin];
      tnvar = dcb_cbnval1[ibin];
      etnvar = dcb_cbenval1[ibin];
      psimassvar = psimass_dcb[ibin];
      epsimassvar = epsimass_dcb[ibin];
    }
    else {
      tscalevar = tot_cbg_scaleval[ibin];
      etscalevar = tot_cbg_escaleval[ibin];
      tnvar = cbg_cbnval1[ibin];
      etnvar = cbg_cbenval1[ibin];
      psimassvar = psimass_cbg[ibin];
      epsimassvar = epsimass_cbg[ibin];
    }
    
     // constraint  on psi  
    std::string psivarname = "psi"+ labels[ibin];
    RooRealVar* psivar = new RooRealVar(psivarname.c_str(), psivarname.c_str(), psimassvar, psimassvar  - epsimassvar , psimassvar + 5*epsimassvar);
    psivar->setConstant(true);
    std::string psiconname = "conpsi" + labels[ibin]; 
    RooGaussian* gpsi = new RooGaussian(psiconname.c_str(),psiconname.c_str(),*psivar,RooConst(psimassvar),RooConst(epsimassvar));
   // psiptr = psivar; 
    //cons.add(*gpsi);     

    std::string scalename = "scalevar" + labels[ibin];
   std::string scaleconname = "scalecon" + labels[ibin];

   std::string n1name = "n1_" + labels[ibin];
   std::string n1namecon = "n1con_" + labels[ibin];
   RooRealVar* n1 = new RooRealVar(n1name.c_str(),n1name.c_str(),tnvar, 0.1,20); n1->setConstant(fixScale);
   RooGaussian* nscale = new RooGaussian( n1namecon .c_str(), n1namecon .c_str(),*n1,RooConst(tnvar),RooConst(etnvar));
   cons.add(*nscale);
   
   RooRealVar* scale = new RooRealVar(scalename.c_str(), scalename.c_str(), tscalevar, 0.75, 1.75);
   //scaleptr = scale ;
   RooGaussian* gscale = new RooGaussian(scaleconname.c_str(),scaleconname.c_str(),*scale,RooConst(tscalevar),RooConst(etscalevar));
   scale->setConstant(fixScale); 
   cons.add(*gscale);
    
    RooRealVar* scaleptr;
    RooAbsPdf* model = makepdf(ibin,m,minmass,maxmass, gRe, gIm, fixRe, fixIm, fixScale,type,btype,psivar,scale,n1);
    pdfs.push_back(model);
    simPdf.addPdf(*model,labels[ibin].c_str());
   
  }

  RooFitResult* fresult;
  if (cons.getSize() != 0) {
    fresult = simPdf.fitTo(data, Extended(), Optimize(false), ExternalConstraints(cons), Save());
  }
  else {
    fresult = simPdf.fitTo(data, Extended(), Optimize(false), Save());
  }

  res->minll = fresult->minNll();

   // save the results etc to workspace
   /* RooWorkspace *w = new RooWorkspace("wspace","workspace") ;
   w->import(*fresult);
   w->import(*model);
   std::string output = "result_" + labels[ibin] +".root";
   w->writeToFile(output.c_str());
   */
   
   const RooArgList& fitPars = fresult->floatParsFinal();
   RooRealVar* fitRe = (RooRealVar*)fitPars.find("gRe");
   RooRealVar* fitIm = (RooRealVar*)fitPars.find("gIm");
   std::string nsigname = "nsig_" + labels[ibin];
   std::string nbgname = "nbg_" + labels[ibin];
   RooRealVar* nsig = (RooRealVar*)fitPars.find(nsigname.c_str());
   RooRealVar* nback = (RooRealVar*)fitPars.find(nbgname.c_str());
   if (nsig) fillValueWithError(&res->nsig,nsig);
   if (nback)  fillValueWithError(&res->nback,nback);
   if (fitRe) fillValueWithError(&res->gRe,fitRe);
   if (fitIm)  fillValueWithError(&res->gIm,fitIm);

   std::cout << "MyResult " << type << " " << btype << res->gRe.value << "+/- " << res->gRe.error << " "
             << res->gIm.value << "+/- " << res->gIm.error << std::endl;
   
  if (makeplots == true){
    for (int ibin = 0; ibin < nbin; ++ibin){
      std::string sname = std::string("slice") + labels[ibin];
      drawSlice(type,btype,blurbs[ibin],labels[ibin], sample, m, data, simPdf ,minmass, maxmass);

      //std::string signame2 = "gaussx1" + label;
      /*
      simpdf.plotOn(frame,Slice(sample,label.c_str()),ProjWData(sample,data), LineColor(4),  LineStyle(3), Components(bgname.c_str())) ;
      simpdf.plotOn(frame,Slice(sample,label.c_str()),ProjWData(sample,data), LineColor(1),  LineStyle(1), Components(signame.c_str())) ;
      simpdf.plotOn(frame,Slice(sample,label.c_str()),ProjWData(sample,data),  LineStyle(1),  LineColor(2)) ;

      frame->SetTitle("");
      TAxis* xachse = frame->GetXaxis(); TAxis* yachse = frame->GetYaxis();
      xachse->SetTitleFont (132);
      yachse->SetTitleFont (132);
      xachse->SetLabelFont (132);
      yachse->SetLabelFont (132); 
      xachse->SetTitle("m(J/#psi #pi^{+} #pi^{#font[122]{ -} }) [GeV/c^{2}]"); 
      std::string aystring = massCandidates(minmass,maxmass,numberofbins); 
      yachse->SetTitle(aystring.c_str()); 
      yachse->SetTitleOffset(0.85); 
      xachse->SetTitleOffset(0.9);
      xachse->SetTitleSize(0.05); 
      yachse->SetTitleSize(0.055);

      frame->Draw();
  
      addBlurb(can1,blurb);
      std::string epsname = "sim_" + type + "_btype_" + label + ".eps";

      can->Print(epsname.c_str());
      std::string pdfname = "epstopdf " + epsname;
      gSystem->Exec(pdfname.c_str());
      */
      
    }
  }

  
  
  return res;

}




Result* datafitbwconsingle(std::string filename, std::string epsname, bool makeplots = true,int ibin = 0, std::string type = "DCB", std::string btype = "polx" ,double re = 38.4e-3, double im = 12e-3, bool fixRe = false, bool fixIm = false,  bool fixScale = true, double mval = 0.186, bool closefiles = false,  double minmass = 3.832, double maxmass = 3.912 , bool makeweights = false)
{

  totScale();

 
 TFile* theFile = new TFile(filename.c_str());
 TTree* tree = (TTree*)theFile->Get("psiCand");
 
 Result* res = new Result();

 RooRealVar m("mgev","mgev",minmass, maxmass);
 m.setBins(100000,"cache"); // crucial line !!!

 
 
 // dataset
 std::stringstream stream;
 stream << mbranch << " >" << minmass << "&&" << mbranch << "< " << maxmass;;
 TCut thecut = stream.str().c_str();
 
 RooDataSet data("DS", "DS", RooArgSet(m), Import(*tree), Cut(thecut));

 nentry[ibin] = data.numEntries();
 // nval[ibin] = data.numEntries()*0.1;
 // bval[ibin] = data.numEntries()*0.9;
   

 // breit wigner stuff
 // RooRealVar* gRe =new RooRealVar("gRe", "gRe", re ,  0,1 );
 RooRealVar* gRe =new RooRealVar("gRe", "gRe",re, 0, 1); gRe->setConstant(fixRe);
 RooRealVar* gIm =new RooRealVar("gIm", "gIm", im, 0,1); gIm->setConstant(fixIm);
 // RooRealVar* gIm =new RooRealVar("gIm", "gIm", 0);
 
 
  

  RooArgSet cons("cons");


 // switch for models
  std::stringstream tdm;
  double tscalevar;  double etscalevar;  double tnvar; double etnvar;  double psimassvar;  double epsimassvar;
   if (type == "DCB"){
     tscalevar = tot_dcb_scaleval[ibin];
     etscalevar =  tot_dcb_escaleval[ibin];
     tnvar = dcb_cbnval1[ibin];
     etnvar = dcb_cbenval1[ibin];
     psimassvar = psimass_dcb[ibin];
     epsimassvar = epsimass_dcb[ibin];
     tdm << psimass_dcb[ibin] << "+DM";
   }
   else {
     tscalevar = tot_cbg_scaleval[ibin];
     etscalevar = tot_cbg_escaleval[ibin];
     tnvar = cbg_cbnval1[ibin];
     etnvar = cbg_cbenval1[ibin];
      psimassvar = psimass_cbg[ibin];
      epsimassvar = epsimass_cbg[ibin];
      tdm << psimass_cbg[ibin] << "+DM";
   }
   


  
  std::string psivarname = "psi"+ labels[ibin];
  RooRealVar* psivar = new RooRealVar(psivarname.c_str(), psivarname.c_str(), psimassvar, psimassvar -5*epsimassvar , psimassvar + 5*epsimassvar);
  psivar->setConstant(true);
  std::string psiconname = "conpsi" + labels[ibin]; 
  RooGaussian* gpsi = new RooGaussian(psiconname.c_str(),psiconname.c_str(),*psivar,RooConst(psimassvar),RooConst(epsimassvar));

   std::string scalename = "scalevar" + labels[ibin];
   std::string scaleconname = "scalecon" + labels[ibin];

   RooRealVar* scale = new RooRealVar(scalename.c_str(), scalename.c_str(), tscalevar, 0.75, 1.75);
   //scaleptr = scale ;
   RooGaussian* gscale = new RooGaussian(scaleconname.c_str(),scaleconname.c_str(),*scale,RooConst(tscalevar),RooConst(etscalevar));
   scale->setConstant(fixScale); 
   std::cout << "Const " << tscalevar << " " << etscalevar << " scale " << fixScale <<std::endl; 


   std::string n1name = "n1_" + labels[ibin];
   std::string n1namecon = "n1con_" + labels[ibin];
   RooRealVar* n1 = new RooRealVar(n1name.c_str(),n1name.c_str(),tnvar, 0.1,20); n1->setConstant(fixScale);
   RooGaussian* nscale = new RooGaussian( n1namecon .c_str(), n1namecon .c_str(),*n1,RooConst(tnvar),RooConst(etnvar));
   
   RooAbsPdf* model = makepdf(ibin,m,minmass,maxmass, gRe, gIm, fixRe,fixIm, fixScale,type, btype,psivar,scale,n1);
  
   RooFitResult* fresult;
   if (fixScale == false) {
     fresult = model->fitTo(data, Extended(), Optimize(true), ExternalConstraints(RooArgSet(*gscale,*nscale)), Save());
     //     fresult = model->fitTo(data, Extended(), Optimize(true), ExternalConstraints(RooArgSet(**nscale)), Save(), NumCPU(2));
   }
   else {
     fresult = model->fitTo(data, Extended(), Optimize(true), Save());
   }
   res->minll = fresult->minNll();
   res->quality = fresult->covQual();
   res->status = fresult->status();
   // save the results etc to workspace
   /* RooWorkspace *w = new RooWorkspace("wspace","workspace") ;
   w->import(*fresult);
   w->import(*model);
   std::string output = "result_" + labels[ibin] +".root";
   w->writeToFile(output.c_str());
   */
   
  
   const RooArgList& fitPars = fresult->floatParsFinal();
   RooRealVar* fitRe = (RooRealVar*)fitPars.find("gRe");
   RooRealVar* fitIm = (RooRealVar*)fitPars.find("gIm");
   std::string nsigname = "nsig_" + labels[ibin];
   std::string nbgname = "nbg_" + labels[ibin];
   RooRealVar* nsig = (RooRealVar*)fitPars.find(nsigname.c_str());
   RooRealVar* nback = (RooRealVar*)fitPars.find(nbgname.c_str());
   if (nsig) fillValueWithError(&res->nsig,nsig);
   if (nback)  fillValueWithError(&res->nback,nback);
   if (fitRe) fillValueWithError(&res->gRe,fitRe);
   if (fitIm)  fillValueWithError(&res->gIm,fitIm);    
   
   if (makeplots == true){

    RooPlot* frame = m.frame(Bins(numberofbins));
    
    TCanvas* can = new TCanvas("can","can",600, 400);
    TPad* pad3 = new TPad("pad1","This is pad1",0.05,0.25,0.95,0.97);
    TPad* pad4 = new TPad("pad2","This is pad2",0.05,0.0,0.95,0.25); // pad for the pull
    pad3->Draw();
    pad4->Draw();
    pad3->cd();
     
    
    data.plotOn(frame);
  
  
    std::string signame = "signal_"+labels[ibin];
    model->plotOn(frame,  LineStyle(2),  LineColor(1), Components(signame.c_str())) ;

    std::string bgname = "bg_"+labels[ibin];
    model->plotOn(frame,  LineStyle(3),  LineColor(4), Components(bgname.c_str())) ;
  //model->plotOn(frame,  LineStyle(2),  LineColor(3), Components("convpdf1_11"));
  //model->plotOn(frame,  LineStyle(2),  LineColor(3), Components("convpdf2_11"));

      model->plotOn(frame,  LineStyle(1),  LineColor(2)) ;

    double chisqdof = frame->chiSquare();
    double ndof = numberofbins - fitPars.getSize();
    double probchi2 = TMath::Prob(chisqdof*ndof,ndof);
    res->probchi2 =  probchi2;
    res->chisqdof = chisqdof;
     
      
    frame->SetTitle("");
    TAxis* xachse = frame->GetXaxis(); TAxis* yachse = frame->GetYaxis();
    xachse->SetTitleFont (132);
    yachse->SetTitleFont (132);
    xachse->SetLabelFont (132);
    yachse->SetLabelFont (132); 
    xachse->SetTitle("m(J/#psi #pi^{+} #pi^{#font[122]{ -} }) [GeV/c^{2}]"); 
    std::string aystring = massCandidates(minmass,maxmass,numberofbins); 
    yachse->SetTitle(aystring.c_str()); 
    yachse->SetTitleOffset(0.85); 
    xachse->SetTitleOffset(0.9);
    xachse->SetTitleSize(0.05); 
    yachse->SetTitleSize(0.055);

  
   frame->Draw();
   addBlurb(can,blurbs[ibin]);

    pad4->cd();
    RooPlot* frame2 = m.frame(Bins(numberofbins));
    RooHist* hpull2 = frame->pullHist();
    frame2->SetMinimum(-3.5);
    frame2->SetMaximum(3.5);
    TAxis* xachse4 = frame2->GetXaxis(); TAxis* yachse4 = frame2->GetYaxis();
    xachse4->SetTitleFont (132);
    yachse4->SetTitleFont (132);
    xachse4->SetLabelFont (132);
    yachse4->SetLabelFont (132);
    yachse4->SetTitleOffset(0.15);
    yachse4->SetTitleSize(0.23); //was 0.3
    yachse4->SetLabelSize(0.165); //was 0.13
    xachse4->SetLabelOffset(0.2);
    yachse4->SetNdivisions(105);
    frame2->SetXTitle("");
    frame2->SetYTitle("Pull");
    frame2->SetTitle("");

    frame2->addPlotable(hpull2,"BX") ;
    frame2->Draw() ;

    TLine* zline2 = new TLine(minmass,0, maxmass, 0);
     zline2->SetLineWidth(2);
    zline2->Draw();

    
    
    can->Print(epsname.c_str());
    std::string pdfname = "epstopdf " + epsname;
    gSystem->Exec(pdfname.c_str());

     std::cout << "Fit chi^2" <<  chisqdof << " prob " << probchi2 <<std::endl; 
   }

   if (makeweights == true){
 
     std::string trailer = "_S.root";
     std::string outputName = createOutputName(filename,trailer);
     TFile* f_out  =new TFile(outputName.c_str(),"RECREATE");
     TTree* smalltree = tree->CopyTree(thecut);

     TTree*  newtree = smalltree->CloneTree(-1);
     std::cout << "copied the tree" << std::endl;

     /*
     std::string nsignamesw = "nsig_" + labels[ibin] + "_sw";
     std::string nbgnamesw = "nbg_" + labels[ibin]+ "_sw"; 
     std::string nsignamef = nsignamesw  + "/F";
     std::string nbgnamef =  nbgnamesw  + "/F"; 
     */
    

     std::string nsignameswn = "nsig_" + labels[ibin] + "_sw";
     std::string nbgnameswn = "nbg_" + labels[ibin]+ "_sw";
     
     std::string nsignamesw = "nsig_sw";
     std::string nbgnamesw = "nbg_sw"; 
     std::string nsignamef = nsignamesw  + "/F";
     std::string nbgnamef =  nbgnamesw  + "/F"; 
     
     RooStats::SPlot* sData = new RooStats::SPlot("sData","An SPlot", data, model, RooArgList(*nsig, *nback) );

     float Nsig_sw; float Nbkg_sw;  //float  L_Nsig; float L_Nbkg;
     TBranch*  b_Nsig_sw = newtree->Branch(nsignamesw.c_str() , &Nsig_sw, nsignamef.c_str() );
     TBranch*  b_Nbkg_sw  = newtree->Branch(nbgnamesw.c_str(), &Nbkg_sw,nbgnamef.c_str()) ;

     for (int i = 0; i < data.numEntries(); ++i) {
        newtree->GetEntry(i);
        const RooArgSet* row = data.get(i);
        Nsig_sw =  row->getRealValue(nsignameswn.c_str());
        Nbkg_sw =  row->getRealValue(nbgnameswn.c_str());
        b_Nsig_sw->Fill();
        b_Nbkg_sw->Fill();
     }

     newtree->Write();
     f_out->Close();
   
   }


   if (closefiles == true){
     theFile->Close();
   }
   
 return res;

}

