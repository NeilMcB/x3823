#include <string>
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooAddModel.h"
#include "RooAbsPdf.h"
#include "RooVoigtian.h"
#include "RooAddPdf.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include "RooFormulaVar.h"
#include "RooBreitWigner.h"
#include "RooRelBreitWigner.h"
#include "RooNumConvPdf.h"
#include "RooChebychev.h"
#include "RooDstD0BG.h"
#include "RooBackPdf.h"
#include "RooFFTConvPdf.h"
//#include "RooFlatte.h"
#include "X3872.h"

std::string massCandidates(double minmass, double maxmass, double nbins){

  double bsize = 1000*(maxmass - minmass)/nbins;
  std::stringstream ax; ax << "Candidates/ (" << bsize << "MeV/c^{2})" ;
  return ax.str();
}

std::string varname(std::string var, std::string header){
  return(header + "_"+ var);
}

double resolution(RooRealVar& s,  RooRealVar& s2,  RooRealVar& f) {
  return(sqrt( s.getVal()*s.getVal()*f.getVal() + (1-f.getVal())* s2.getVal()*s2.getVal()  ));
} 

RooAddPdf* cbGaussMC(RooRealVar&m ,double minmass, double maxmass, std::string header ="Model" ){

 std::cout << "Model cbGaussMC" << std::endl;
 RooRealVar* mx = new RooRealVar(varname("mx",header).c_str(),varname("mx",header).c_str(), 0.5*(minmass+maxmass), minmass,maxmass);
 RooRealVar* s2 = new RooRealVar(varname("sigma2",header).c_str(),varname("sigma2",header).c_str(),7,3, 10); // s.setConstant(true);
 RooGaussian* gauss = new RooGaussian(varname("gauss",header).c_str(),varname("gauss(x,mean,sigma)",header).c_str(),m, *mx,*s2);

 RooRealVar* a = new RooRealVar(varname("a",header).c_str(), varname("a",header).c_str(), 2.9, 0., 4); 
 RooRealVar* n = new RooRealVar(varname("n",header).c_str() ,varname("n",header).c_str(), 1., 0., 10); //n->setConstant(true);
 RooRealVar* s = new RooRealVar(varname("sigma",header).c_str(),varname("sigma",header).c_str(),2,0,5 ); // s.setConstant(true);
 RooCBShape* cb = new RooCBShape(varname("cb",header).c_str(),varname("cb",header).c_str(),m, *mx,*s,*a,*n);

 RooRealVar* f = new RooRealVar(varname("f",header).c_str(),varname("f",header).c_str(),0.75, 0.,1.); //f.setConstant(true);
 RooAddPdf* sigmodel= new RooAddPdf(header.c_str(),header.c_str(),RooArgList(*cb,*gauss), *f); 
 return sigmodel;
}


RooAddPdf* cbGauss(RooRealVar&m ,double minmass, double maxmass, std::string header ="Model", 
                   std::string res="X" , double scaleval = 1,  bool fixScale = true, bool fixFrac = true ){

  header += "_" + res;

  double aval,nval,sval,sval2,fval, mval;
  if (res=="X"){
    aval = 2.13;
    nval = 1.53;
    sval = 3.01;
    sval2 = 7.01;
    fval = 0.92;
    mval = 3872.2;
  }
  else {
    aval = 2.49;
    nval = 1.08;
    sval = 2.25;
    sval2 = 5.24;
    fval = 0.89;
    mval = 3686.11;
  }

 std::stringstream form1 ; form1 << varname("scale",header) <<"*" << sval; 
 std::stringstream form2 ; form2 << varname("scale",header) <<"*" << sval2;
 RooRealVar* scale = new RooRealVar(varname("scale",header).c_str(), varname("scale",header).c_str(), scaleval,0,2);  
 if (fixScale == true) scale->setConstant(true);

 RooRealVar* mx = new RooRealVar(varname("mx",header).c_str(),varname("mx",header).c_str(), mval, minmass,maxmass);
 RooRealVar* a = new RooRealVar(varname("a",header).c_str(), varname("a",header).c_str(), aval); a->setConstant(true);
 RooRealVar* n = new RooRealVar(varname("n",header).c_str() ,varname("n",header).c_str(), nval); n->setConstant(true);
 RooFormulaVar* s = new RooFormulaVar(varname("sigma",header).c_str(), varname("sigma",header).c_str(), form1.str().c_str(),*scale);
 RooCBShape* cb = new RooCBShape(varname("cb",header).c_str(),varname("cb",header).c_str(),m, *mx,*s,*a,*n);

 RooFormulaVar* s2 = new RooFormulaVar(varname("sigma2",header).c_str(), varname("sigma2",header).c_str() ,form2.str().c_str(),*scale);
 RooGaussian* gauss = new RooGaussian(varname("gauss",header).c_str(),varname("gauss(x,mean,sigma)",header).c_str(),m, *mx,*s2);

 RooRealVar* f = new RooRealVar(varname("f",header).c_str(),varname("f",header).c_str(),fval, 0.,1.);
 if (fixFrac == true) f->setConstant(true);
 RooAddPdf* sigmodel= new RooAddPdf(header.c_str(),header.c_str(),RooArgList(*cb,*gauss), *f); 
 return sigmodel;
}

RooAddPdf* cbVoight(RooRealVar&m ,double minmass, double maxmass, std::string header ="Model" ,
                    std::string res="X" , double scaleval = 1, bool fixScale = true, bool fixFrac = true){

 header += "_" + res;

 double aval,nval,sval,sval2,fval,wval,mval;
  if (res=="X"){
    aval = 1.7;
    nval = 1.73;
    sval = 2.67;
    sval2 = 4.97;
    fval = 0.69;
    wval = 0.317;
    mval = 3872.2;
  }
  else {
    aval = 1.85;
    nval = 1.47;
    sval = 2.08;
    sval2 = 4.12;
    fval = 0.76;
    wval = 0.304;
    mval = 3686.11;
  }

 std::stringstream form1 ; form1 << varname("scale",header) <<"*" << sval; 
 std::stringstream form2 ; form2 << varname("scale",header) <<"*" << sval2;
 RooRealVar* scale = new RooRealVar(varname("scale",header).c_str(), varname("scale",header).c_str(), scaleval, 0,2);  
 if (fixScale == true) scale->setConstant(true); 

 RooRealVar* mx = new RooRealVar(varname("mx",header).c_str(),varname("mx",header).c_str(), mval, minmass,maxmass);
 RooFormulaVar* s = new RooFormulaVar(varname("sigma",header).c_str(), varname("sigma",header).c_str(), form1.str().c_str(),*scale);
 RooRealVar* width =new RooRealVar(varname("width",header).c_str(), varname("width",header).c_str(), wval, -100., 100); width->setConstant(true);
 RooVoigtian* voight= new  RooVoigtian(varname("voight",header).c_str(),varname("voight(x,mean,sigma)",header).c_str(),m, *mx , *width,*s); 


 RooRealVar* a = new RooRealVar(varname("a",header).c_str(), varname("a",header).c_str(), 2.9, 0., 4);  a->setConstant(true);
 RooRealVar* n = new RooRealVar(varname("n",header).c_str() ,varname("n",header).c_str(), 1., 0., 10); n->setConstant(true);
 RooFormulaVar* s2 = new RooFormulaVar(varname("sigma2",header).c_str(), varname("sigma2",header).c_str(), form2.str().c_str(),*scale);
 RooCBShape* cb = new RooCBShape(varname("cb",header).c_str(),varname("cb",header).c_str(),m, *mx,*s2,*a,*n);

 RooRealVar* f = new RooRealVar(varname("f",header).c_str(),varname("f",header).c_str(),fval, 0.,1.); if (fixFrac) f->setConstant(true); //f.setConstant(true);
 RooAddPdf* sigmodel= new RooAddPdf(header.c_str(),header.c_str(),RooArgList(*voight,*cb), *f); 

 return sigmodel;
}


RooAddPdf* cbVoightMC(RooRealVar&m ,double minmass, double maxmass, std::string header ="Model" , std::string res = "X"){

  double wval; double mval;
  if (res=="X"){
    wval = 0.317;
    mval = 3871.5;
  }
  else {
    wval = 0.304;
    mval = 3686.11;
  }


 RooRealVar* mx = new RooRealVar(varname("mx",header).c_str(),varname("mx",header).c_str(), mval, mval - 5,mval+5);
 RooRealVar* s = new RooRealVar(varname("sigma",header).c_str(),varname("sigma",header).c_str(),2,0, 5); // s.setConstant(true);
 RooRealVar* width =new RooRealVar(varname("width",header).c_str(), varname("width",header).c_str(), 0.0, -2., 1.); width->setConstant(true);
 RooVoigtian* voight= new  RooVoigtian(varname("voight",header).c_str(),varname("voight(x,mean,sigma)",header).c_str(),m, *mx , *width,*s); 


 RooRealVar* a = new RooRealVar(varname("a",header).c_str(), varname("a",header).c_str(), 2.9, 0., 4); 
 RooRealVar* n = new RooRealVar(varname("n",header).c_str() ,varname("n",header).c_str(), 1., 0., 10); //n->setConstant(true);
 RooRealVar* s2 = new RooRealVar(varname("sigma2",header).c_str(),varname("sigma2",header).c_str(),6,3,10 ); // s.setConstant(true);
 RooCBShape* cb = new RooCBShape(varname("cb",header).c_str(),varname("cb",header).c_str(),m, *mx,*s2,*a,*n);

 RooRealVar* f = new RooRealVar(varname("f",header).c_str(),varname("f",header).c_str(),0.75, 0.,1.); //f.setConstant(true);
 RooAddPdf* sigmodel= new RooAddPdf(header.c_str(),header.c_str(),RooArgList(*voight,*cb), *f); 
 return sigmodel;
}


RooBreitWigner* chi0_2P(RooRealVar&m ,double minmass, double maxmass){

  //3918.4 +/- 1.9, 20+/-5 
  RooRealVar* m_chi0_2P = new RooRealVar("m_chi0_2P", "m_chi0_2P", 3918, minmass,maxmass); m_chi0_2P->setConstant(true);
  RooRealVar* s_chi0_2P = new RooRealVar("s_chi0_2P", "s_chi0_2P",20, minmass,maxmass); s_chi0_2P->setConstant(true);
  return new RooBreitWigner("chi0_2P", "chi0_2P", m, *m_chi0_2P, *s_chi0_2P);

}

RooBreitWigner* psi_3770(RooRealVar&m ,double minmass, double maxmass){

  //3918.4 +/- 1.9, 20+/-5 
  RooRealVar* m_chi0_2P = new RooRealVar("m_psi_3770", "m_psi_3770", 3773.15, minmass,maxmass); m_chi0_2P ->setConstant(true);
  RooRealVar* s_chi0_2P = new RooRealVar("s_psi_3770", "m_psi_3770",27.2, minmass,maxmass); s_chi0_2P->setConstant(true);
  return new RooBreitWigner("psi_3770", "psi_3770", m, *m_chi0_2P, *s_chi0_2P);

}


RooAddPdf* cbVoightMC2(RooRealVar&m ,double minmass, double maxmass, std::string header ="Model" , std::string res = "X", bool fixWidth = false, double fwidth = 0.0){


 double wval; double mval;
  if (res=="X"){
    wval = 0.0317;
    //wval = 0.3;
    mval = 3871.5;
  }
  else {
    wval = 0.304;
    mval = 3686.11;
  }

  if (fixWidth == true){
    wval = fwidth;
  }

  RooRealVar* mx = new RooRealVar(varname("mx",header).c_str(),varname("mx",header).c_str(), mval, mval-3,mval+3);
  // RooRealVar* s = new RooRealVar(varname("sigma2",header).c_str(),varname("sigma2",header).c_str(),2.79,1,10 ); s->setConstant(true);
 //MC psi// 
  // RooRealVar* s = new RooRealVar(varname("sigma",header).c_str(),varname("sigma",header).c_str(),2.23,0, 5);  s->setConstant(true);
 // data psi
  // RooRealVar* s = new RooRealVar(varname("sigma",header).c_str(),varname("sigma",header).c_str(),2.18,0, 5);  s->setConstant(true);
  // RooRealVar* s = new RooRealVar(varname("sigma",header).c_str(),varname("sigma",header).c_str(),2.33,0, 5);  s->setConstant(true);
 // data X
 // data X
 RooRealVar* s = new RooRealVar(varname("sigma",header).c_str(),varname("sigma",header).c_str(),2.82,0, 5);  s->setConstant(true);
 RooRealVar* width =new RooRealVar(varname("width",header).c_str(), varname("width",header).c_str(), wval, -3, 3); 
 if (fixWidth == true) width->setConstant(true);
 RooVoigtian* voight= new  RooVoigtian(varname("voight",header).c_str(),varname("voight(x,mean,sigma)",header).c_str(),m, *mx , *width,*s); 

 RooRealVar* a = new RooRealVar(varname("a",header).c_str(), varname("a",header).c_str(), 2, 0., 4);a->setConstant(true);
 RooRealVar* n = new RooRealVar(varname("n",header).c_str() ,varname("n",header).c_str(), 1., 0., 10); n->setConstant(true);

 // MC x 
 //RooRealVar* s2 = new RooRealVar(varname("sigma2",header).c_str(),varname("sigma2",header).c_str(),4.79,1,10 ); s2->setConstant(true);
 //MC psi 
 //RooRealVar* s2 = new RooRealVar(varname("sigma2",header).c_str(),varname("sigma2",header).c_str(),3.67,1,10 ); s2->setConstant(true);
 // data psi 
 //RooRealVar* s2 = new RooRealVar(varname("sigma2",header).c_str(),varname("sigma2",header).c_str(),3.88,1,10 ); s2->setConstant(true);
 // data X 
 RooRealVar* s2 = new RooRealVar(varname("sigma2",header).c_str(),varname("sigma2",header).c_str(),4.64,1,10 ); s2->setConstant(true);
 RooCBShape* cb = new RooCBShape(varname("cb",header).c_str(),varname("cb",header).c_str(),m, *mx,*s2,*a,*n);

 RooRealVar* f = new RooRealVar(varname("f",header).c_str(),varname("f",header).c_str(),0.79, 0.,1.); f->setConstant(true);
 RooAddPdf* sigmodel= new RooAddPdf(header.c_str(),header.c_str(),RooArgList(*voight,*cb), *f); 
 return sigmodel;
}


RooRelBreitWigner* bw(RooRealVar&m ,double minmass, double maxmass, std::string header ="Model" , std::string res = "X", bool fixWidth = false, double fwidth = 0.317 ){

  double wval; double mval;
  if (res=="X"){
    wval = 0.317;
    mval = 3871.5;
  }
  else {
    wval = 0.304;
    mval = 3686.11;
  }

  if (fixWidth == true){
    wval = fwidth;
  }

  RooRealVar* mx = new RooRealVar(varname("mx",header).c_str(),varname("mx",header).c_str(), mval, mval-1,mval+1);
  RooRealVar* ma = new RooRealVar(varname("ma",header).c_str(),varname("ma",header).c_str(), 3096.916);
  RooRealVar* mb = new RooRealVar(varname("mb",header).c_str(),varname("mb",header).c_str(),507);
  RooRealVar* width =new RooRealVar(varname("width",header).c_str(), varname("width",header).c_str(), wval,  0.01, 5*wval); 
  RooRealVar* radius =new RooRealVar(varname("radius",header).c_str(), varname("radius",header).c_str(),3e-3); 
  RooRealVar* spin =new RooRealVar(varname("spin",header).c_str(), varname("spin",header).c_str(), 0); 
  if (fixWidth == true) width->setConstant(true); //mx.setConstant(true);
  return new RooRelBreitWigner(header.c_str(),header.c_str(), m, *mx, *width,*radius,*ma,*mb,*spin);

}

RooRelBreitWigner* bw2(RooRealVar&m ,double minmass, double maxmass, RooRealVar* mx,  RooRealVar* width,  std::string header ="Model" , std::string res = "X", bool fixWidth = true){

 
  width->setConstant(fixWidth);

  RooRealVar* ma = new RooRealVar(varname("ma",header).c_str(),varname("ma",header).c_str(), 3096.916);
  RooRealVar* mb = new RooRealVar(varname("mb",header).c_str(),varname("mb",header).c_str(),348);

  RooRealVar* radius =new RooRealVar(varname("radius",header).c_str(), varname("radius",header).c_str(),3e-3); 
  RooRealVar* spin =new RooRealVar(varname("spin",header).c_str(), varname("spin",header).c_str(), 0); 
  //if (fixWidth == true) width->setConstant(true); //mx.setConstant(true);
  return new RooRelBreitWigner(header.c_str(),header.c_str(), m, *mx, *width,*radius,*ma,*mb,*spin);

}

RooRelBreitWigner* bwg2(RooRealVar&m ,double minmass, double maxmass, RooAbsReal* mx,  RooRealVar* width,  std::string header ="Model" , double avalue = 518.0e-3, bool fixWidth = true){

 
  width->setConstant(fixWidth);

  RooRealVar* ma = new RooRealVar(varname("ma",header).c_str(),varname("ma",header).c_str(), 3096.916e-3);
  RooRealVar* mb = new RooRealVar(varname("mb",header).c_str(),varname("mb",header).c_str(),avalue);

  RooRealVar* radius =new RooRealVar(varname("radius",header).c_str(), varname("radius",header).c_str(),3e-3); 
  RooRealVar* spin =new RooRealVar(varname("spin",header).c_str(), varname("spin",header).c_str(), 0); 
  //if (fixWidth == true) width->setConstant(true); //mx.setConstant(true);
  return new RooRelBreitWigner(header.c_str(),header.c_str(), m, *mx, *width,*radius,*ma,*mb,*spin);

}



RooRelBreitWigner* bwg(RooRealVar&m ,double minmass, double maxmass, std::string header ="Model" , std::string res = "X", bool fixWidth = false, double fwidth = 0.317e-3 ){

  double wval; double mval;
  if (res=="X"){
    wval = 0.317e-3;
    mval = 3871.5e-3;
  }
  else {
    wval = 0.304e-3;
    mval = 3686.11e-3;
  }

  if (fixWidth == true){
    wval = fwidth;
  }

  RooRealVar* mx = new RooRealVar(varname("mx",header).c_str(),varname("mx",header).c_str(), mval, mval-1,mval+1);
  RooRealVar* ma = new RooRealVar(varname("ma",header).c_str(),varname("ma",header).c_str(), 3.096916);
  RooRealVar* mb = new RooRealVar(varname("mb",header).c_str(),varname("mb",header).c_str(),0.348);
  RooRealVar* width =new RooRealVar(varname("width",header).c_str(), varname("width",header).c_str(), wval,  1e-5, 10*wval); 
  RooRealVar* radius =new RooRealVar(varname("radius",header).c_str(), varname("radius",header).c_str(),3); 
  RooRealVar* spin =new RooRealVar(varname("spin",header).c_str(), varname("spin",header).c_str(), 0); 
  if (fixWidth == true) width->setConstant(true); //mx.setConstant(true);
  return new RooRelBreitWigner(header.c_str(),header.c_str(), m, *mx, *width,*radius,*ma,*mb,*spin);

}


RooChebychev* cheb(RooRealVar&m , std::string header = "cheb",  
                   bool fixA = true, bool fixB = true, bool fixC = true,
                   double aval = 0.01, double bval = 0.01, double cval = 0.01  ){

 RooRealVar* a = new RooRealVar(varname("a",header).c_str(), varname("a",header).c_str(),aval, -0.5,0.8); a->setConstant(fixA); 
 RooRealVar* b= new RooRealVar(varname("b",header).c_str(), varname("b",header).c_str(), bval, -1,1);  b->setConstant(fixB);
 RooRealVar* c=new RooRealVar(varname("c",header).c_str(), varname("c",header).c_str(), cval , -0.1,1); c->setConstant(fixC);
 return new RooChebychev(header.c_str(),header.c_str(),m,RooArgList(*a,*b,*c));
 
}


RooChebychev* cheb2(RooRealVar&m , std::string header = "cheb",  
                   bool fixA = true, bool fixB = true,
		   double aval = 0.05, double bval = -0.01, 
		    double eamin = -1.5, double eamax = 1.5,
		    double ebmin = -1.5, double ebmax = 1.5
   ){

 RooRealVar* a = new RooRealVar(varname("a",header).c_str(), varname("a",header).c_str(),aval, eamin,eamax); a->setConstant(fixA); 
 RooRealVar* b= new RooRealVar(varname("b",header).c_str(), varname("b",header).c_str(), bval, ebmin,ebmax);  b->setConstant(fixB);

 return new RooChebychev(header.c_str(),header.c_str(),m,RooArgList(*a,*b));
 
}


RooChebychev* cheb1(RooRealVar&m , std::string header = "cheb",  
                   bool fixA = true, 
		   double aval = 0.1,
		   double eamin = -0.1, double eamax = 0.5
                 
   ){

 RooRealVar* a = new RooRealVar(varname("a",header).c_str(), varname("a",header).c_str(),aval, eamin,eamax); a->setConstant(fixA); 
 
 return new RooChebychev(header.c_str(),header.c_str(),m,RooArgList(*a));
 
}



RooDstD0BG* dstar(RooRealVar&m , std::string header = "dstar",  
                  bool fixA = true, bool fixB = true, bool fixC = true, double aval = 0.61, double bval = -1.12, double cval = 3){


 RooRealVar* dm = new RooRealVar(varname("dm",header).c_str(), varname("dm",header).c_str(),3.376); 
 RooRealVar* a = new RooRealVar(varname("a",header).c_str(), varname("a",header).c_str(), aval,-5,5); a->setConstant(fixA); 
 RooRealVar* b = new RooRealVar(varname("b",header).c_str(), varname("b",header).c_str(), bval, -10,10); b->setConstant(fixB);
 RooRealVar* c = new RooRealVar(varname("c",header).c_str(), varname("c",header).c_str(), cval , 2.5,4.5); c->setConstant(fixC);

 return new RooDstD0BG(header.c_str(),header.c_str(),m,*dm,*c,*a,*b);

}

RooBackPdf* backPdf(RooRealVar&m , std::string header = "bg", 
                    double c0val = 2.84, double c1val = -0.02,
                    double c2val = 1e-5, 
                    bool fixC0 = true, bool fixC1 = true, bool fixC2 = true){

 RooRealVar* c0 = new RooRealVar(varname("c0",header).c_str(), varname("c0",header).c_str(),c0val ,-20,20); c0->setConstant(fixC0);
 RooRealVar* c1 = new RooRealVar(varname("c1",header).c_str(), varname("c1",header).c_str(),  c1val ,-25.5,25.5);c1->setConstant(fixC1);
 RooRealVar* c2 = new RooRealVar(varname("c2",header).c_str(), varname("c2",header).c_str(),c2val , -0.1,0.1); c2->setConstant(fixC2);
 RooRealVar* mr = new RooRealVar(varname("dm",header).c_str(), varname("dm",header).c_str(), 3376.0e-3, 2000,4000); mr->setConstant(true);
 return new RooBackPdf(header.c_str(),header.c_str(),m, *mr, *c0, *c1, *c2);
}


RooNumConvPdf* convolutionPdf(RooRealVar&m, RooRelBreitWigner* bWigner, std::string name = "ConvPdf" , double sigmaVal = 2.8, bool fixSigma = false){

  RooRealVar* meanres = new RooRealVar(varname("meanres", name).c_str(),varname("meanres", name).c_str(), 0);// mx.setConstant(true); 
  RooRealVar* s = new RooRealVar(varname("sigma",name).c_str(),varname("sigma",name).c_str(),sigmaVal, 0,20); if (fixSigma) s->setConstant(true);
  RooGaussian* gauss = new RooGaussian(varname("gauss",name).c_str(),varname("gauss",name).c_str() ,m, *meanres,*s);
  return new RooNumConvPdf(name.c_str(),name.c_str(),  m, *bWigner,*gauss); 
 
}


RooFormulaVar* tailModel(RooRealVar& sigma, std::string name= "atail",
                         std::string par0 = "2.61", std::string par1 = "1.17", 
                         std::string par2 = "0.585" ){

 std::string ss = sigma.GetName() ; 
 std::string formString = par0+"*pow("+par2+"*" + ss + ","+par1+")/(1 + pow("+par2 + "*" + ss + "," +par1+"))";
 return new RooFormulaVar(name.c_str(), name.c_str(),  formString.c_str(), sigma);
}



RooFormulaVar* tailModelScale(RooRealVar& scale, std::string sigma, 
                              std::string name= "atail", 
                              std::string par0 = "2.61", std::string par1 = "1.17", 
                              std::string par2 = "0.585"){

 std::string sString = scale.GetName() ;
 std::string ss = sString + "*" + sigma; 
 std::string formString = par0+"*pow("+par2+"*" + ss + ","+par1+")/(1 + pow("+par2 + "*" + ss + "," +par1+"))";
 return new RooFormulaVar(name.c_str(), name.c_str(),  formString.c_str(), scale);
}

std::string scaledSigma(std::string sigma, std::string name = "scale"){
  //std::cout << "scaled sigma "<< sigma + "*" + name << std::endl;
  return sigma + "*" + name;
}

std::string scaledSigma(double sval, std::string name = "scale"){
  std::stringstream sigma; sigma << sval;
  //std::cout << "scaled sigma "<< sigma.str() + "*" + name << std::endl;
  return sigma.str() + "*" + name;
}

/*
RooFlatte* flat(RooRealVar&m ,double minmass, double maxmass, RooRealVar* mx,  RooRealVar* width, bool fixWidth,  std::string header = "FlatModel"){
 
  width->setConstant(fixWidth);
 
  double a0 = 0.4;
  double a1 = 0.03;

  RooRealVar* f0 = new RooRealVar(varname("f0",header).c_str(),varname("f0",header).c_str(), a0/(a0+a1));
  RooRealVar* f1 = new RooRealVar(varname("f1",header).c_str(),varname("f1",header).c_str(), a1/(a0+a1));

  RooFormulaVar*  g0 =  new RooFormulaVar("g0", "g0", "width *f0",RooArgSet(*width,*f0)); 
  RooFormulaVar*  g1 =  new RooFormulaVar("g1", "g1", "width *f1",RooArgSet(*width,*f1)); 

  RooRealVar* m0a = new RooRealVar(varname("m0a",header).c_str(),varname("m0a",header).c_str(), 3096.96e-3);
  RooRealVar* m0b = new RooRealVar(varname("m0b",header).c_str(),varname("m0b",header).c_str(), 0.717);
  RooRealVar* m1a = new RooRealVar(varname("m1a",header).c_str(),varname("m1a",header).c_str(), 1864e-3);
  RooRealVar* m1b = new RooRealVar(varname("m1b",header).c_str(),varname("m1b",header).c_str(), 2007e-3);


  return new RooFlatte(header.c_str(),header.c_str(), m, *mx, *g0, *m0a, *m0b, *g1, *m1a, *m1b  );

}
*/

X3872::Hanhard* hanhard(RooRealVar&m ,double minmass, double maxmass, RooRealVar* mx, RooRealVar* g, RooRealVar* f_rho,  RooRealVar* f_omega , RooRealVar* f_g0,  std::string header = "HanhardModel"){

  return new X3872::Hanhard( header.c_str(), header.c_str(), m, *mx, *g, *f_rho, *f_omega, *f_g0 );
}

RooFFTConvPdf* fftpdf(std::string name, RooRealVar& m, RooAbsPdf* bw ,RooAbsPdf* pdf){
 RooFFTConvPdf* tfun1 = new RooFFTConvPdf(name.c_str(),name.c_str(),  m, *bw,*pdf);
 tfun1->setBufferStrategy(RooFFTConvPdf::Extend);
 tfun1->setBufferFraction(0.5);
 return tfun1;
}
