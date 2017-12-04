/***************************************************************************** 
 * Project: RooFit                                                           * 
 *                                                                           * 
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/ 

// Your description goes here... 

#include "Riostream.h" 

#include "RooATanBack.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include <math.h> 
#include "TMath.h" 

ClassImp(RooATanBack) 

 RooATanBack::RooATanBack(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _a,
                        RooAbsReal& _b,
                        RooAbsReal& _c,
                        RooAbsReal& _mr) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   a("a","a",this,_a),
   b("b","b",this,_b),
   c("c","c",this,_c),
   mr("mr","mr",this,_mr)
 { 
 } 


 RooATanBack::RooATanBack(const RooATanBack& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   a("a",this,other.a),
   b("b",this,other.b),
   c("c",this,other.c),
   mr("mr",this,other.mr)
 { 
 } 



 Double_t RooATanBack::evaluate() const 
 { 
   // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE 
   return std::atan((x - mr)/c) + b*((x/mr)-1) + a; 
 } 



Int_t RooATanBack::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* ) const 
{
  if (matchArgs(allVars, analVars, x)) return 1;
  return 0;
}

Double_t  RooATanBack::analyticalIntegral(Int_t code, const char* rangeName) const 
{
  R__ASSERT(code==1) ;

  double sum = theNorm(x.max(rangeName)) -  theNorm(x.min(rangeName));    
  return sum;  
  
}

double RooATanBack::theNorm(const double value) const{

  double relmass = (value-mr)/c;
  double term1= mr*std::atan(-relmass) + value*atan(relmass);
  double term2 = -0.5*c*log(pow(c,2) + pow(mr-value,2));
  double term3 = b*value*(0.5*(value/mr) - value);
  double term4 = a*value;
  return term1 + term2 + term3 + term4;
} 

