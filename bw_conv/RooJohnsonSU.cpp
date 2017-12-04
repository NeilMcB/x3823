/***************************************************************************** 
 * Project: RooFit                                                           * 
 *                                                                           * 
 * Copyright (c) 2000-2005, Regents of the University of California          * 
 *                          and Stanford University. All rights reserved.    * 
 *                                                                           * 
 * Redistribution and use in source and binary forms,                        * 
 * with or without modification, are permitted according to the terms        * 
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             * 
 *****************************************************************************/ 

// -- CLASS DESCRIPTION [PDF] -- 
// file: RooJohnsonSU.cxx
// author: Maurizio Martinelli (Nikhef)
// email: maurizio.martinelli@cern.ch
// description: this distribution can have highly asymmetric tails and is 
//   therefore helpful in fitting the mass difference between D* and D0 mass.
// 

#include "Riostream.h" 

#include "RooJohnsonSU.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h"

#include "TMath.h"

ClassImp(RooJohnsonSU) 
  
RooJohnsonSU::RooJohnsonSU(const char *name, const char *title, 
			   RooAbsReal& _x,
			   RooAbsReal& _xMed,
			   RooAbsReal& _sigx,
			   RooAbsReal& _delta,
			   RooAbsReal& _gamma) :
RooAbsPdf(name,title), 
  x("x","x",this,_x),
  xMed("xMed","xMed",this,_xMed),
  sigx("sigx","sigx",this,_sigx),
  delta("delta","delta",this,_delta),
  gamma("gamma","gamma",this,_gamma)
{ 
} 


RooJohnsonSU::RooJohnsonSU(const RooJohnsonSU& other, const char* name) :  
  RooAbsPdf(other,name), 
  x("x",this,other.x),
  xMed("xMed",this,other.xMed),
  sigx("sigx",this,other.sigx),
  delta("delta",this,other.delta),
  gamma("gamma",this,other.gamma)
{ 
} 



Double_t RooJohnsonSU::evaluate() const 
{ 
  // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE 
  double residual = (x-xMed)/sigx;
  double JSU = delta*TMath::Exp(-0.5*TMath::Power(gamma+delta*TMath::ASinH(residual),2))/
    (TMath::Sqrt(1+TMath::Power(residual,2))*sigx);

  double val= JSU;

  return val ; 
} 



