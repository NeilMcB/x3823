/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
 * Copyright (c) 2000-2007, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

// file: RooJohnsonSU.h
// author: Maurizio Martinelli (Nikhef)
// email: maurizio.martinelli@cern.ch
// description: this distribution can have highly asymmetric tails and is 
//   therefore helpful in fitting the mass difference between D* and D0 mass.

#ifndef ROOJOHNSONSU
#define ROOJOHNSONSU

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
 
class RooJohnsonSU : public RooAbsPdf {
public:
  RooJohnsonSU() {} ; 
  RooJohnsonSU(const char *name, const char *title,
	      RooAbsReal& _x,
	      RooAbsReal& _xMed,
	      RooAbsReal& _sigx,
	      RooAbsReal& _delta,
	       RooAbsReal& _gamma);
  RooJohnsonSU(const RooJohnsonSU& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooJohnsonSU(*this,newname); }
  inline virtual ~RooJohnsonSU() { }

protected:

  RooRealProxy x ;
  RooRealProxy xMed ;
  RooRealProxy sigx ;
  RooRealProxy delta ;
  RooRealProxy gamma ;
  
  Double_t evaluate() const ;

private:

  ClassDef(RooJohnsonSU,1) // Your description goes here...
};
 
#endif
