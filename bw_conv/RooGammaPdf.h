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

#ifndef ROOGAMMAPDF
#define ROOGAMMAPDF

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
 
class RooGammaPdf : public RooAbsPdf {
public:
  RooGammaPdf(const char *name, const char *title,
	      RooAbsReal& _x,
	      RooAbsReal& _c,
	      RooAbsReal& _s,
              RooAbsReal& _loc);
  RooGammaPdf(const RooGammaPdf& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooGammaPdf(*this,newname); }
  inline virtual ~RooGammaPdf() { }

protected:

  RooRealProxy x ;
  RooRealProxy c ;
  RooRealProxy s ;
  RooRealProxy loc;
  
  Double_t evaluate() const ;

private:

  ClassDef(RooGammaPdf,0) // Your description goes here...
};
 
#endif
