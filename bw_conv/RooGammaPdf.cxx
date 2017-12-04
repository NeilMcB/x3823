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
 // Your description goes here... 

 #include <iostream> 

 #include "RooGammaPdf.h" 
 #include "RooAbsReal.h" 
 #include "TMath.h"

 ClassImp(RooGammaPdf) 

 RooGammaPdf::RooGammaPdf(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _c,
                	  RooAbsReal& _s,
                        RooAbsReal& _loc)
          :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   c("c","c",this,_c),
   s("s","s",this,_s),
   loc("loc","loc", this,_loc)
 { 
 } 


 RooGammaPdf::RooGammaPdf(const RooGammaPdf& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   c("c",this,other.c),
   s("s",this,other.s),
   loc("loc",this,other.loc)
 { 
 } 



 Double_t RooGammaPdf::evaluate() const 
 { 
   // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE 
   double val = 0;
   if (x>loc) val = TMath::GammaDist(x ,c , loc , s ); 
   return val; 
 } 



