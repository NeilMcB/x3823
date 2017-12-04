#ifndef X3872_BRAATAN_H 
#define X3872_BRAATAN_H 1
// ============================================================================
// Include files 
// ============================================================================
#include <cmath>
#include <complex>
// ============================================================================
// ROOT 
// ============================================================================
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooRealProxy.h"
// ============================================================================
namespace X3872
{
  // ==========================================================================
  /** @class Braatan
   *  @see arXiv:0907.3167
   */
  class Braatan : public RooAbsPdf 
  {
  public:
    // ========================================================================
    ClassDef(X3872::Braatan,1);
    // ========================================================================
  public:
    // ========================================================================
    Braatan(); 
    Braatan
    ( const char *name   , 
      const char *title  ,
      RooAbsReal& x      ,
      RooAbsReal& gRe    ,
      RooAbsReal& gIm    ,
      RooAbsReal& thres  ) ;
    // ========================================================================
    Braatan ( const Braatan& other, const char* name=0) ;
    virtual Braatan* clone(const char* newname) const ;
    virtual ~Braatan();
    // ========================================================================    
  public:
    // ========================================================================    
    double EX     () const;    
    double GammaX () const;    
    // double fwhm   () const;
    // ========================================================================    
  protected:
    // ========================================================================
    RooRealProxy m_x     ;
    RooRealProxy m_gRe   ;
    RooRealProxy m_gIm   ;
    RooRealProxy m_thres ;
    // ========================================================================    
    Double_t evaluate() const ;
    // ========================================================================    
  private:
    // ========================================================================    
    double mu  ()         const;    
    double fE2 (double E) const;
    // ========================================================================
  };
  // ==========================================================================
} //                                                 the end of namespace X3872 
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // X3872_BRAATAN_H
// ============================================================================
