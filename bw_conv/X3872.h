#ifndef OSTAP_X3872_H 
#define OSTAP_X3872_H 1
// ============================================================================
// Include files 
// ============================================================================
#include <cmath>
#include <complex>
// ============================================================================

// ============================================================================
// ROOT 
// ============================================================================
#include "Math/ChebyshevApprox.h"
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooRealProxy.h"
// ============================================================================
namespace X3872
{
  // ==========================================================================
  /** @class Hanhard
   *  Flatte parametrization of X(3872) shape 
   *  @see Phys.Rev. D76 (2007) 034007  
   *  @see http://arxiv.org/abs/arXiv:0704.0605 
   *  @see http://dx.doi.org/10.1103/PhysRevD.76.034007
   *
   *  Gamma_0 term is added according to 
   *  Yu.S. Kalashnikova, A.V. Nefediev, "Nature of X(3872) from data",
   *  Phys.Rev. D80 (2009) 074004 
   *  http://arxiv.org/abs/arXiv:0907.4901
   *  http://dx.doi.org/10.1103/PhysRevD.80.074004
   */
  class Hanhard : public RooAbsPdf 
  {
  public:
    // ======================================================================
    ClassDef(X3872::Hanhard, 1) ;
    // ======================================================================
  public:
    // ======================================================================
    /// constructor from all parameters 
    Hanhard ( const char*          name      , 
              const char*          title     ,
              RooAbsReal&          x         ,
              RooAbsReal&          m0        ,   // the mass 
              RooAbsReal&          g         ,   // couplings to D
              RooAbsReal&          f_rho     ,   // coupling  to J/psi pi pi
              RooAbsReal&          f_omega   ,   // coupling  to J/psi pi pi pi  
              RooAbsReal&          gamma_0   ) ; // other decays: psi' gamma, ..
    /// "copy" constructor 
    Hanhard ( const Hanhard& , const char* name = 0 ) ;
    /// virtual destructor 
    virtual ~Hanhard(){} ;
    /// clone 
    virtual Hanhard* clone ( const char* name ) const ; 
    // =======================================================================
  public: // some fake functionality
    // =======================================================================
    // fake default contructor, needed just for proper (de)serialization 
    Hanhard() ;
    // =======================================================================
  public:
    // =======================================================================
    // the actual evaluation of function 
    virtual Double_t evaluate() const ;
    // =======================================================================
  public:
    // ========================================================================
    /// get Gamma(J/psi rho)
    double gamma_rho   () const ;
    /// get Gamma(J/psi omega)
    double gamma_omega () const ;    
    // ========================================================================
  protected:
    // ========================================================================
    /// get D(M) factor 
    std::complex<double> Dfunc ( const double E    ,
                                 const double g2pi , 
                                 const double g3pi , 
                                 const double g0   ) const ;
    // ========================================================================
  public:
    // ========================================================================
    /// the shape in J/psi pipi    channel 
    double gamma_rho   ( const double M ) const ;
    /// the shape in J/psi pipipi0 channel 
    double gamma_omega ( const double M ) const ;
    // ========================================================================
  public:
    // ========================================================================
    /// integrated shape in J/psi pipi    channel (modulo constant factors)
    double R_rho   () const ;
    /// integrated shape in J/psi pipipi0 channel (modulo constant factors)
    double R_omega () const ;    
    // ========================================================================
  protected:
    // =======================================================================
    RooRealProxy m_x      ;
    RooRealProxy m_m0     ;
    RooRealProxy m_g      ; // couplings to D
    RooRealProxy m_frho   ; // coupling  to J/psi pi pi 
    RooRealProxy m_fomega ; // coupling  to J/psi pi pi pi 
    RooRealProxy m_gamma0 ; // effective width for other decays
    // ========================================================================
  private:
    // ========================================================================
    ROOT::Math::ChebyshevApprox m_gamma_2pi  ;
    ROOT::Math::ChebyshevApprox m_gamma_3pi  ;
    // ========================================================================
  } ; 
  // ==========================================================================
}
// ============================================================================
#endif // OSTAP_X3872_H
// ============================================================================
// The END 
// ============================================================================
