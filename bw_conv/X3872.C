// $Id:$ 
// ============================================================================
// Include files 
// ============================================================================
#include <cmath>
#include <complex>
// ============================================================================

// ============================================================================
// ROOT 
// ============================================================================
#include "Math/Integrator.h"
#include "Math/GSLIntegrator.h"
#include "Math/ChebyshevApprox.h"
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooRealProxy.h"
// ============================================================================
// Local
// ============================================================================
#include "X3872.h"
#include "X3872_constants.h"
// ============================================================================
namespace X3872
{
  // ==========================================================================
  // units 
  //const double GeV   =      1 ;  // unit
  //const double MeV   = 0.001 * GeV ;
  
  // rho0
  // const double m_rho   = 775.49    * MeV ;
  //const double g_rho   = 149.1     * MeV ;
  
  // omega 
  //const double m_omega = 782.65    * MeV ;
  //  const double g_omega =   8.49    * MeV ;
  
  // D0 & D*0 
  //const double m_D0     = 1864.83  * MeV ;
  //  const double m_D0star = 2006.85  * MeV ;
  

  // D0 & D*0 
  ///const double m_Dp     = 1869.58  * MeV ;
  //const double m_Dpstar = 2010.26  * MeV ;
  
  // J/psi 
  //const double m_Jpsi   = 3096.900 * MeV ;

  // pion 
  //const double m_pion   =  139.570 * MeV ;
  // const double m_pi0    =  134.977 * MeV ;

  // reduced masses
  //const double mu1      = m_D0 * m_D0star / ( m_D0 + m_D0star ) ;
  //const double mu2      = m_Dp * m_Dpstar / ( m_Dp + m_Dpstar ) ;
  
  // delta 
  //const double delta    = ( m_Dp + m_Dpstar ) - ( m_D0 + m_D0star ) ;
  // ==========================================================================
  
  // ==========================================================================
  /** momentum of di-pion/tri-pion system in X(3872) rest frame 
   *  @param m  mass of dipion/tripion system
   *  @param M  mass of  X3872 resonance 
   */
  inline double  q ( const double m , 
                     const double M ) 
  {
    const double m_1 = m + m_Jpsi ;
    if ( m_1 >= M ) {  return 0 ; }          // RETURN
    //
    const double M2  = M * M ;
    const double m_2 = m - m_Jpsi ;
    const double r   = ( M2 - m_1 * m_1 ) * ( M2 - m_2 * m_2 ) ;
    //
    return r <= 0 ? 0 : 0.5 * std::sqrt ( r ) / M ;
    //return 1 ;
  }
  
  /// rho0 shape 
  inline double rho_line ( const double m , void* params ) 
  {
    const double M  =  *( (double*) params ) ;
    const double qm  = q ( m , M ) ;
    if ( qm <= 0 ) { return 0 ; }              // REUTRN
    //
    const double dm = m - m_rho ;
    const double fun = qm * g_rho / ( dm * dm + 0.25 * g_rho * g_rho ) ;
    return fun / ( 2 * M_PI ) ;  
  }
  
  /// rho0 shape 
  inline double omega_line ( const double m , void* params ) 
  {
    const double M  =  *( (double*) params ) ;
    const double qm  = q ( m , M ) ;
    if ( qm <= 0 ) { return 0 ; }              // RETURN
    //
    const double dm  = m - m_omega ;
    const double fun = qm * g_omega / ( dm * dm + 0.25 * g_omega * g_omega ) ;
    return fun / ( 2 * M_PI ) ;  
  }
  // ==========================================================================  
  /// integrator 
  ROOT::Math::GSLIntegrator  s_integrator 
  ( ROOT::Math::Integration::kADAPTIVE , 
    ROOT::Math::Integration::kGAUSS31  , 1.e-6 , 1.e-4 ) ;
  // ==========================================================================  
  /// gamma_2pi 
  double gamma_2pi ( double M , void* p = 0 ) 
  {
    void* M_void = &M ;
    return 
      M - m_Jpsi <= 2*m_pion ? 0 : 
      s_integrator.Integral ( &X3872::rho_line   , M_void , 2*m_pion , M - m_Jpsi ) ;
  }
  // ==========================================================================  
  /// gamma_3pi 
  double gamma_3pi ( double M , void* p = 0 ) 
  {
    void* M_void = &M ;
    return  M - m_Jpsi <= 2*m_pion + m_pi0 ? 0 : 
      s_integrator.Integral ( &X3872::omega_line   , M_void , 2*m_pion , M - m_Jpsi ) ;
  }
  // ==========================================================================
  double diff_R_rho   ( double M , void* p ) 
  {
    const X3872::Hanhard* f = (const X3872::Hanhard*) p ;
    return f->gamma_rho   ( M ) ;
  }  
  double diff_R_omega ( double M , void* p ) 
  {
    const X3872::Hanhard* f = (const X3872::Hanhard*) p ;
    return f->gamma_omega ( M ) ;
  }
}

// ============================================================================
ClassImp(X3872::Hanhard);
// ============================================================================
// constructor from all parameters 
// ============================================================================
X3872::Hanhard::Hanhard 
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         ,
  RooAbsReal&          m0        ,   // the mass 
  RooAbsReal&          g         ,   // couplings to D
  RooAbsReal&          f_rho     ,   // coupling  to J/psi pi pi
  RooAbsReal&          f_omega   ,   // coupling  to J/psi pi pi pi  
  RooAbsReal&          gamma_0   ) // pther decays; psi' gamma, ...
  //
  : RooAbsPdf( name , title ) 
  , m_x      ( "x"      , "Observable"          , this , x       ) 
  , m_m0     ( "m0"     , "Peak"                , this , m0      ) 
  , m_g      ( "g"      , "Coupling to D"       , this , g       )
  , m_frho   ( "frho"   , "Coupling to rho"     , this , f_rho   )
  , m_fomega ( "fomega" , "Coupling to omega"   , this , f_omega )
  , m_gamma0 ( "gamma0" , "Other decays"        , this , gamma_0 )
    //
  , m_gamma_2pi  (  &X3872::gamma_2pi , 0 , 3.8 * GeV , 4.0 * GeV ,  20 ) 
  , m_gamma_3pi  (  &X3872::gamma_3pi , 0 , 3.8 * GeV , 4.0 * GeV , 100 ) 
{}
// ============================================================================
// "copy" constructor 
// ============================================================================
X3872::Hanhard::Hanhard
( const X3872::Hanhard& right , 
  const char*           name  ) 
  : RooAbsPdf ( right , name ) 
    //
  , m_x      ( "x"       , this , right.m_x      ) 
  , m_m0     ( "m0"      , this , right.m_m0     ) 
  , m_g      ( "g"       , this , right.m_g      ) 
  , m_frho   ( "frho"    , this , right.m_frho   ) 
  , m_fomega ( "fomega"  , this , right.m_fomega ) 
  , m_gamma0 ( "gamma0"  , this , right.m_gamma0 ) 
    //
  , m_gamma_2pi  ( &X3872::gamma_2pi , 0 , 3.8 * GeV , 4.0 * GeV ,  20 ) 
  , m_gamma_3pi  ( &X3872::gamma_3pi , 0 , 3.8 * GeV , 4.0 * GeV , 100 ) 
{}
// ============================================================================
// fake default contructor, needed just for proper (de)serialization 
// ============================================================================
X3872::Hanhard::Hanhard() 
  : RooAbsPdf    () 
  , m_gamma_2pi  ( &X3872::gamma_2pi , 0 , 3.8 * GeV , 4.0 * GeV ,  20 ) 
  , m_gamma_3pi  ( &X3872::gamma_3pi , 0 , 3.8 * GeV , 4.0 * GeV , 100 ) 
{}
// ============================================================================
// clone it!
// ============================================================================
X3872::Hanhard*
X3872::Hanhard::clone( const char* name ) const 
{ return new X3872::Hanhard ( *this , name ) ; }
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t X3872::Hanhard::evaluate() const 
{
  // the mass of X(3872) resonance 
  const double M = m_x ;
  return gamma_rho ( M ) ;
}
// ============================================================================
// get D(M) 
// ============================================================================
std::complex<double> X3872::Hanhard::Dfunc 
( const double M    , 
  const double g2pi , 
  const double g3pi , 
  const double g0   ) const 
{
  static const std::complex<double> im (0,1) ; // i
  //
  std::complex<double> D = ( M - m_m0 ) + 0.5 * im * ( g2pi * m_frho  + g3pi * m_fomega + std::abs ( g0 ) ) ;
  //
  const double E = M - m_D0 - m_D0star ;
  if      ( E <= 0     )
  {
    const double kappa1 = std::sqrt( -2 * mu1 *   E            ) ;
    const double kappa2 = std::sqrt( -2 * mu2 * ( E - delta )  ) ; 
    const double g      = m_g ;
    D +=  0.5 * (-kappa1 -kappa2 ) * g ;  
  }
  else if ( E <= delta ) 
  {
    const double k1     = std::sqrt(  2 * mu1 *   E            ) ;
    const double kappa2 = std::sqrt( -2 * mu2 * ( E - delta )  ) ;
    const double g      = m_g ;
    D +=  0.5 * ( im * k1 - kappa2 ) * g ;  
  }
  else 
  {
    const double k1     = std::sqrt(  2 * mu1 *   E            ) ;
    const double k2     = std::sqrt(  2 * mu2 * ( E - delta )  ) ;
    const double g      = m_g ;
    D +=  0.5 * im  * ( k1 + k2 ) * g  ;  
  }
  return D ;
}
// ============================================================================
double X3872::Hanhard::gamma_rho   ( const double M ) const 
{
  if ( M <= m_Jpsi + 2 * m_pion ) { return 0 ; }  // RETURN
  const double E = M - m_D0 - m_D0star ;
  
  // direct integration, time consuming
  // const double g2pi = gamma_2pi ( M ) ; 
  // const double g3pi = gamma_3pi ( M ) ; 
  //
  // Chebyshev approximation, much faster
  const double g2pi = m_gamma_2pi ( M ) ; 
  const double g3pi = m_gamma_3pi ( M ) ; 
  const double g0   = m_gamma0          ;
  //
  const std::complex<double> D = Dfunc ( M , g2pi , g3pi , g0 ) ;
  //
  return g2pi / std::norm ( D ) * m_frho ;
}
// ============================================================================
double X3872::Hanhard::gamma_omega ( const double M ) const 
{
  if ( M <= m_Jpsi + 2 * m_pion ) { return 0 ; }  // RETURN
  const double E = M - m_D0 - m_D0star ;
  
  // direct integration, time consuming
  // const double g2pi = gamma_2pi ( M ) ; 
  // const double g3pi = gamma_3pi ( M ) ; 
  //
  // Chebyshev approximation, much faster
  const double g2pi = m_gamma_2pi ( M ) ; 
  const double g3pi = m_gamma_3pi ( M ) ; 
  const double g0   = m_gamma0          ;
  //
  const std::complex<double> D = Dfunc ( M , g2pi , g3pi , g0 ) ;
  //
  return g3pi / std::norm ( D ) * m_fomega ;
}
// ============================================================================
double X3872::Hanhard::R_rho   () const 
{
  const double low  = m_D0 + m_D0star - 20 *MeV ;
  const double high = m_D0 + m_D0star + 20 *MeV ;
  //
  void* params = const_cast<X3872::Hanhard*>(this) ;
  return 
    s_integrator.Integral ( &X3872::diff_R_rho  , params , low , high  ) ;
}
// ============================================================================
double X3872::Hanhard::R_omega () const 
{
  const double low  = m_D0 + m_D0star - 20 *MeV ;
  const double high = m_D0 + m_D0star + 20 *MeV ;
  //
  void* params = const_cast<X3872::Hanhard*>(this) ;
  return 
    s_integrator.Integral ( &X3872::diff_R_omega  , params , low , high  ) ;
}
// ============================================================================


// ============================================================================
// The END 
// ============================================================================
