// Include files
// ============================================================================
// local
// ============================================================================
#include "TMath.h"
#include "Braatan.h"
#include "X3872_constants.h"
// ============================================================================
ClassImp(X3872::Braatan) ;
// ============================================================================
X3872::Braatan::Braatan (){}
// ============================================================================
X3872::Braatan::~Braatan(){}
// ============================================================================
X3872::Braatan* X3872::Braatan::clone(const char* newname) const
{ return new X3872::Braatan(*this,newname); }
// ============================================================================
X3872::Braatan::Braatan
( const char *name  , 
  const char *title , 
  RooAbsReal& x     ,
  RooAbsReal& gRe   ,
  RooAbsReal& gIm   ,
  RooAbsReal& thres ) 
  : RooAbsPdf(name,title)
  , m_x     ("x"    ,"x"    ,this , x     )
  , m_gRe   ("gRe"  ,"gRe"  ,this , gRe   )
  , m_gIm   ("gIm"  ,"gIm"  ,this , gIm   )
  , m_thres ("thres","thres",this , thres )
{}
// ============================================================================
X3872::Braatan::Braatan
( const X3872::Braatan& other , 
  const char*           name  ) 
  : RooAbsPdf(other,name)
  , m_x     ( "x"     , this , other.m_x     )
  , m_gRe   ( "gRe"   , this , other.m_gRe   )
  , m_gIm   ( "gIm"   , this , other.m_gIm   )
  , m_thres ( "thres" , this , other.m_thres )
{} 
// ============================================================================
double X3872::Braatan::mu() const { return X3872::mu1; }
// ============================================================================
double X3872::Braatan::fE2 ( double E ) const
{
  const std::complex<double> g(-m_gRe, -m_gIm);
  const std::complex<double> c(-2*mu()*E,-X3872::GammaStar*mu());
  const std::complex<double> sc  = std::sqrt(c);
  const std::complex<double> sum = std::complex<double>(1) /(g + sc);
  return std::norm(sum);
}
// ============================================================================
double X3872::Braatan::EX() const 
{
  double top = std::pow(m_gRe,2) - std::pow(m_gIm, 2);
  return top/(2*mu());
}
// ============================================================================
double X3872::Braatan::GammaX() const
{ return X3872::GammaStar + ((2*m_gRe*m_gIm)/mu()); }
// ============================================================================
// double  RooBraatan::fwhm() const{
//   double term1= 4*gRe*gIm;
//   double term2= 2*mu()*Braatan::GammaStar;
//   double diff = gRe*gRe - gIm*gIm;
//   double diff3 =  3*gRe*gRe - gIm*gIm;
//   double term3 = std::pow(gRe/gIm,3) * (diff/diff3) * mu()*mu()*Braatan::GammaStar*Braatan::GammaStar;  
//   return (term1 + term2 + term3)/(2*mu());
// }
// ============================================================================
Double_t X3872::Braatan::evaluate() const 
{ 
  // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE
  const double val = m_x -  m_thres ;
  const double f = fE2(val);
  const double top = mu()*mu()*GammaX()*f;
  const double bottom = 2*TMath::Pi()*(std::pow(m_gRe,2) + std::pow(m_gIm, 2));
  return top/bottom ; 
} 
// ============================================================================
//  The END 
// ============================================================================


