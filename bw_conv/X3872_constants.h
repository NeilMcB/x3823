#ifndef X3872_X3872_CONSTANTS_H 
#define X3872_X3872_CONSTANTS_H 1
// ============================================================================
// Include files
// ============================================================================
namespace X3872
{
  // ==========================================================================
  // units 
  const double GeV   = 1.000       ;  // unit
  const double MeV   = 0.001 * GeV ;
  
  // rho0
  const double m_rho   = 775.49    * MeV ;
  const double g_rho   = 149.1     * MeV ;
  
  // omega 
  const double m_omega = 782.65    * MeV ;
  const double g_omega =   8.49    * MeV ;
  
  // D0 & D*0 
  const double m_D0     = 1864.83  * MeV ;
  const double m_D0star = 2006.85  * MeV ;
  
  
  // D0 & D*0 
  const double m_Dp     = 1869.61  * MeV ;
  const double m_Dpstar = 2010.27  * MeV ;
  
  // J/psi 
  const double m_Jpsi   = 3096.9 * MeV ;
  
  // pion 
  const double m_pion   =  139.570 * MeV ;
  const double m_pi0    =  134.977 * MeV ;
  
  // reduced masses
  const double mu1      = m_D0 * m_D0star / ( m_D0 + m_D0star ) ;
  const double mu2      = m_Dp * m_Dpstar / ( m_Dp + m_Dpstar ) ;
  
  // delta 
  const double delta    = ( m_Dp + m_Dpstar ) - ( m_D0 + m_D0star ) ;

  // width of D*0
  const double GammaStar = 65.5e-3 *MeV;
  
  // ==========================================================================
} // end of namespace X3872
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // X3872_X3872_CONSTANTS_H
// ============================================================================
