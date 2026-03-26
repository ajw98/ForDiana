/***************************************************************
*
* Base class for energy model
* 
* DRT -- Wed, 11 Apr 2017
*
****************************************************************/

#include "Model_Energy_Base.h"

Model_Energy_Base :: Model_Energy_Base( Operators_Base & opobj, int ncomp, int nicomp ):
_OpObj( &opobj ),
_Ncomp( ncomp ),
_NIcomp( nicomp )
{ // {{{
} // }}}

Model_Energy_Base :: ~Model_Energy_Base()
{ // {{{
} // }}}

void Model_Energy_Base :: Test_Energy_Model()
{ // {{{

  SmartFieldVec phi(_NIcomp);
  SmartFieldVec mu(_NIcomp);
  SmartFieldVec mu_lin_exp(_NIcomp);
  SmartFieldVec mu_lin_imp(_NIcomp);
  SmartFieldOpMat mu_lin_op(_NIcomp, _NIcomp);

  phi[0].readfield("phi1.in");
  phi[1].readfield("phi2.in");
  phi.fft();

  for (int i=0; i<_NIcomp; i++)
  {
    std::cout << "phi" << i << " \n";
    std::cout << phi[i].getelement(0) << " \n";
  }
  
  calc_mu( phi, mu );

  for (int i=0; i<_NIcomp; i++)
  {
    std::cout << "mu " << i << " \n";
    std::cout << mu[i].getelement(0) << " \n";
  }

  calc_mu_lin_exp( phi, mu_lin_exp );

  for (int i=0; i<_NIcomp; i++)
  {
    std::cout << "mu_lin_exp " << i << " \n";
    std::cout << mu_lin_exp[i].getelement(0) << " \n";
  }

  calc_mu_lin_imp( phi, mu_lin_op );
  mu_lin_op.OpDot( phi, mu_lin_imp );

  for (int i=0; i<_NIcomp; i++)
  {
    std::cout << "mu_lin_imp " << i << " \n";
    std::cout << mu_lin_imp[i].getelement(0) << " \n";
  }

} // }}}

