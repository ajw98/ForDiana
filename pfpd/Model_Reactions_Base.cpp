/***************************************************************
*
* Base class for energy model
* 
* DRT -- Wed, 11 Apr 2017
*
****************************************************************/

#include "Model_Reactions_Base.h"

Model_Reactions_Base :: Model_Reactions_Base( int ncomp, int nicomp ):
_Ncomp( ncomp ),
_NIcomp( nicomp )
{ // {{{
} // }}}

Model_Reactions_Base :: ~Model_Reactions_Base()
{ // {{{
} // }}}

void Model_Reactions_Base :: Test_Reactions_Model()
{ // {{{

  SmartFieldVec phi(_NIcomp);
  SmartFieldVec rxn(_NIcomp);
  SmartFieldVec rxn_lin_exp(_NIcomp);
  SmartFieldVec rxn_lin_imp(_NIcomp);
  SmartFieldOpMat Rxn_lin_op(_NIcomp, _NIcomp);

  phi[0].readfield("phi1.in");
  phi[1].readfield("phi2.in");
  phi[2].readfield("phi3.in");
  phi.fft();

  for (int i=0; i<_NIcomp; i++)
  {
    std::cout << "phi" << i << " \n";
    std::cout << phi[i].getelement(0) << " \n";
  }
  
  calc_reactions( phi, rxn );

  for (int i=0; i<_NIcomp; i++)
  {
    std::cout << "rxn " << i << " \n";
    std::cout << rxn[i].getelement(0) << " \n";
  }

  calc_reactions_lin_exp( phi, rxn_lin_exp );

  for (int i=0; i<_NIcomp; i++)
  {
    std::cout << "rxn_lin_exp " << i << " \n";
    std::cout << rxn_lin_exp[i].getelement(0) << " \n";
  }

  calc_reactions_lin_imp( phi, Rxn_lin_op );
  Rxn_lin_op.OpDot( phi, rxn_lin_imp );

  for (int i=0; i<_NIcomp; i++)
  {
    std::cout << "rxn_lin_imp " << i << " \n";
    std::cout << rxn_lin_imp[i].getelement(0) << " \n";
  }

} // }}}

