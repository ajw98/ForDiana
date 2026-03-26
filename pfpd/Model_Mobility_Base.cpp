/***************************************************************
*
* Base class for mobility model
* 
* DRT -- Wed, 11 Apr 2017
*
****************************************************************/

#include "Model_Mobility_Base.h"

Model_Mobility_Base :: Model_Mobility_Base( int ncomp, int nicomp, Model_Visc_Base* visc_model ):
_MobilityModelFlag(-1),
_Ncomp( ncomp ),
_NIcomp( nicomp ),
_ViscModel( visc_model)
{ // {{{
} // }}}

Model_Mobility_Base :: ~Model_Mobility_Base()
{ // {{{
} // }}}

Vec2dFieldType Model_Mobility_Base :: get_mobility_max( const SmartFieldMat & Mobility )
{ // {{{

  if (Mobility.getflag_inrealspace() != true)
  {
    std::cout << "*** Error in Model_Mobility_Base::get_mobility_max ***\n";
    std::cout << "*** Mobility must be in Real space ***\n";
  }

  Vec2dFieldType M_max( _NIcomp, Vec1dFieldType( _NIcomp, 0.) );
  SmartField M_norm;
  UInt max_indx;

  Mobility.norm2(M_norm, true); // Using norm^2, which is faster and does not change max c.f. norm()
  max_indx = M_norm.maxidx();
  M_max = Mobility.getelement(max_indx);

  return M_max;

} // }}}

