/***************************************************************
*
* Class defining an exponential viscosity model similar to the
* form by Vogel-Fulcher-Tamman-Hesse (VFTH)
* 
* JUG --- 8 November 2018
*
****************************************************************/

#ifndef _MODEL_VISC_VFTH
#define _MODEL_VISC_VFTH

#include "global.h"
#include "SmartField.h"
#include "SmartFieldVec.h"
#include "Model_Visc_Base.h"

class Model_Visc_VFTH : public Model_Visc_Base
{

  public:
    Model_Visc_VFTH( int ncomp, int nicomp );
    ~Model_Visc_VFTH();

    void set_params();
    void set_eta(const SmartFieldVec &phi, SmartField &eta);
    RealType  phi_of_Eta ( const RealType eta );
    void capField ( SmartField &F, const RealType k );

  private:

    //int _ViscModelFlag;
    int _ViscCapFlag; //flag for viscosity cap
    RealType _EtaR; // reference viscosity 
    RealType _EtaS; // solvent/nonsolvent viscosity 
    RealType _Gamma; // exponent coefficient 
    RealType _PhiDagger; // glass transition concentration 
    RealType _EtaCap; // viscosity cap 
};

#endif //_MODEL_VISC_VFTH
