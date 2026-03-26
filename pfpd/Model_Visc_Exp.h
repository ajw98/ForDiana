/***************************************************************
*
* Class defining a sigmoidal viscosity model
* 
* DRT -- Wed, 18 Aug 2015
*
****************************************************************/

#ifndef _MODEL_VISC_EXP_H
#define _MODEL_VISC_EXP_H

#include "global.h"
#include "SmartField.h"
#include "SmartFieldVec.h"
#include "Model_Visc_Base.h"

class Model_Visc_Exp : public Model_Visc_Base
{

  public:
    Model_Visc_Exp( int ncomp, int nicomp );
    ~Model_Visc_Exp();

    void set_params();
    void set_eta(const SmartFieldVec &phi, SmartField &eta);

  private:

    // sigmoidal function for nonlinear viscosity
    void f_exp(SmartField &phi, RealType phiT, RealType w);

    RealType _EtaR; // viscosity of ref
    Vec1dReal _Eta; // viscosity of all Ncomp
    Vec1dReal _PhiT; // threshold concentrations
    Vec1dReal _W; // transition width
    Vec1dInt _ExpOn; // ==1 for comp i to be a sigmoid, ==0, linear
  
};

#endif //_MODEL_VISC_EXP_H
