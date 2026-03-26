/***************************************************************
*
* Base class for free energy models
* 
* DRT -- Wed, 18 Aug 2015
*
****************************************************************/
#ifndef _MODEL_ENERGY_BASE_H
#define _MODEL_ENERGY_BASE_H

#include "global.h"
#include "SmartFieldVec.h"
#include "SmartFieldOpMat.h"
#include "Operators_Base.h"

class Model_Energy_Base
{

  public:
    Model_Energy_Base( Operators_Base & OpObj, int ncomp, int nicomp );
    virtual ~Model_Energy_Base();

    virtual void set_params() = 0;

    // calculate chemical potentials
    // (All inputs/outputs are in Fourier Space!)
    virtual void calc_mu(const SmartFieldVec & phi, SmartFieldVec & mu) = 0; //chemical potential
    virtual void calc_mu_lin_exp(const SmartFieldVec & phi, SmartFieldVec & mu_lin) = 0; //explicit, linear chem. pot.
    virtual void calc_mu_lin_imp(const SmartFieldVec & phi, SmartFieldOpMat & F) = 0; //implicit linear chem. pot.
    virtual void calc_del_mu(const SmartFieldVec & phi, SmartFieldMat & del_mu) = 0;

    // scale for the chemical potential
    // (e.g. used in Model H)
    virtual RealType mu_scale() = 0;

    void Test_Energy_Model();
 
  protected:
    Operators_Base*  _OpObj;

    int _EnergyModelFlag; // Model # Flag
    int _Ncomp; // number of components
    int _NIcomp; // number of *independent* components

};

#endif //_MODEL_ENERGY_BASE_H
