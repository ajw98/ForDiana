/***************************************************************
*
* Class defining a Flory Huggins free energy model
* 
* DRT -- Wed, 18 Aug 2015
*
****************************************************************/

#ifndef _MODEL_ENERGY_LINEAR_H
#define _MODEL_ENERGY_LINEAR_H

#include "global.h"
#include "Field.h"
#include "SmartFieldVec.h"
#include "SmartFieldOpMat.h"
#include "Model_Energy_Base.h"

class Model_Energy_Linear : public Model_Energy_Base
{

  public:
    Model_Energy_Linear( Operators_Base & OpObj, int ncomp, int nicomp );
    ~Model_Energy_Linear();

    void set_params();

    // calculate chemical potentials
    // (All inputs/outputs are in Fourier Space!)
    void calc_mu(const SmartFieldVec & phi, SmartFieldVec & mu); //chemical potential
    void calc_mu_lin_exp(const SmartFieldVec & phi, SmartFieldVec & mu_lin); //explicit linear chem. pot.
    void calc_mu_lin_imp(const SmartFieldVec & phi, SmartFieldOpMat & F); // implicit linear chem pot.
    void calc_del_mu(const SmartFieldVec & phi, SmartFieldMat & del_mu);

    // scale for the chemical potential
    // (e.g. used in Model H)
    inline RealType mu_scale() { return 1.; };

  private:

    // model parameters
    Vec2dFieldType _Hess; // Hessian matrix
    Vec1dReal _Kappa; // Gradient coefficients
  
};

#endif //_MODEL_ENERGY_LINEAR_H
