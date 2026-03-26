/***************************************************************
*
* Class defining a Flory Huggins free energy model with Kappa 3
*
* DRT -- Wed, 18 Aug 2015
* RA  -- Fri, 20 Dec 2019
*
****************************************************************/

#ifndef _MODEL_ENERGY_FHG_K3_H
#define _MODEL_ENERGY_FHG_K3_H

#include "global.h"
#include "SmartFieldVec.h"
#include "SmartFieldMat.h"
#include "SmartFieldOpMat.h"
#include "Model_Energy_Base.h"
#include "Operators_Base.h"

class Model_Energy_FHG_K3 : public Model_Energy_Base
{

  public:
    Model_Energy_FHG_K3( Operators_Base & OpObj, int ncomp, int nicomp );
    ~Model_Energy_FHG_K3();

    void set_params();

    // calculate chemical potentials
    // (All inputs/outputs are in Fourier Space!)
    void calc_mu(const SmartFieldVec & phi, SmartFieldVec & mu); //chemical potential
    void calc_mu_lin_exp(const SmartFieldVec & phi, SmartFieldVec & mu_lin); //explicit linear chem. pot.
    void calc_mu_lin_imp(const SmartFieldVec & phi, SmartFieldOpMat & F); // implicit linear chem pot.
    void calc_del_mu(const SmartFieldVec &phi, SmartFieldMat & del_mu);

    // scale for the chemical potential
    // (e.g. used in Model H)
    inline RealType mu_scale() {return _Nr;};

    // supporting member functions
    void calc_Hessian( const SmartFieldVec & phi, SmartFieldMat & Hess );
    Vec2dFieldType calc_Hmax( const SmartFieldVec & phi );

  private:

    // model parameters
    RealType _Nr; // Molecular weight
    RealType _CReg, _CRegN, _Delta; // Regularization parameters

    Vec1dReal _N; // degree of polymerization
    Vec1dReal _alpha; // relative degree of polymerization to Nr
    Vec2dReal _Chi; // Flory interaction parameters
    Vec2dReal _ChiN; // Flory interaction parameters * Nr
    Vec2dReal _ChiN_ijM; // = ChiN_ij - ChiN_iM - ChiN_jM (used repeatedly)
    Vec2dReal _Kappa; // Gradient coefficients

    // These variables help so we only have to do the calculation once
    // This helps considerably with speed.
    bool _PrecalcMuExp; // have we precalculated variables for explicit linear part?
    bool _PrecalcMuImp; // have we precalculated variables for implicit linear part?
    Vec2dFieldType _Hmax; // Maximum Hessian
    SmartFieldVec _KappaGrad2Phi;
    UInt _Dim;
  
};

#endif //_MODEL_ENERGY_FHG_K3_H
