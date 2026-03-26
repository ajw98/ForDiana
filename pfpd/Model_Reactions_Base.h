/***************************************************************
*
* Base class for reaction models
* 
* MDT -- Wed, 2 Jan 2018
*
****************************************************************/
#ifndef _MODEL_REACTIONS_BASE_H
#define _MODEL_REACTIONS_BASE_H

#include "global.h"
#include "SmartFieldVec.h"
#include "SmartFieldOpMat.h"
#include "Field.h"

class Model_Reactions_Base
{

  public:
    Model_Reactions_Base( int ncomp, int nicomp );
    virtual ~Model_Reactions_Base();

    virtual void set_params() = 0;

    // calculate chemical reactions
    // (All inputs/outputs are in Fourier Space!)
    virtual void calc_reactions(const SmartFieldVec & phi, SmartFieldVec & mu) = 0; //chemical potential
    virtual void calc_reactions_lin_exp(const SmartFieldVec & phi, SmartFieldVec & mu_lin) = 0; //explicit, linear chem. pot.
    virtual void calc_reactions_lin_imp(const SmartFieldVec & phi, SmartFieldOpMat & F) = 0; //implicit linear chem. pot.

    // scale for the chemical potential
    // (e.g. used in Model H)
    virtual RealType rxn_scale() = 0;

    int GetReactionsFlag(){ return _ReactionsModelFlag; };

    void Test_Reactions_Model();
  
  protected:
    int _ReactionsModelFlag; // Model # Flag
    int _Ncomp; // number of components
    int _NIcomp; // number of *independent* components

};

#endif //_MODEL_REACTIONS_BASE_H
