/***************************************************************
*
* Class defining reaction terms
* 
* DRT -- Wed, 18 Aug 2015
*
****************************************************************/

#ifndef _MODEL_REACTIONS_NONE_H
#define _MODEL_REACTIONS_NONE_H

#include "global.h"
#include "SmartFieldVec.h"
#include "SmartFieldMat.h"
#include "SmartFieldOpMat.h"
#include "Model_Reactions_Base.h"

class Model_Reactions_None : public Model_Reactions_Base
{

  public:
    Model_Reactions_None( int ncomp, int nicomp );
    ~Model_Reactions_None();

    void set_params();

    // calculate reaction terms d/dt( phi )
    // (All inputs/outputs are in Fourier Space!)
    void calc_reactions(const SmartFieldVec & phi, SmartFieldVec & rxn); //full reactions
    void calc_reactions_lin_exp(const SmartFieldVec & phi, SmartFieldVec & rxn_lin); //explicit linear reactions
    void calc_reactions_lin_imp(const SmartFieldVec & phi, SmartFieldOpMat & Rxn); // implicit linear reactions

    // scale for the chemical potential
    // (e.g. used in Model H)
    inline RealType rxn_scale() {return 1.;};

  private:

};

#endif //_MODEL_REACTIONS_NONE_H
