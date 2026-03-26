/***************************************************************
*
* Class defining a Flory Huggins free energy model
* 
* DRT -- Wed, 18 Aug 2015
*
****************************************************************/

#ifndef _MODEL_MOBILITY_CONST_H
#define _MODEL_MOBILITY_CONST_H

#include "global.h"
#include "SmartFieldVec.h"
#include "SmartFieldMat.h"
#include "Model_Mobility_Base.h"
#include "Model_Visc_Base.h"

class Model_Mobility_Const : public Model_Mobility_Base
{

  public:
    Model_Mobility_Const( int ncomp, int nicomp, Model_Visc_Base* visc_model );
    ~Model_Mobility_Const();

    void set_params();
    virtual void set_mobility(const SmartFieldVec &phi, SmartFieldMat &M);

  private:

    // model parameters
    Vec2dReal _Mobility; // Constant mobility matrix

};

#endif //_MODEL_MOBILITY_CONST_H
