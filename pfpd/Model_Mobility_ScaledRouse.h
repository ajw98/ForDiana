/***************************************************************
*
* Class defining a Mobility model scaled by an arbitrary viscosity
* DRT -- Wed, 3 Jul 2019
*
****************************************************************/

#ifndef _MODEL_MOBILITY_SCALEDROUSE_H
#define _MODEL_MOBILITY_SCALEDROUSE_H

#include "global.h"
#include "SmartFieldVec.h"
#include "SmartFieldMat.h"
#include "Model_Mobility_Base.h"
#include "Model_Visc_Base.h"

class Model_Mobility_ScaledRouse : public Model_Mobility_Base
{

  public:
    Model_Mobility_ScaledRouse( int ncomp, int nicomp, Model_Visc_Base* visc_model );
    ~Model_Mobility_ScaledRouse();

    void set_params();
    virtual void set_mobility(const SmartFieldVec &phi, SmartFieldMat &M);

  private:

    Vec2dReal _Const; // Constant mobility matrix prefactor
};

#endif //_MODEL_MOBILITY_SCALEDROUSE_H



                               
