/***************************************************************
*
* Base class and factory for the mobility model
* 
* DRT -- Thurs, 20 Aug 2015
*
****************************************************************/
#ifndef _MODEL_MOBILITY_BASE_H
#define _MODEL_MOBILITY_BASE_H

#include "global.h"
#include "Field.h"
#include "SmartFieldVec.h"
#include "SmartFieldMat.h"
#include "Model_Visc_Base.h"

class Model_Mobility_Base
{

  public:
    Model_Mobility_Base( int Ncomp, int NIcomp, Model_Visc_Base* visc_model );
    virtual ~Model_Mobility_Base();

    virtual void set_params() = 0;
    virtual void set_mobility(const SmartFieldVec &phi, SmartFieldMat &M) = 0;

    Vec2dFieldType get_mobility_max( const SmartFieldMat & Mobility );
  
  protected:

    int _MobilityModelFlag; // Model # Flag
    int _Ncomp; // number of components
    int _NIcomp; // number of *independent* components
    Model_Visc_Base* _ViscModel; // Viscosity model object (some models need to be able to calculate the viscosity)

};

#endif //_MODEL_MOBILITY_BASE_H
