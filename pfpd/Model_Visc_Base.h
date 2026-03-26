/***************************************************************
*
* Base class for the viscosity model
* 
* DRT -- Wed, 18 Aug 2015
*
****************************************************************/
#ifndef _MODEL_VISC_BASE_H
#define _MODEL_VISC_BASE_H

#include "global.h"
#include "Field.h"
#include "SmartFieldVec.h"

class Model_Visc_Base
{

  public:
    Model_Visc_Base( int ncomp, int nicomp ):_Ncomp(ncomp), _NIcomp(nicomp), _ViscModelFlag(-1) {};
    virtual ~Model_Visc_Base() {};

    virtual void set_params() = 0;
    virtual void set_eta(const SmartFieldVec &phi, SmartField &eta) = 0;
  
  protected:

    int _ViscModelFlag; // Model # Flag
    int _Ncomp; // number of components
    int _NIcomp; // number of *independent* components
};

#endif //_MODEL_VISC_BASE_H
