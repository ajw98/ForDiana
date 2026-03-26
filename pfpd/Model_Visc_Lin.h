/***************************************************************
*
* Class defining a linear viscosity model
* 
* DRT -- Wed, 18 Aug 2015
*
****************************************************************/

#ifndef _MODEL_VISC_LIN_H
#define _MODEL_VISC_LIN_H

#include "global.h"
#include "Field.h"
#include "SmartFieldVec.h"
#include "Model_Visc_Base.h"

class Model_Visc_Lin : public Model_Visc_Base
{

  public:
    Model_Visc_Lin( int ncomp, int nicomp );
    ~Model_Visc_Lin();

    void set_params();
    void set_eta(const SmartFieldVec &phi, SmartField &eta);

  private:

    RealType _EtaR; // viscosity of ref
    Vec1dReal _Eta; // viscosity of all Ncomp
  
};

#endif //_MODEL_VISC_LIN_H
