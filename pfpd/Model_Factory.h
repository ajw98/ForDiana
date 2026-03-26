/***************************************************************
*
* Base class and factory for the viscosity model
* 
* DRT -- Wed, 18 Aug 2015
* JUG -- Sun, 29 Jan 2017
* DRT -- Wed, 05 Apr 2017
* RA  -- Fri, 20 Dec 2019
*
****************************************************************/
#ifndef _MODEL_FACTORY_H
#define _MODEL_FACTORY_H

#include "global.h"
#include "Operators_Base.h"
#include "Model_Energy_Base.h"
#include "Model_Energy_Linear.h"
#include "Model_Energy_FHG.h"
#include "Model_Energy_DUN.h"
#include "Model_Energy_A_B_ABlinear.h"
#include "Model_Energy_FHG_K3.h"
#include "Model_Mobility_Base.h"
#include "Model_Mobility_Const.h"
#include "Model_Mobility_Rouse.h"
#include "Model_Mobility_ScaledRouse.h"
#include "Model_Visc_Base.h"
#include "Model_Visc_Lin.h"
#include "Model_Visc_Sig.h"
#include "Model_Visc_Exp.h"
#include "Model_Visc_VFTH.h"
#include "Model_Reactions_Base.h"
#include "Model_Reactions_Binary.h"
#include "Model_Reactions_None.h"

class Model_Factory
{

  public:
    Model_Factory();
    ~Model_Factory();

    Model_Energy_Base*   make_energy_model( int key, Operators_Base & OpObj, int Ncomp, int NIcomp );
    Model_Visc_Base*     make_visc_model( int key, int Ncomp, int NIcomp );
    Model_Mobility_Base* make_mobility_model( int key, int Ncomp, int NIcomp, Model_Visc_Base* visc_model );
    Model_Reactions_Base* make_reactions_model( int key, int Ncomp, int NIcomp ); 
 
};

#endif //_MODEL_FACTORY_H
