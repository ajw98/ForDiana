/***************************************************************
*
* Basee class and factory for the time integration scheme
* 
* DRT -- Mon 16 May 2016
*
****************************************************************/

#ifndef _TIMEINT_FACTORY_H
#define _TIMEINT_FACTORY_H

#include "global.h"
#include "FFTlayout.h"
#include "pllhandler.h"
#include "Field.h"
#include "SmartFieldVec.h"
#include "SmartFieldMat.h"
#include "SmartField.h"
#include "TimeInt_Base.h"
#include "TimeInt_ModelB.h"
#include "TimeInt_ModelH.h"
#include "TimeInt_ModelB_Explicit.h"
#include "TimeInt_ModelH_GridDef.h"
#include "TimeInt_ModelB_DelMu.h"
#include "TimeInt_ModelH_DelMu.h"
#include "Operators_Base.h"

class TimeInt_Factory
{

  public:

    // --- constructor/destructor ---
    TimeInt_Factory();
    ~TimeInt_Factory();

    // --- TimeInt factory --- 
    TimeInt_Base* make_timeint( int key, 
                                int Ncomp, 
                                int NIcomp, 
                                Grid* CurrGrid );

};


#endif //_TIMEINT_FACTORY_H
 
