/***************************************************************
*
* Factory class for the boundary conditions
* 
* DRT -- Mon 16 May 2016
*
****************************************************************/

#ifndef _BCS_FACTORY_H
#define _BCS_FACTORY_H

#include "global.h"
#include "BCs_Base.h"
#include "BCs_Periodic.h"
#include "BCs_Const_Phi.h"
#include "BCs_Const_Vel.h"
#include "BCs_TimeDep_Phi.h"
#include "BCs_TimeDep_Vel.h"
#include "BCs_VIPS.h"

class BCs_Factory
{

  public:

    // --- constructor/destructor ---
    BCs_Factory();
    ~BCs_Factory();

    // --- BCs factory --- 
    BCs_Base* make_BCs( int key, int phivelkey, std::string filename,
                        int NIcomp, Grid* CurrGrid );

};


#endif //_BCS_FACTORY_H
 
