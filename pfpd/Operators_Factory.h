/***************************************************************
*
* Basee class and factory for the time integration scheme
* 
* DRT -- Thurs 18 Jan 2017
*
****************************************************************/

#ifndef _OPERATORS_FACTORY_H
#define _OPERATORS_FACTORY_H

#include "global.h"
#include "Grid.h"
#include "Operators_Base.h"
#include "Operators_PS.h"
#include "Operators_FD.h"

class Operators_Factory
{

  public:

    // --- constructor/destructor ---
    Operators_Factory();
    ~Operators_Factory();

    // --- Operators factory --- 
    Operators_Base* make_operators( Grid * gridptr, BCs_Base * bcsptr );
    void recycle_operators( Operators_Base* ThisOp );

};


#endif //_OPERATORS_FACTORY_H
 
