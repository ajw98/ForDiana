/***************************************************************
*
* Class defining a Flory Huggins free energy model
* 
* DRT -- Wed, 18 Aug 2015
*
****************************************************************/

#ifndef _BCS_BASE_H
#define _BCS_BASE_H

#include "global.h"
#include "Grid.h"
#include "SmartFieldOp.h"
#include "SmartFieldVec.h"
#include "SmartFieldOpMat.h"

class BCs_Base
{

  public:
    BCs_Base( int nicomp, Grid * gridarg ) :
      _NIcomp(nicomp), _CurrGrid( gridarg ) {};

    virtual ~BCs_Base() {};

    virtual void solve_BCs( SmartFieldOpMat & T, 
                            SmartFieldVec & rhs,
                            SmartFieldVec & x,
                            RealType t ) = 0;

    inline int BCFlag() { return _BCFlag; };

  protected:

    UInt _NIcomp; // # of components
    Grid* _CurrGrid;
    int _BCFlag; // the kind of BC

};

#endif //_BCS_BASE_H
