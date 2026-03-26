/***************************************************************
*
* Const Vel boundary condition solver
* 
* DRT -- Mon, 17 Apr 2017
*
****************************************************************/

#ifndef _BCS_CONST_VEL_H
#define _BCS_CONST_VEL_H

#include "global.h"
#include "Grid.h"
#include "SmartFieldOp.h"
#include "SmartFieldOpMat.h"
#include "MatInterFieldMat.h"
#include "BCs_Base.h"
#include "FieldStack.h"

class BCs_Const_Vel : public BCs_Base
{

  public:
    BCs_Const_Vel( std::string filename, int nicomp, Grid * gridarg );
    ~BCs_Const_Vel();

    void solve_BCs( SmartFieldOpMat & T,
                    SmartFieldVec & rhs, 
                    SmartFieldVec & x,
                    RealType t );

    void calc_BCs( MatInterFieldMat & T, VecInterFieldVec & b );

    inline int BCFlag() {return _BCFlag;};

  private:

    Vec1dReal _Vel0, _VelL;
    int _Dim;

};

#endif //_BCS_CONST_VEL_H
