/***************************************************************
*
* Const Phi boundary condition solver
* 
* DRT -- Sat, 15 Apr 2017
*
****************************************************************/

#ifndef _BCS_TIMEDEP_VEL_H
#define _BCS_TIMEDEP_VEL_H

#include "global.h"
#include "Grid.h"
#include "SmartFieldOp.h"
#include "SmartFieldOpMat.h"
#include "MatInterFieldMat.h"
#include "BCs_Base.h"
#include "FieldStack.h"

class BCs_TimeDep_Vel : public BCs_Base
{

  public:
    BCs_TimeDep_Vel( std::string filename, int nicomp, Grid * gridarg );
    ~BCs_TimeDep_Vel();

    void solve_BCs( SmartFieldOpMat & T,
                    SmartFieldVec & rhs, 
                    SmartFieldVec & x,
                    RealType t );

    void calc_BCs( MatInterFieldMat & T, VecInterFieldVec & b, RealType t );

    inline int BCFlag() {return _BCFlag;};

  private:

    Vec2dReal _Vel0, _VelL;
    Vec1dReal _BC_times;
    UInt _Nt;   // # of time points in BC file
    int _Dim;
};

#endif //_BCS_TIMEDEP_VEL_H
