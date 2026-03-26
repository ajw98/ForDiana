/***************************************************************
*
* Const Phi boundary condition solver
* 
* DRT -- Sat, 15 Apr 2017
*
****************************************************************/

#ifndef _BCS_CONST_PHI_H
#define _BCS_CONST_PHI_H

#include "global.h"
#include "Grid.h"
#include "SmartFieldOp.h"
#include "SmartFieldOpMat.h"
#include "MatInterFieldMat.h"
#include "BCs_Base.h"
#include "FieldStack.h"

class BCs_Const_Phi : public BCs_Base
{

  public:
    BCs_Const_Phi( std::string filename, int nicomp, Grid * gridarg );
    ~BCs_Const_Phi();

    void solve_BCs( SmartFieldOpMat & T,
                    SmartFieldVec & rhs,
                    SmartFieldVec & x,
                    RealType t );

    void calc_BCs( MatInterFieldMat & T, VecInterFieldVec & b );

    inline int BCFlag() {return _BCFlag;};

  private:

    int _GhostNodeFlag;
    UInt _NBCs; // # of BCs (related to the order of 
                // the derivatives of the operator)
    Vec2dInt _Order_xeq0, _Order_xeqL; // order of derivative of BCs
    Vec2dReal _BC_xeq0, _BC_xeqL; // value of BCs

};

#endif //_BCS_CONST_PHI_H
