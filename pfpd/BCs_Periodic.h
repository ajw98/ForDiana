/***************************************************************
*
* Periodic boundary condition solver
* 
* DRT -- Sat 15 Apr 2017
*
****************************************************************/

#ifndef _BCS_H
#define _BCS_H

#include "global.h"
#include "Grid.h"
#include "SmartFieldOp.h"
#include "SmartFieldOpMat.h"
#include "MatInterFieldMat.h"
#include "BCs_Base.h"

class BCs_Periodic : public BCs_Base
{

  public:
    BCs_Periodic( std::string filename, int nicomp, 
                  Grid * gridarg );
    ~BCs_Periodic();

    void solve_BCs( SmartFieldOpMat & T, 
                    SmartFieldVec & rhs, 
                    SmartFieldVec & x,
                    RealType t);

    void solve_SMW( VecInterFieldVec & x, VecInterFieldVec & y, 
                    MatInterFieldMat & Z, MatInterFieldMat & V );

  private:

};

#endif //_BCS_H
