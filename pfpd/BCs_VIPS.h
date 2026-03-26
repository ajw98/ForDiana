/***************************************************************
*
* VIPS boundary condition solver
* 
* JUG -- Wed 28 Nov 2018
*
****************************************************************/

#ifndef _BCS_VIPS_H
#define _BCS_VIPS_H

#include "global.h"
#include "Grid.h"
#include "SmartFieldOp.h"
#include "SmartFieldOpMat.h"
#include "MatInterFieldMat.h"
#include "BCs_Base.h"
#include "FieldStack.h"

class BCs_VIPS : public BCs_Base
{

  public:
    BCs_VIPS( std::string filename, int nicomp, Grid * gridarg );
    ~BCs_VIPS();

    void solve_BCs( SmartFieldOpMat & T,
                    SmartFieldVec & rhs,
                    SmartFieldVec & x,
                    RealType t );

    void calc_BCs( MatInterFieldMat & T, VecInterFieldVec & b );

    inline int BCFlag() {return _BCFlag;};

  private:

    int _GhostNodeFlag;
    int _DdxAcc;
    UInt _NBCs; // # of BCs (related to the order of 
                // the derivatives of the operator)
    Vec2dInt _Order_xeq0, _Order_xeqL; // order of derivative of BCs
    Vec2dReal _BC_xeq0, _BC_xeqL; // value of BCs 
    void check_input_Lx();
    void convert_BCs2VIPS ( SmartFieldVec & phi, int Nx );
    void calc_g ( RealType & g, RealType & phiP, RealType & phiN );
 
    // thermodynamic model parameters
    RealType _Nr; // Molecular weight
    RealType _CReg, _CRegN, _Delta; // Regularization parameters

    Vec1dReal _N; // degree of polymerization
    Vec1dReal _alpha; // relative degree of polymerization to Nr
    Vec2dReal _Chi; // Flory interaction parameters
    Vec2dReal _ChiN; // Flory interaction parameters * Nr
    Vec2dReal _ChiN_ijM; // = ChiN_ij - ChiN_iM - ChiN_jM (used repeatedly)
    Vec1dReal _Kappa; // Gradient coefficients
    int _EnergyModelFlag; // Model # Flag
    int _Ncomp; // number of components
    Vec2dFieldType _Hmax; // Maximum Hessian
    RealType _Hpp, _Hnn, _Hpn;
    void calc_Hess( RealType & phiP0, RealType & phiN0 );
    void print_jp ( RealType phiP0, RealType phiN0, RealType dp, RealType d3p,
		    RealType dn, RealType d3n );


};

#endif //_BCS_VIPS_H
