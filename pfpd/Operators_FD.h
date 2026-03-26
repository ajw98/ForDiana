/***************************************************************
*
* Class for the finite difference operators
* 
* DRT -- Fri, 13 Mar 2015
*
****************************************************************/

#ifndef _OPERATORS_FD_BASE_H
#define _OPERATORS_FD_BASE_H

#include "global.h"
#include "FFTlayout.h"
#include "pllhandler.h"
#include "Field.h"
#include "SmartField.h"
#include "SmartFieldVec.h"
#include "SmartFieldMat.h"
#include "Operators_Base.h"
#include "SmartFieldOpMat.h"

class Operators_FD : public Operators_Base, private CellListener
{

  public:

    // --- constructor/destructor ---

    Operators_FD( Grid * gridarg, BCs_Base * bcsarg );
    ~Operators_FD();

    // --- explicit operators ---

    // *** All operators act on fields in k-space ***
    // ***      and return fields in k-space      ***

    // Recall: f is (_NIcomp x 1), graidents are (_Ndim x 1)
    void Del_f_ex( const SmartField &f, SmartFieldVec & Del_f );
    void Del_f_ex( const SmartFieldVec &f, SmartFieldMat & Del_f );
    void Div_f_ex( const SmartFieldVec &f, SmartField & Div_f );
    void Div_f_ex( const SmartFieldMat &f, SmartFieldVec & Div_f );
    void Del2_f_ex( const SmartField &f, SmartField & Del2_f );
    void Del2_f_ex( const SmartFieldVec &f, SmartFieldVec & Del2_f );
    void Del2_f_ex_inplace( SmartField &f );
    void Coulomb_f_ex( const SmartField &f, SmartField & Coul_f );
    void Del_Coulomb_f_ex( const SmartField &f, SmartFieldVec & DelCoul_f );
    void Del3_f_ex( const SmartFieldVec &f, SmartFieldMat & Del3_f );
    void Del4_f_ex( const SmartFieldVec &f, SmartFieldVec & Del4_f );
    void A_Del_f_ex ( const SmartFieldVec &f, 
                                    const SmartFieldMat &A, 
                                          SmartFieldMat & Adf );
    void Del_A_Del_f_ex ( const SmartFieldVec &f, 
                          const SmartFieldMat &A, 
                          SmartFieldVec & dAdf );
    void Del_A_Del_f_ex ( const SmartFieldMat &Del_f, 
                      const SmartFieldMat &A, 
                      SmartFieldVec & dAdf );   
    void Del_A_Del3_f_ex( const SmartFieldVec &f, 
                          const SmartFieldMat &A, 
                          SmartFieldVec & dAd3f );
    void transverse_f_ex ( const SmartFieldVec & f,
                                SmartFieldVec & T_f );

    // --- implicit operators (linear only) ---

    void Eye_im( SmartFieldOpMat & T );
    void Del2_im( SmartFieldOp & T );
    void Coulomb_im( SmartFieldOp & T );
    void A_Del2_im( const Vec2dFieldType &A, SmartFieldOpMat & T );
    void A_Coulomb_im( const Vec2dFieldType &A, SmartFieldOpMat & T );
    void A_Del4_im( const Vec2dFieldType &A, SmartFieldOpMat & T );
    void A_Del2_F_im( const Vec2dFieldType &A, 
                      const SmartFieldOpMat &F, 
                            SmartFieldOpMat & T );

    // --- solvers ---

    void Solve_im_matrix( SmartFieldOpMat & A, 
                          SmartFieldVec & rhs, 
                          SmartFieldVec & x,
                          RealType t);

    void Solve_im_matrix( SmartFieldOp & A, 
                          SmartField & rhs, 
                          SmartField & x,
                          RealType t);

    // --- Finite difference formulas ---

    // Finite differences evaluated at nodes (safe ends)
    void d1f( Vec2dUInt & stencil, Vec2dReal & coeffs );
    void d2f( Vec2dUInt & stencil, Vec2dReal & coeffs );
    void d3f( Vec2dUInt & stencil, Vec2dReal & coeffs );
    void d4f( Vec2dUInt & stencil, Vec2dReal & coeffs );

    // Finite differences evaluated at half-nodes
    void d1f_half( Vec2dReal & stencil, Vec2dReal & coeffs );
    void d3f_half( Vec2dReal & stencil, Vec2dReal & coeffs );

    // Finite difference schemes, periodic (note different structure)
    void d1f_per( Vec1dInt & stencil, Vec1dFieldType & coeffs );
    void d2f_per( Vec1dInt & stencil, Vec1dFieldType & coeffs );
    void d3f_per( Vec1dInt & stencil, Vec1dFieldType & coeffs );
    void d4f_per( Vec1dInt & stencil, Vec1dFieldType & coeffs );

    // useful stencil operations
    void FD_der( const Vec2dUInt & stencil, const Vec2dReal & coeffs, 
                 const SmartField & f, SmartField & df );
    Vec1dUInt add_stencils( Vec1dReal & stencil_a, Vec1dReal & stencil_b );
    Vec1dUInt add_stencils( Vec1dReal & stencil_a, RealType stencil_b );

  private:
    void cacheKvectors();
    inline virtual void CellUpdated(){cacheKvectors();};

  protected:

    // --- Member Variables --
    
    // Note _K2 and _Del2 are not FFTable since real field with C2C layout.
    // The SmartField class does not give me the flexibility to specify this.
    // For now, ignore, but keep this in mind if there are problems later.

    SmartField _K2; // k-modes squared
    SmartField _Del2; // Laplacian operator del . del
    SmartFieldVec _Del; //gradient operator
    SmartFieldMat _Tij; // transverse projection operator

};


#endif // _OPERATORS_FD_BASE_H
 
