/***************************************************************
*
* Class for the psuedospectral operators
* 
* DRT -- Fri, 13 Mar 2015
*
****************************************************************/

#ifndef _OPERATORS_PS_H
#define _OPERATORS_PS_H

#include "global.h"
#include "FFTlayout.h"
#include "pllhandler.h"
#include "Field.h"
#include "SmartField.h"
#include "SmartFieldVec.h"
#include "SmartFieldMat.h"
#include "Operators_Base.h"
#include "SmartFieldOpMat.h"

class Operators_PS : public Operators_Base, private CellListener
{

  public:

    // --- constructor/destructor ---

    Operators_PS( Grid * gridarg, BCs_Base * bcsarg );
    ~Operators_PS();

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
                                SmartFieldVec & dAd3f);
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
    void A_Del4_F_im( const Vec2dFieldType &A, 
                      const SmartFieldOpMat &F, 
                            SmartFieldOpMat & T );

    void Solve_im_matrix( SmartFieldOpMat & A, 
                          SmartFieldVec & rhs, 
                          SmartFieldVec & x,
                          RealType t );

    void Solve_im_matrix( SmartFieldOp & A, 
                          SmartField & rhs, 
                          SmartField & x,
                          RealType t );

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
    SmartFieldVec _DelWithNyquist;
    SmartField _Coul; // Coulomb operator 1/k^2
    SmartFieldMat _Tij; // transverse projection operator

};


#endif // _OPERATORS_PS_H
 
