/***************************************************************
*
* Base class for the Operator Class
* 
* DRT -- 25 Jan 2017
*
****************************************************************/

#ifndef _OPERATORS_BASE_H
#define _OPERATORS_BASE_H

#include "global.h"
#include "SmartField.h"
#include "SmartFieldVec.h"
#include "SmartFieldMat.h"
#include "Grid.h"
#include "SmartFieldOp.h"
#include "SmartFieldOpMat.h"
#include "BCs_Base.h"

class Operators_Base
{

  public:

    // --- constructor/destructor ---

    Operators_Base( Grid * gridarg, BCs_Base * bcsarg );
    virtual ~Operators_Base();
  
    // for setting up BCs
    void setup_BCs( std::string filename, int BCNcomp );
    void cleanup_BCs();

    // --- explicit operators ---

    // *** All operators act on fields in k-space ***
    // ***     and return fields in k-space       ***

    // Recall: f is (_NIcomp x 1), graidents are (_Ndim x 1)
    virtual void Del_f_ex( const SmartField & f, 
                                 SmartFieldVec & Del_f ) = 0;
    virtual void Del_f_ex( const SmartFieldVec & f, 
                                 SmartFieldMat & Del_f ) = 0;
    virtual void Div_f_ex( const SmartFieldVec & f,
                                 SmartField & Div_f ) = 0;
    virtual void Div_f_ex( const SmartFieldMat & f,
                                 SmartFieldVec & Div_f ) = 0;
    virtual void Del2_f_ex( const SmartField & f, 
                                  SmartField & Del2_f ) = 0;
    virtual void Del2_f_ex( const SmartFieldVec & f, 
                                  SmartFieldVec & Del2_f ) = 0;
    virtual void Del2_f_ex_inplace( SmartField &f ) = 0;
    virtual void Coulomb_f_ex( const SmartField & f, 
                                  SmartField & Coul_f ) = 0;
    virtual void Del_Coulomb_f_ex( const SmartField & f, 
                                  SmartFieldVec & DelCoul_f ) = 0;
    virtual void Del3_f_ex( const SmartFieldVec & f, 
                                  SmartFieldMat & Del3_f ) = 0;
    virtual void Del4_f_ex( const SmartFieldVec & f, 
                                  SmartFieldVec & Del4_f ) = 0;
    virtual void A_Del_f_ex ( const SmartFieldVec &f,
                              const SmartFieldMat &A,
                               SmartFieldMat & Adf ) = 0;
    virtual void Del_A_Del_f_ex ( const SmartFieldVec & f, 
                                  const SmartFieldMat & A, 
                                        SmartFieldVec & dAdf ) = 0;
    virtual void Del_A_Del_f_ex ( const SmartFieldMat & Del_f, 
                              const SmartFieldMat & A, 
                                    SmartFieldVec & dAdf ) = 0;
    virtual void Del_A_Del3_f_ex( const SmartFieldVec & f, 
                                  const SmartFieldMat & A, 
                                        SmartFieldVec & d3Adf) = 0;
    virtual void transverse_f_ex ( const SmartFieldVec & f,
                                        SmartFieldVec & T_f ) = 0;

    // --- implicit operators (linear only) ---

    virtual void Eye_im( SmartFieldOpMat & T ) = 0;
    virtual void Del2_im( SmartFieldOp & T ) = 0;
    // virtual void Del2_im( SmartFieldOpMat & T ) = 0;
    virtual void Coulomb_im( SmartFieldOp & T ) = 0;
    virtual void A_Del2_im( const Vec2dFieldType &A, 
                                  SmartFieldOpMat & T ) = 0;
    virtual void A_Coulomb_im( const Vec2dFieldType &A, 
                                  SmartFieldOpMat & T ) = 0;
    virtual void A_Del4_im( const Vec2dFieldType &A, 
                                  SmartFieldOpMat & T ) = 0;
    virtual void A_Del2_F_im( const Vec2dFieldType &A, 
                              const SmartFieldOpMat &F, 
                                    SmartFieldOpMat & T ) = 0;

    virtual void Solve_im_matrix( SmartFieldOpMat & A, 
                                  SmartFieldVec & rhs, 
                                  SmartFieldVec & x,
                                  RealType t ) = 0;

    virtual void Solve_im_matrix( SmartFieldOp & A, 
                                  SmartField & rhs, 
                                  SmartField & x,
                                  RealType t ) = 0;

    // getters/setters
    inline Grid* GetGrid() { return _CurrGrid; };
    inline BCs_Base* GetBCs() { return _CurrBCs; };

    // Testing routine
    void TestOperators();

  protected:

    Grid *_CurrGrid; // Need this variable or the grid goes out of scope!
    BCs_Base *_CurrBCs; // Each operator needs has a set of BCs to fully 
                        // define operators and to execute solves

};

#endif //_OPERATORS_BASE_H
 
