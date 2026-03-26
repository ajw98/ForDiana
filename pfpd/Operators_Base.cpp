/***************************************************************
*
* Base class for the operators
* 
* DRT -- Mon, 20 Mar 2017
*
****************************************************************/

#include "Operators_Base.h"

Operators_Base :: Operators_Base( Grid * gridarg, BCs_Base * bcsarg ) :
_CurrGrid( gridarg ),
_CurrBCs( bcsarg )
{ // {{{
} // }}}

Operators_Base::~Operators_Base() 
{ // {{{
} // }}}

void Operators_Base::TestOperators()
{ // {{{

  // test the Operators class
//  {
//    SmartFieldVec phi(2);
//    SmartFieldVec del2_phi(2);
//    SmartFieldVec del4_phi(2);
//
//    SmartFieldMat A(2, 2);
//    A[0][0] = FieldType(1.);
//    A[0][1] = FieldType(2.);
//    A[1][0] = FieldType(3.);
//    A[1][1] = FieldType(4.);
//
//    phi[0].readfield("phi1.in");
//    phi[1].readfield("phi2.in");
//    
//    phi.fft();
//
//    this->Del2_f_ex( phi, del2_phi );
//    this->Del4_f_ex( phi, del4_phi );
//
//    phi.ifft();
//    del2_phi.ifft();
//    del4_phi.ifft();
//
//    phi[0].writerealfield("phi_test.dat");
//    del2_phi[0].writerealfield("del2phi.dat");
//    del4_phi[0].writerealfield("del4phi.dat");
//
//    // implicit operators
//    Vec2dFieldType B(2, Vec1dFieldType(2, 0.));
//    B[0][0] = FieldType(1.);
//    B[0][1] = FieldType(0.);
//    B[1][0] = FieldType(0.);
//    B[1][1] = FieldType(1.);
//
//    SmartFieldOpMat T1(2, 2);
//    SmartFieldOpMat Del2_Phi(2, 2);
//    SmartFieldVec del2_im(2);
//    SmartFieldVec del4_im(2);
//    SmartFieldVec del2_f_im(2);
//
//    phi.fft();
//
//    // --- del^2
//    T1.zero();
//    T1.setflag_inrealspace( false );
//    this->A_Del2_im(B, T1);
//    T1.OpDot( phi, del2_im );
//    del2_im.ifft();
//    del2_im[0].writerealfield("del2_im.dat");
//
//    // --- del^4
//    T1.zero();
//    T1.setflag_inrealspace( false );
//    this->A_Del4_im(B, T1);
//    T1.OpDot( phi, del4_im );
//    del4_im.ifft();
//    del4_im[0].writerealfield("del4_im.dat");
//
//    // --- del^2( del^2 phi )
//    Del2_Phi.zero();
//    Del2_Phi.setflag_inrealspace( false );
//    this->A_Del2_im(B, Del2_Phi );
//    
//    T1.zero();
//    T1.setflag_inrealspace( false );
//    this->A_Del2_F_im(B, Del2_Phi, T1);
//
//    T1.OpDot( phi, del2_f_im );
//    del2_f_im.ifft();
//    del2_f_im[0].writerealfield("del2_f_im.dat");
//
//  }
//
//  // testing implicit solver 
//  {
//    SmartFieldVec rhs(2);
//    SmartFieldVec lhs(2);
//    SmartFieldOpMat T1(2, 2);
//
//    rhs[0].readfield("phi1.in");
//    rhs[1].readfield("phi2.in");
//    rhs.fft();
//
//    Vec2dFieldType B(2, Vec1dFieldType(2, 0.));
//    B[0][0] = FieldType(1.);
//    B[0][1] = FieldType(0.);
//    B[1][0] = FieldType(0.);
//    B[1][1] = FieldType(1.);
//
//    // --- del^2
//    T1.zero();
//    T1.setflag_inrealspace( false );
//    this->Eye_im(T1);
//    this->A_Del2_im(B, T1);
//
//    lhs.zero();
//    lhs.setflag_inrealspace( false );
//    T1.OpDot( rhs, lhs );
//
//    lhs.zero();
//    lhs.setflag_inrealspace( false );
//    this->Solve_im_matrix(T1, rhs, lhs);
//
//    lhs.ifft();
//    lhs[0].writerealfield("lhs1.dat");
//    lhs[1].writerealfield("lhs2.dat");
//  }

  // test divergence operator
  {
    SmartFieldVec vel(2);

    vel[0].readfield("vx.in");
    vel[1].readfield("vy.in");
    vel.fft();

    SmartField div_vel;
    div_vel.setflag_inrealspace(false);
    this->Div_f_ex( vel, div_vel );
    div_vel.ifft();
    div_vel.writerealfield("div_vel.dat");

    SmartFieldMat grad_vel(2,2);
    grad_vel.setflag_inrealspace(false);
    this->Del_f_ex( vel, grad_vel );
    grad_vel.ifft();
    grad_vel(0,0).writerealfield("grad_vel_xx.dat");
    grad_vel(0,1).writerealfield("grad_vel_xy.dat");
    grad_vel(1,0).writerealfield("grad_vel_yx.dat");
    grad_vel(1,1).writerealfield("grad_vel_yy.dat");

    // get del2 from tensor
    SmartFieldVec div_del_vel(2);
    grad_vel.fft();
    this->Div_f_ex( grad_vel, div_del_vel );
    div_del_vel.ifft();
    div_del_vel[0].writerealfield("div_del_vel_x.dat");
    div_del_vel[1].writerealfield("div_del_vel_y.dat");

    // get del2 from implicit operator
    SmartFieldOp Del2_op;
    Del2_op.zero();
    Del2_op.setflag_inrealspace(false);
    this->Del2_im( Del2_op );

    SmartFieldVec del2_vel(2);
    del2_vel.zero();
    del2_vel.setflag_inrealspace(false);
    Del2_op.OpDot( vel[0], del2_vel[0] );
    Del2_op.OpDot( vel[1], del2_vel[1] );
    del2_vel.ifft();
    del2_vel[0].writerealfield("del2_im_vel_x.dat");
    del2_vel[1].writerealfield("del2_im_vel_y.dat");
  
    // get del2 from explicit operator
    SmartFieldVec del2_vel_ex(2);
    del2_vel_ex.zero();
    del2_vel_ex.setflag_inrealspace(false);
    Del2_f_ex( vel, del2_vel_ex );
    del2_vel_ex.ifft();
    del2_vel_ex[0].writerealfield("del2_ex_vel_x.dat");
    del2_vel_ex[1].writerealfield("del2_ex_vel_y.dat");

    // test transverse operator
    SmartFieldVec tr_vel(2);
    transverse_f_ex( vel, tr_vel );
    tr_vel.ifft();
    tr_vel[0].writerealfield("tr_vel_x.dat");
    tr_vel[1].writerealfield("tr_vel_y.dat");
  
  }

} // }}}

