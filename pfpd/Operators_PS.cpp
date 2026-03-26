/***************************************************************
*
* Class for the time integration scheme
* 
* DRT -- Fri, 13 Mar 2015
*
****************************************************************/

#include "Operators_PS.h"

// ----------------- Constructor/Destructor -----------------

Operators_PS :: Operators_PS( Grid * gridarg, BCs_Base * bcsarg ):
Operators_Base(gridarg, bcsarg), 
_Del( gridarg->Dim() ),
_DelWithNyquist( gridarg->Dim() ),
_Tij( gridarg->Dim(), gridarg->Dim() )
{ // {{{

  std::cout << "\n * Initiating Pseudospectral (PS) operators" << std::endl;

  cacheKvectors();
  // Cause cached K vectors to update if the simulation cell every changes
  _CurrGrid->GetLayout()->getCell().RegisterListener(this,"Operators_PS");
}

void Operators_PS::cacheKvectors()
{
  // --- get k-modes for derivatives ---
  // (Ideally this will be part of the field class
  //  sometime down the road. It would be nice to
  //  have a field class with derivatives and with
  //  support for vector fields.)

  // Gather K^2 (if using MPI, this is just the local process' mesh) into non-Field storage
  UInt ML(_CurrGrid->GetLayout()->getNPWlocal());
  Vec1dFieldType localmesh(ML);
  for(UInt m=0; m<ML; m++)
    localmesh[m] = _CurrGrid->GetLayout()->getK2Map()[m];

  // Fields for second derivatives.
  // square of k-modes
  _K2.setflag_inrealspace( false );
  _K2 = localmesh;
  // Laplacian operator in k-space: del^2 = -k^2
  _Del2.setflag_inrealspace( false );
  _Del2 = _K2;
  _Del2 *= -1.;
  // Coulomb operator in k-space: coul = 1/k^2
  _Coul.setflag_inrealspace( false );
  _Coul = 1.;
  _Coul /= _K2;
  _Coul.zero_k0();

  // Fields for first derivatives.
  // kx, ky, kz
  // Note that for odd-order derivatives the Nyquist mode should be eliminated.
  // However, some operations (e.g., del.(A(r) grad f(r)) should not have the Nyquist mode eliminated,
  // and we therefore maintain two operators.
  _Del.setflag_inrealspace( false );
  _DelWithNyquist.setflag_inrealspace(false);
  std::vector<int> offset(_CurrGrid->Dim());
  for(UInt idir=0; idir<_CurrGrid->Dim(); ++idir)
  {
    for(UInt m=0; m<ML; ++m)
      localmesh[m] = std::complex<double>(0., _CurrGrid->GetLayout()->getKiMap(idir)[m]); // Del_idir = i.k_idir
    _DelWithNyquist[idir] = localmesh; // Copy to Field storage

    // Zero out Nyquist mode for operators involving odd derivatives
    for(UInt m=0; m<ML; ++m)
    {
      _CurrGrid->GetLayout()->MapFromFFTindx(offset, m, true, true, false);
      if(offset[idir] == _CurrGrid->GetLayout()->getGridDim(idir)/2)
        localmesh[m] = 0.;
    }
    _Del[idir] = localmesh; // Copy to Field storage
  }

  // Transverse projection operator: Tij = 1/k^2 * (I - kk/k^2)
  _Tij.setflag_inrealspace( false );
  SmartField K2inv(1.0);
  K2inv.setflag_inrealspace(false);
  K2inv /= _K2;
  for( UInt i=0; i<_CurrGrid->Dim(); ++i )
  {
    // Diagonals
    _Tij(i,i) = 2./3.;
    _Tij(i,i) *= K2inv;
    // Off-diagonals
    for( UInt j=i+1; j<_CurrGrid->Dim(); ++j )
    {
      _Tij(i,j) = _Del[i];
      _Tij(i,j) *= _Del[j];
      _Tij(i,j) *= K2inv;
      _Tij(i,j) += 1.0;
      _Tij(i,j) *= K2inv;

      _Tij(j,i) = _Tij(i,j); // Matrix symmetry
    }
  }

} // }}}

Operators_PS :: ~Operators_PS()
{ // {{{
  _CurrGrid->GetLayout()->getCell().DeregisterListener(this);
} // }}}

// ----------------- Explicit Operators -----------------

void Operators_PS::Del2_f_ex_inplace( SmartField &f )
{ // {{{

  if ( f.getflag_inrealspace() )
    codeerror_abort("Cannot perform gradient operations on a field in real space",__FILE__,__LINE__);

  f *= _Del2;

} // }}}

void Operators_PS::Del_f_ex( const SmartField &f, SmartFieldVec & Del_f )
{ // {{{

  // f: scalar
  // Del_f: Rows = _Dim

  if ( f.getflag_inrealspace() )
    codeerror_abort("Cannot perform gradient operations on a field in real space",__FILE__,__LINE__);

  Del_f.setflag_inrealspace(false);
  Del_f = _Del;
  Del_f *= f;

} // }}}

void Operators_PS::Del_f_ex( const SmartFieldVec &f, SmartFieldMat & Del_f )
{ // {{{

  // f: Rows = _NIcomp
  // Del_f: Rows = _NIcomp, _Cols = _Ndim

  if ( f.getflag_inrealspace() )
    codeerror_abort("Cannot perform gradient operations on a field in real space",__FILE__,__LINE__);

  Del_f.setflag_inrealspace(false);
  Del_f.outer(f, _Del);

} // }}}

void Operators_PS::Div_f_ex( const SmartFieldVec &f, SmartField & Div_f )
{ // {{{

  // f: Rows = _Dim
  // Div_f: (scalar)

  if ( f.getflag_inrealspace() )
    codeerror_abort("Cannot perform gradient operations on a field in real space",__FILE__,__LINE__);

  Div_f.setflag_inrealspace(false);
//  SmartFieldVec tmp(f.get_Nelem());
//  tmp.setflag_inrealspace(false);
//  tmp = _Del;
//  dot( tmp, f, Div_f );
  dot( _Del, f, Div_f );

} // }}}

void Operators_PS::Div_f_ex( const SmartFieldMat &f, SmartFieldVec & Div_f )
{ // {{{

  // f: Rows = _Dim, Columns = _Dim
  // Div_f: Rows = _Dim

  if ( f.getflag_inrealspace() )
    codeerror_abort("Cannot perform gradient operations on a field in real space",__FILE__,__LINE__);

  Div_f.setflag_inrealspace(false);
  Div_f.dot( _Del, f );

} // }}}

void Operators_PS::Del2_f_ex( const SmartField &f, SmartField & Del2_f )
{ // {{{

  if ( f.getflag_inrealspace() )
    codeerror_abort("Cannot perform gradient operations on a field in real space",__FILE__,__LINE__);

  Del2_f = f;
  Del2_f *= _Del2;

} // }}}

void Operators_PS::Del2_f_ex( const SmartFieldVec &f, SmartFieldVec & Del2_f )
{ // {{{

  // f: Rows = _NIcomp
  // Del2_f: Rows = _NIcomp

  if ( f.getflag_inrealspace() )
    codeerror_abort("Cannot perform gradient operations on a field in real space",__FILE__,__LINE__);

  Del2_f = f;
  Del2_f *= _Del2;

} // }}}

void Operators_PS::Coulomb_f_ex( const SmartField &f, SmartField & Coul_f )
{ // {{{

  // f: Rows = _NIcomp
  // Del2_f: Rows = _NIcomp

  if ( f.getflag_inrealspace() )
    codeerror_abort("Cannot perform gradient operations on a field in real space",__FILE__,__LINE__);

  Coul_f = f;
  Coul_f *= _Coul;

} // }}}

void Operators_PS::Del_Coulomb_f_ex( const SmartField &f, SmartFieldVec & DelCoul_f )
{ // {{{

  // f: Rows = _NIcomp
  // Del2_f: Rows = _NIcomp

  if ( f.getflag_inrealspace() )
    codeerror_abort("Cannot perform gradient operations on a field in real space",__FILE__,__LINE__);

  DelCoul_f.setflag_inrealspace(false);
  DelCoul_f = _Del;
  DelCoul_f *= f;
  DelCoul_f *= _Coul;

} // }}}

void Operators_PS::Del3_f_ex( const SmartFieldVec &f, SmartFieldMat & Del3_f )
{ // {{{

  // f: Rows = _NIcomp
  // Del_f: Rows = _NIcomp, _Cols = _Ndim

  if ( f.getflag_inrealspace() )
    codeerror_abort("Cannot perform gradient operations on a field in real space",__FILE__,__LINE__);

  Del3_f.setflag_inrealspace(false);
  Del3_f.outer(f, _Del);
  Del3_f *= _Del2;

} // }}}

void Operators_PS::Del4_f_ex( const SmartFieldVec &f, SmartFieldVec & Del4_f )
{ // {{{

  // f: Rows = _NIcomp
  // Del2_f: Rows = _NIcomp

  if ( f.getflag_inrealspace() )
    codeerror_abort("Cannot perform gradient operations on a field in real space",__FILE__,__LINE__);

  Del4_f = f;
  Del4_f *= _Del2;
  Del4_f *= _Del2;

} // }}}

void Operators_PS::A_Del_f_ex ( const SmartFieldVec &f, 
                                    const SmartFieldMat &A, 
                                          SmartFieldMat & Adf )
{ // {{{

  // f: Rows = _NIcomp
  // dAdf: Rows = _NIcomp

  SmartFieldMat tmp( f.get_Nelem(), _CurrGrid->Dim() ); // Nelem x Dim

  if ( f.getflag_inrealspace() )
  {
    std::cout << "***(Error) Cannot perform gradient operations on a field in real space.***\n";
    exit(1);
  }

  if ( ! A.getflag_inrealspace() )
  {
    std::cout << "***(Error) A matrix must be in real space***\n";
    exit(1);
  }

  tmp.setflag_inrealspace(false);
  tmp.outer(f, _Del); // tmp = NIcomp x Dim

  tmp.ifft(); // go to real space for multiplication
  Adf.setflag_inrealspace(true);
  Adf.dot(A, tmp); // (NIcomp x NIcomp) . (NIcomp x Dim)
  Adf.fft(); // back to k-space

} // }}}

void Operators_PS::Del_A_Del_f_ex ( const SmartFieldVec &f, 
                                    const SmartFieldMat &A, 
                                          SmartFieldVec & dAdf )
{ // {{{

  // f: Rows = _NIcomp
  // dAdf: Rows = _NIcomp

  SmartFieldMat tmp( dAdf.get_Nelem(), _CurrGrid->Dim() ); // Nelem x Dim
  SmartFieldMat tmp2( dAdf.get_Nelem(), _CurrGrid->Dim() ); // Nelem x Dim

  if ( f.getflag_inrealspace() )
    codeerror_abort("Cannot perform gradient operations on a field in real space",__FILE__,__LINE__);

  if ( ! A.getflag_inrealspace() )
  {
    std::cout << "***(Error) A matrix must be in real space***\n";
    exit(1);
  }

  tmp.setflag_inrealspace(false);
  tmp.outer(f, _DelWithNyquist); // tmp = NIcomp x Dim

  tmp.ifft(); // go to real space for multiplication
  tmp2.setflag_inrealspace(true);
  tmp2.dot(A, tmp); // (NIcomp x NIcomp) . (NIcomp x Dim)
  tmp2.fft(false); // back to k-space. The 1/M is omitted for EACH of the Dim*Nelem elements

  dAdf.setflag_inrealspace(false);
  dAdf.dot(tmp2, _DelWithNyquist); // (NIcomp x Dim) . Dim
  // Each element of the result needs the 1/M scaling
  dAdf *= RealType(1./_CurrGrid->GetLayout()->getNPWglobal());

} // }}}

void Operators_PS::Del_A_Del_f_ex ( const SmartFieldMat &Del_f, 
                                const SmartFieldMat &A, 
                                      SmartFieldVec & dAdf )
{ // {{{

  // f: Rows = _NIcomp
  // dAdf: Rows = _NIcomp

  //SmartFieldMat tmp( dAdf.get_Nelem(), _CurrGrid->Dim() ); // Nelem x Dim
  SmartFieldMat tmp2( dAdf.get_Nelem(), _CurrGrid->Dim() ); // Nelem x Dim

  if ( ! Del_f.getflag_inrealspace() )
    codeerror_abort("Del_A_Del_f_ex with Del_f input needs Del_f in real space",__FILE__,__LINE__);

  if ( ! A.getflag_inrealspace() )
  {
    std::cout << "***(Error) A matrix must be in real space***\n";
    exit(1);
  }

  //tmp.setflag_inrealspace(false);
  //tmp = Del_f; // tmp = NIcomp x Dim

  //Del_f.ifft(); // go to real space for multiplication
  tmp2.setflag_inrealspace(true);
  tmp2.dot(A, Del_f); // (NIcomp x NIcomp) . (NIcomp x Dim)
  tmp2.fft(false); // back to k-space. The 1/M is omitted for EACH of the Dim*Nelem elements

  dAdf.setflag_inrealspace(false);
  dAdf.dot(tmp2, _Del); // (NIcomp x Dim) . Dim
  // Each element of the result needs the 1/M scaling
  dAdf *= RealType(1./_CurrGrid->GetLayout()->getNPWglobal());

} // }}}

void Operators_PS::Del_A_Del3_f_ex( const SmartFieldVec &f, 
                                    const SmartFieldMat &A, 
                                          SmartFieldVec & dAd3f)
{ // {{{

  // f: Rows = _NIcomp
  // dAdf: Rows = _NIcomp

  SmartFieldMat tmp( dAd3f.get_Nelem(), _CurrGrid->Dim() ); // Nelem x Dim
  SmartFieldMat tmp2( dAd3f.get_Nelem(), _CurrGrid->Dim() ); // Nelem x Dim

  if ( f.getflag_inrealspace() )
    codeerror_abort("Cannot perform gradient operations on a field in real space",__FILE__,__LINE__);

  if ( ! A.getflag_inrealspace() )
  {
    std::cout << "***(Error) A matrix must be in real space***\n";
    exit(1);
  }

  tmp.setflag_inrealspace(false);
  tmp.outer(f, _DelWithNyquist); // tmp = NIcomp x Dim
  tmp *= _Del2;

  tmp.ifft(); // go to real space for multiplication
  tmp2.setflag_inrealspace(true);
  tmp2.dot(A, tmp); // (NIcomp x NIcomp) . (NIcomp x Dim)
  tmp2.fft(); // back to k-space

  dAd3f.setflag_inrealspace(false);
  dAd3f.dot(tmp2, _DelWithNyquist); // (NIcomp x Dim) . Dim

} // }}}

void Operators_PS::transverse_f_ex ( const SmartFieldVec &f, 
                                           SmartFieldVec & T_f )
{ // {{{

  //
  // transverse = del^2 - del del
  // ** this may help with computing the velocity
  //    without so much round off error in Model H **
  //
  // ( del^2 - del del ) f = 
  //   [  dyy + dzz  del2 - dxy  del2 - dxz ] [ fx ] 
  //   [ del2 - dyx   dxx + dzz  del2 - dyz ] [ fy ]
  //   [ del2 - dzx  del2 - dzy   dyy + dzz ] [ fz ]
  //
  //   [  (dyy + dzz) fx + (del2 - dxy) fy + (del2 - dxz) fz ]
  //   [ (del2 - dyx) fx +  (dxx + dzz) fy + (del2 - dyz) fz ]
  //   [ (del2 - dzx) fx + (del2 - dzy) fy +  (dyy + dzz) fz ]
  //

  if ( f.getflag_inrealspace() )
    codeerror_abort("Cannot perform gradient operations on a field in real space",__FILE__,__LINE__);

  SmartFieldMat T( f.get_Nelem(), f.get_Nelem() );
  //T.zero();
  //T.setflag_inrealspace( false );
  T.outer( _Del, _Del );

  //T_f.zero();
  //T_f.setflag_inrealspace( false );
  T_f.dot( T, f );

  //SmartFieldVec tmpvec( f.get_Nelem() );
  //tmpvec.zero();
  //tmpvec.setflag_inrealspace( false );
  //tmpvec = f;
  SmartFieldVec tmpvec( f );
  tmpvec *= _Del2;

  T_f.axpy_inplace(tmpvec,RealType(-1.));
//  T_f *= RealType(-1.);
//  T_f += tmpvec;

} // }}}

// ----------------- Implicit Operators -----------------

void Operators_PS::Eye_im( SmartFieldOpMat & T )
{ // {{{

  //SmartFieldMat Eye( T.Nrow(), T.Ncol() );
  //Eye.setflag_inrealspace(false);
  //Eye.eye( FieldType(1.) );
  //T.AddBand(0, Eye);

  Vec2dFieldType Eye( T.Nrow(), Vec1dFieldType( T.Ncol(), FieldType(0.)) );
  for(UInt i=0; i<T.Nrow(); ++i)
    if(i < T.Ncol())
      Eye[i][i] = 1.;
  T.AddBand(0, Eye);

} // }}}

void Operators_PS::Del2_im( SmartFieldOp & T )
{ // {{{

  if (T.getflag_inrealspace() != false)
  {
    std::cout << "*** Error in Operators_PS::Del2_im ***\n";
    std::cout << "T is not in Fourier space.\n";
    exit(1);
  }
   
  T.AddBand(0, _Del2);

} // }}}

void Operators_PS::Coulomb_im( SmartFieldOp & T )
{ // {{{

  if (T.getflag_inrealspace() != false)
  {
    std::cout << "*** Error in Operators_PS::Coulomb_im ***\n";
    std::cout << "T is not in Fourier space.\n";
    exit(1);
  }
   
  T.AddBand(0, _Coul);

} // }}}

void Operators_PS::A_Del2_im( const Vec2dFieldType &A, SmartFieldOpMat & T )
{ // {{{

  SmartFieldMat A_Del2( T.Nrow(), T.Ncol() );

  A_Del2.setflag_inrealspace(false);
  A_Del2 = A;
  A_Del2 *= _Del2;
  T.AddBand(0, A_Del2);

} // }}}

void Operators_PS::A_Coulomb_im( const Vec2dFieldType &A, SmartFieldOpMat & T )
{ // {{{

  SmartFieldMat A_Coul( T.Nrow(), T.Ncol() );

  A_Coul.setflag_inrealspace(false);
  A_Coul = A;
  A_Coul *= _Coul;
  T.AddBand(0, A_Coul);

} // }}}

void Operators_PS::A_Del4_im( const Vec2dFieldType &A, SmartFieldOpMat & T )
{ // {{{

  SmartFieldMat A_Del4( T.Nrow(), T.Ncol() );

  A_Del4.setflag_inrealspace(false);
  A_Del4 = A;
  A_Del4 *= _Del2;
  A_Del4 *= _Del2;
  T.AddBand(0, A_Del4);

} // }}}

void Operators_PS::A_Del2_F_im( const Vec2dFieldType &A,
                                const SmartFieldOpMat &F,
                                      SmartFieldOpMat & T )
{ // {{{

  SmartFieldMat tmp( T.Nrow(), T.Ncol() );
  SmartFieldMat A_Del2_F( T.Nrow(), T.Ncol() );

  tmp.setflag_inrealspace(false);
  F.GetBand(0, tmp);
  tmp *= _Del2;
  A_Del2_F.setflag_inrealspace(false);
  A_Del2_F.dot( A, tmp );
  T.AddBand(0, A_Del2_F);

} // }}}

void Operators_PS::A_Del4_F_im( const Vec2dFieldType &A, 
                                const SmartFieldOpMat &F, 
                                      SmartFieldOpMat & T )
{ // {{{

  SmartFieldMat tmp( T.Nrow(), T.Ncol() );
  SmartFieldMat A_Del4_F( T.Nrow(), T.Ncol() );

  tmp.setflag_inrealspace(false);
  F.GetBand(0, tmp);
  tmp *= _Del2;
  tmp *= _Del2;
  A_Del4_F.setflag_inrealspace(false);
  A_Del4_F.dot( A, tmp );
  T.AddBand(0, A_Del4_F);

} // }}}

// ----------------- Solvers -----------------

void Operators_PS::Solve_im_matrix( SmartFieldOpMat & A, 
                                    SmartFieldVec & rhs, 
                                    SmartFieldVec & x,
                                    RealType t )
{ // {{{

  // solve A.x = rhs
  // where A = a I
  // t is time (needed for time-dependent BCs)

  SmartFieldMat a( A.Nrow(), A.Ncol() );
  a.setflag_inrealspace(false);
  A.GetBand(0, a);
  x.linsolve(a, rhs);

} // }}}

void Operators_PS::Solve_im_matrix( SmartFieldOp & A, 
                                    SmartField & rhs, 
                                    SmartField & x,
                                    RealType t )
{ // {{{

  // solve A.x = rhs
  // t is time (needed for time-dependent BCs)

  SmartField a;
  a.setflag_inrealspace(false);
  A.GetBand(0, a);
  x = rhs; 
  x /= a;

} // }}}

