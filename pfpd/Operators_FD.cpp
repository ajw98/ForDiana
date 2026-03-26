/***************************************************************
*
* Class for the time integration scheme
* 
* DRT -- Fri, 13 Mar 2015
*
****************************************************************/

#include "Operators_FD.h"

// ----------------- Constructor/Destructor -----------------

Operators_FD :: Operators_FD( Grid * gridarg, BCs_Base * bcsarg ):
Operators_Base(gridarg, bcsarg),
_Del( gridarg->PSDim() ),
_Tij( gridarg->PSDim(), gridarg->PSDim() )
{ // {{{

  std::cout << "\n * Initiating Hybrid FD/PS operators" << std::endl;
  //
  cacheKvectors();
  _CurrGrid->GetLayout()->getCell().RegisterListener(this,"Operators_FD");
}


void Operators_FD::cacheKvectors()
{
  // --- Initialize things we need for gradients ---

  // --- get k-modes for derivatives ---
  // (Ideally this will be part of the field class
  //  sometime down the road. It would be nice to 
  //  have a field class with derivatives and with
  //  support for vector fields.)
  
  Vec1dFieldType k2(_CurrGrid->NPW(), FieldType(0.));
  Vec2dFieldType kx(  _CurrGrid->PSDim(), 
                      Vec1dFieldType( _CurrGrid->NPW(), 
                                      FieldType(0.)    )); // PSDim x NPW
  Vec1dFieldType tij(_CurrGrid->NPW(), FieldType(0.));

#ifdef __MPI__

  // For the MPI code, information about the grid is dispersed
  // on all of the nodes. Aggregate this info here to initialize
  // the derivative fields.

  int buf( _CurrGrid->GetLayout()->getBufferSize() );
  std::vector<double> k2_local( _CurrGrid->GetLayout()->getK2Map() ); // buf x 1
  std::vector< std::vector<double> > kx_local( _CurrGrid->GetLayout()->getKVecMap() ); // PSDim x buf

  std::vector<double> k2_global(_CurrGrid->NPW(), 0.); // NPW x 1
  std::vector< std::vector<double> > kx_global( _CurrGrid->PSDim(), 
                                                std::vector<double>(_CurrGrid->NPW(), 0.)); // PSDim x NPW

  MPI_Allgather( &k2_local[0], buf, MPI_DOUBLE, &k2_global[0], buf, MPI_DOUBLE, MPI_COMM_WORLD );
  for( UInt i=0; i<_CurrGrid->PSDim(); i++ )
  {
    MPI_Allgather( &kx_local[i][0], buf, MPI_DOUBLE, &kx_global[i][0], buf, MPI_DOUBLE, MPI_COMM_WORLD );
  }

#endif

  // need to cast the doubles from the field 
  // class to fieldtype variables for assignment
  // k2
  for( UInt m=0; m<_CurrGrid->NPW(); ++m)
  {
#ifdef __MPI__
    k2[m] = k2_global[m];
#else
    k2[m] = _CurrGrid->GetLayout()->getK2Map()[m];
#endif
  }

  // ky, kz
  UInt div = 1;
  for( int i=_CurrGrid->PSDim()-1; i>=0; i-- )
  {
    for( UInt m=0; m<_CurrGrid->NPW(); ++m )
    {
      // set Nyquist mode to zero in the relevant dimension
      if( ((m/div) % _CurrGrid->GetLayout()->getGridDim(i)) == _CurrGrid->GetLayout()->getGridDim(i)/2 )
      {
        kx[i][m] = 0.;
      }
      else
      {
#ifdef __MPI__
        kx[i][m] = kx_global[i][m];
#else
        kx[i][m] = _CurrGrid->GetLayout()->getKVecMap()[i][m];
#endif
      }
    }
    div *= _CurrGrid->GetLayout()->getGridDim(i);  
  }

  // Now we are finally ready to initialize
  // the fields for derivatives, etc.

  // square of k-modes
  _K2.setflag_inrealspace( false );
  _K2 = k2;

  // Laplacian operator in k-space: del^2 = -k^2
  _Del2.setflag_inrealspace( false );
  _Del2 = k2;
  _Del2 *= -1.;

  // Gradient operator in k-space: del = 1i*ky
  _Del.setflag_inrealspace( false );
  for( UInt i=0; i<_CurrGrid->PSDim(); i++ )
  {
    _Del[i] = kx[i];
    _Del[i] *= std::complex<double>(0.0, 1.0);
  }

  // For now, neglect hydrodynamics
  // Transverse projection operator: Tij = 1/k^2 * (I - kk/k^2)
  _Tij.setflag_inrealspace( false );
  for( UInt i=0; i<_CurrGrid->PSDim(); ++i )
  {
    for( UInt j=0; j<_CurrGrid->PSDim(); ++j )
    {
      for( UInt m=0; m<_CurrGrid->NPW(); ++m )
      {
        if (m == 0)
        { 
          tij[m] = 0.; //zero the zero k-mode (let bulk velocity equal 0)
        }
        else if (i == j)
        {
          tij[m] = 1./k2[m] * ( 1. - kx[i][m]*kx[j][m]/k2[m] );
        } 
        else 
        {
          tij[m] = 1./k2[m] * ( -kx[i][m]*kx[j][m]/k2[m] );
        }
      }
      _Tij(i,j) = tij;
    }
  }

} // }}}

Operators_FD :: ~Operators_FD()
{ // {{{
} // }}}

// ----------------- Explicit Operators -----------------
void Operators_FD::Del2_f_ex_inplace( SmartField &f )
{ // {{{

  codeerror_abort("FILL DEL2_F_EX_INPLACE OPERATOR",__FILE__,__LINE__);
} // }}}

void Operators_FD::Del_f_ex( const SmartField &f, SmartFieldVec & Del_f )
{ // {{{

  // f: (scalar)
  // Del_f: Rows = _Dim

  if ( f.getflag_inrealspace() )
  {
    std::cout << "***(Error) Operators_FD::Del_f_ex***\n";
    std::cout << "Cannot perform gradient operations on a field in real space.\n";
    exit(1);
  }
  
  // For FD coeffs: (stencil width x Nx)
  Vec2dReal coeffs(0, Vec1dReal(0, 0.));
  Vec2dUInt stencil(0, Vec1dUInt(0, 0));
  d1f( stencil, coeffs );

  SmartField tmp;
  tmp.zero();
  tmp.setflag_inrealspace(false);

  // del f = [d/dx f, del_yz f]
  Del_f.zero();
  Del_f.setflag_inrealspace(false);

  // x-component (FD)
  for (int j = 0; j < stencil.size(); j++)
  {
    tmp.subscript( f, stencil[j] ); // tmp = f[i][stencil[j]];
    Del_f[0].xpby_inplace( tmp, coeffs[j] );
  }

  // y-/z-component (PS)
  for(int j=0; j<_CurrGrid->PSDim(); j++) // outer product over dimension
  {
    Del_f[j+1] = f;
    Del_f[j+1] *= _Del[j];
  }

} // }}}

void Operators_FD::Del_f_ex( const SmartFieldVec &f, SmartFieldMat & Del_f )
{ // {{{

  // f: Rows = _NIcomp
  // Del_f: Rows = _NIcomp, _Cols = _Ndim

  if ( f.getflag_inrealspace() )
  {
    std::cout << "***(Error) Operators_FD::Del_f_ex***\n";
    std::cout << "Cannot perform gradient operations on a field in real space.\n";
    exit(1);
  }
  
  Del_f.setflag_inrealspace(false);

  // For FD coeffs: (stencil width x Nx)
  Vec2dReal coeffs(0, Vec1dReal(0, 0.));
  Vec2dUInt stencil(0, Vec1dUInt(0, 0));
  d1f( stencil, coeffs );

  SmartField tmp;
  tmp.setflag_inrealspace(false);

  // loop over components and get derivatives
  for (int i=0; i<f.get_Nelem(); i++) // loop over # components
  {

    // del f = [d/dx f, del_yz f]
    for(int j=0; j<_CurrGrid->PSDim()+1; j++) // zero all components x and yz
    {
      Del_f(i,j) = FieldType(0.);
    }

    // x-component (FD)
    for (int j = 0; j < stencil.size(); j++)
    {
      tmp.subscript( f[i], stencil[j] ); // tmp = f(i,stencil[j]);
      Del_f(i,0).xpby_inplace( tmp, coeffs[j] );
    }

    // y-/z-component (PS)
    for(int j=0; j<_CurrGrid->PSDim(); j++) // outer product over dimension
    {
      Del_f(i,j+1) = f[i];
      Del_f(i,j+1) *= _Del[j];
    }

  }

} // }}}

void Operators_FD::Div_f_ex( const SmartFieldVec &f, SmartField & Div_f )
{ // {{{

  // f: Rows = _Dim
  // Del_f: (scalar)

  if ( f.getflag_inrealspace() )
  {
    std::cout << "***(Error) Operators_FD::Div_f_ex***\n";
    std::cout << "Cannot perform divergence on a field in real space.\n";
    exit(1);
  }
  
  // Get FD coeffs: (stencil width x Nx)
  Vec2dReal coeffs(0, Vec1dReal(0, 0.));
  Vec2dUInt stencil(0, Vec1dUInt(0, 0));
  d1f( stencil, coeffs );

  Div_f.setflag_inrealspace(false);
  Div_f.zero();

  // div f = d/dx f + del_yz.f
  SmartField tmp;
  tmp.setflag_inrealspace(false);
  tmp.zero();

  // x-component (FD)
  for (int j = 0; j < stencil.size(); j++)
  {
    tmp.subscript( f[0], stencil[j] ); // tmp = f[0][stencil[j]];
    Div_f.xpby_inplace( tmp, coeffs[j] );
  }

  // y-/z-component (PS)
  // inner product over remaining dimensions
  for(int j=0; j<_CurrGrid->PSDim(); j++)
  {
    tmp = f[j+1]; // f has dimesions 0=x, 1=y, and 2=z
    tmp *= _Del[j]; // _Del has dimesions 0=y, 1=z
    Div_f += tmp;
  }

} // }}}

void Operators_FD::Div_f_ex( const SmartFieldMat &f, SmartFieldVec & Div_f )
{ // {{{

  // f: Rows = _Dim, Cols = _Dim
  // Del_f: Rows = _Dim

  if ( f.getflag_inrealspace() )
  {
    std::cout << "***(Error) Operators_FD::Div_f_ex***\n";
    std::cout << "Cannot perform divergence on a field in real space.\n";
    exit(1);
  }
  
  // Get FD coeffs: (stencil width x Nx)
  Vec2dReal coeffs(0, Vec1dReal(0, 0.));
  Vec2dUInt stencil(0, Vec1dUInt(0, 0));
  d1f( stencil, coeffs );

  // div f = d/dx f + del_yz.f
  Div_f.setflag_inrealspace(false);
  Div_f.zero();

  SmartField tmp;
  tmp.setflag_inrealspace(false);
  tmp.zero();

  // loop over columns 
  // - div . f is a left multiply, produces a row vector
  for (int j=0; j<f.get_Ncol(); j++)
  {

    // x-component (FD)
    for (int i = 0; i < stencil.size(); i++)
    {
      // tmp = f(0,j)[stencil[j]];
      tmp.subscript( f(0,j), stencil[i] );
      Div_f[j].xpby_inplace( tmp, coeffs[i] );
    }

    // y-/z-component (PS)
    // inner product over remaining dimensions
    for(int i=0; i<_CurrGrid->PSDim(); i++)
    {
      tmp = f(i+1,j); // f has dimesions 0=x, 1=y, and 2=z
      tmp *= _Del[i]; // _Del has dimesions 0=y, 1=z
      Div_f[j] += tmp;
    }

  }

} // }}}

void Operators_FD::Del2_f_ex( const SmartField &f, SmartField & Del2_f )
{ // {{{

  // f: Rows = _NIcomp
  // Del2_f: Rows = _NIcomp

  if ( f.getflag_inrealspace() )
  {
    std::cout << "***(Error) Operators_FD::Del2_f_ex***\n";
    std::cout << "Cannot perform gradient operations on a field in real space.\n";
    exit(1);
  }

  // del^2 f = d^2/dx^2 f + del_yz^2 f
  Del2_f = FieldType(0.);
  Del2_f.setflag_inrealspace(false);
  
  // d^2/dx^2 f
  Vec2dReal coeffs(0, Vec1dReal(0, 0.));
  Vec2dUInt stencil(0, Vec1dUInt(0, 0));
  d2f( stencil, coeffs );

  SmartField tmpf;
  tmpf.setflag_inrealspace(false);
  for (int j = 0; j < stencil.size(); j++)
  {
    tmpf.subscript( f, stencil[j] ); // tmp = f[stencil[j]];
    Del2_f.xpby_inplace( tmpf, coeffs[j] );
  }
    
  // del_yz^2 f
  SmartField tmp;
  tmp = f;
  tmp *= _Del2;
  Del2_f += tmp;

} // }}}

void Operators_FD::Del2_f_ex( const SmartFieldVec &f, SmartFieldVec & Del2_f )
{ // {{{

  // f: Rows = _NIcomp
  // Del2_f: Rows = _NIcomp

  if ( f.getflag_inrealspace() )
  {
    std::cout << "***(Error) Operators_FD::Del2_f_ex***\n";
    std::cout << "Cannot perform gradient operations on a field in real space.\n";
    exit(1);
  }

  // del^2 f = d^2/dx^2 f + del_yz^2 f
  Del2_f = FieldType(0.);
  Del2_f.setflag_inrealspace(false);
  
  // d^2/dx^2 f
  Vec2dReal coeffs(0, Vec1dReal(0, 0.));
  Vec2dUInt stencil(0, Vec1dUInt(0, 0));
  d2f( stencil, coeffs );

  SmartFieldVec tmpvec( f.get_Nelem() );
  tmpvec.setflag_inrealspace(false);
  for (int j = 0; j < stencil.size(); j++)
  {
    tmpvec.subscript( f, stencil[j] ); // tmp = f[stencil[j]];
    Del2_f.xpby_inplace( tmpvec, coeffs[j] );
  }
    
  // del_yz^2 f
  SmartFieldVec tmp( f.get_Nelem() );
  tmp = f;
  tmp *= _Del2;
  Del2_f += tmp;

} // }}}

void Operators_FD::Coulomb_f_ex( const SmartField &f, SmartField & Coul_f )
{ // {{{

  // f: Rows = _NIcomp
  // Coul_f: Rows = _NIcomp

  std::cout << "*** Caveat emptor: Operators_FD::Coulomb_f_ex() was called ***\n";
  std::cout << "***     It is not tested and the output may be garbage     ***\n";

  if ( f.getflag_inrealspace() )
  {
    std::cout << "***(Error) Operators_FD::Coulomb_f_ex***\n";
    std::cout << "Cannot perform gradient operations on a field in real space.\n";
    exit(1);
  }

  // del^2 f = d^2/dx^2 f + del_yz^2 f
  Coul_f = FieldType(0.);
  Coul_f.setflag_inrealspace(false);

  SmartField rhs; // Need this because f is const
  rhs = f;
  
  // d^2/dx^2 operator
  SmartFieldOp Del2;
  Del2_im( Del2 );

  // Solve: Del2 Coul_f = -f
  Solve_im_matrix( Del2, rhs, Coul_f, 0.);
  Coul_f *= -1.;

} // }}}

void Operators_FD::Del_Coulomb_f_ex( const SmartField &f, SmartFieldVec & DelCoul_f )
{ // {{{

  // f: Rows = _NIcomp
  // Coul_f: Rows = _NIcomp

  std::cout << "*** Caveat emptor: Operators_FD::Del_Coulomb_f_ex() was called ***\n";
  std::cout << "***       It is not tested and the output may be garbage       ***\n";

  if ( f.getflag_inrealspace() )
  {
    std::cout << "***(Error) Operators_FD::Del_Coulomb_f_ex***\n";
    std::cout << "Cannot perform gradient operations on a field in real space.\n";
    exit(1);
  }

  DelCoul_f.setflag_inrealspace(false);

  // Find Coulomb (f) first
  SmartField Coul_f;
  Coul_f.setflag_inrealspace(false);

  // Get gradient of Coul_f
  Del_f_ex( Coul_f, DelCoul_f);

} // }}}

void Operators_FD::Del3_f_ex( const SmartFieldVec &f, SmartFieldMat & Del3_f )
{ // {{{

  // f: Rows = _NIcomp
  // Del_f: Rows = _NIcomp, _Cols = _Ndim
  
  if ( f.getflag_inrealspace() )
  {
    std::cout << "***(Error) Operators_FD::Del3_f_ex***\n";
    std::cout << "Cannot perform gradient operations on a field in real space.\n";
    exit(1);
  }

  Del3_f.setflag_inrealspace(false);

  // For FD coeffs: (stencil width x Nx)
  Vec2dReal d1_coeffs (0, Vec1dReal(0, 0.));
  Vec2dUInt d1_stencil(0, Vec1dUInt(0, 0 ));
  Vec2dReal d2_coeffs (0, Vec1dReal(0, 0.));
  Vec2dUInt d2_stencil(0, Vec1dUInt(0, 0 ));
  Vec2dReal d3_coeffs (0, Vec1dReal(0, 0.));
  Vec2dUInt d3_stencil(0, Vec1dUInt(0, 0 ));

  d1f( d1_stencil, d1_coeffs );
  d2f( d2_stencil, d2_coeffs );
  d3f( d3_stencil, d3_coeffs );

  SmartField tmp1, tmp2;

  // loop over components and get derivatives
  for (int i=0; i<f.get_Nelem(); i++) // loop over # components
  {

    // del del^2 f = [d^3/dx^3 f + d/dx del_yz^2 f, 
    //                d^2/dx^2 del_yz f + del_yz del_yz^2 f]

    // --- x-component ---
    Del3_f(i,0) = FieldType(0.);

    // d^3/dx^3 f
    tmp1.setflag_inrealspace(false);
    tmp1 = 0.;
    for (int j = 0; j < d3_stencil.size(); j++)
    {
      tmp1.subscript( f[i], d3_stencil[j] ); // f[i][stencil[j]];
      Del3_f(i,0).xpby_inplace( tmp1, d3_coeffs[j] );
    }

    // d/dx del_yz^2
    
    tmp1.setflag_inrealspace(false);
    tmp2.setflag_inrealspace(false);
    tmp1 = 0.;
    tmp2 = 0.;
    for (int j = 0; j < d1_stencil.size(); j++)
    {
      tmp1.subscript( f[i], d1_stencil[j] ); // f[i][stencil[j]];
      tmp2.xpby_inplace( tmp1, d1_coeffs[j] );
    }
    tmp2 *= _Del2;

    Del3_f(i,0) += tmp2;

    // --- y-/z-component ---

    // del_yz d^2/dx^2 f
    tmp1.setflag_inrealspace(false);
    tmp2.setflag_inrealspace(false);
    tmp1 = 0.;
    tmp2 = 0.;
    for (int j = 0; j < d2_stencil.size(); j++)
    {
      tmp1.subscript( f[i], d2_stencil[j] ); // f[i][stencil[j]];
      tmp2.xpby_inplace( tmp1, d2_coeffs[j] );
    }

    for(int j=0; j<_CurrGrid->PSDim(); j++) // outer product over dimension
    {
      Del3_f(i,j+1) = tmp2;
      Del3_f(i,j+1) *= _Del[j];
    }

    // del_yz del_yz^2 f]
    for(int j=0; j<_CurrGrid->PSDim(); j++) // outer product over dimension
    {
      tmp1  = f[i];
      tmp1 *= _Del[j];
      tmp1 *= _Del2;

      Del3_f(i,j+1) += tmp1;
    }

  }

} // }}}

void Operators_FD::Del4_f_ex( const SmartFieldVec &f, SmartFieldVec & Del4_f )
{ // {{{

  // f: Rows = _NIcomp
  // Del2_f: Rows = _NIcomp

  if ( f.getflag_inrealspace() )
  {
    std::cout << "***(Error) Operators_FD::Del4_f_ex***\n";
    std::cout << "Cannot perform gradient operations on a field in real space.\n";
    exit(1);
  }

  // del^2 del^2 f =  d^4/dx^4 f + 2*d^2/dx^2 del_yz^2 f + del_yz^4 f]
  Del4_f.setflag_inrealspace(false);
  Del4_f = FieldType(0.);

  SmartFieldVec tmp1(f.get_Nelem());
  SmartFieldVec tmp2(f.get_Nelem());

  Vec2dReal d4coeffs (0, Vec1dReal(0, 0.));
  Vec2dUInt d4stencil(0, Vec1dUInt(0, 0 ));

  // d^4/dx^4 f
  tmp1.setflag_inrealspace(false);
  tmp1 = FieldType(0.);

  d4f( d4stencil, d4coeffs );

  for (int j = 0; j < d4stencil.size(); j++)
  {
    tmp1.subscript( f, d4stencil[j] ); // f[stencil[j]];
    Del4_f.xpby_inplace( tmp1, d4coeffs[j] );
  }

  // 2*d^2/dx^2 del_yz^2 f
  tmp1.setflag_inrealspace(false);
  tmp2.setflag_inrealspace(false);
  tmp1 = FieldType(0.);
  tmp2 = FieldType(0.);

  Vec2dReal d2coeffs (0, Vec1dReal(0, 0.));
  Vec2dUInt d2stencil(0, Vec1dUInt(0, 0 ));

  d2f( d2stencil, d2coeffs );

  for (int j = 0; j < d2stencil.size(); j++)
  {
    tmp1.subscript( f, d2stencil[j] ); // f[i][stencil[j]];
    tmp2.xpby_inplace( tmp1, d2coeffs[j] );
  }
  tmp2 *= _Del2;
  tmp2 *= FieldType(2.);

  Del4_f += tmp2;

  // del_yz^4 f
  tmp1 = f;
  tmp1 *= _Del2;
  tmp1 *= _Del2;
  Del4_f += tmp1;

} // }}}

void Operators_FD::A_Del_f_ex ( const SmartFieldVec &f, 
                                    const SmartFieldMat &A, 
                                          SmartFieldMat & Adf )
{ // {{{

  // f: Rows = _NIcomp
  // dAdf: Rows = _NIcomp
  // This is just a place holder for my intermediate testing.
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

void Operators_FD::Del_A_Del_f_ex ( const SmartFieldVec &f, 
                                    const SmartFieldMat &A, 
                                          SmartFieldVec & dAdf )
{ // {{{

  // f: Rows = _NIcomp
  // A: Rows = _NIcomp, Cols = _NIcomp
  // dAdf: Rows = _NIcomp

  if ( f.getflag_inrealspace() )
  {
    std::cout << "***(Error) Operators_FD::Del_A_Del_f_ex***\n";
    std::cout << "Cannot perform gradient operations on a field in real space.\n";
    exit(1);
  }
 
  if ( ! A.getflag_inrealspace() )
  {
    std::cout << "***(Error) Operators_FD::Del_A_Del_f_ex***\n";
    std::cout << "A matrix must be in real space\n";
    exit(1);
  }

  UInt NIcomp = f.get_Nelem();

  // del.(A del f) = d/dx( A df/dx) + del_yz . (A del_yz f)

  dAdf.setflag_inrealspace(false);
  dAdf = FieldType(0.);

  // d/dx( A df/dx) -- FD
  { // {{{

    SmartFieldMat A_half( NIcomp, NIcomp ); 
    SmartFieldVec df( NIcomp );
    SmartFieldVec tmpvec( NIcomp );
    SmartFieldMat tmpmat( NIcomp, NIcomp );
    tmpmat.zero();

    // For FD coeffs: (stencil width x Nx)
    Vec2dReal coeffs_outer(0, Vec1dReal(0, 0.));
    Vec2dReal coeffs_inner(0, Vec1dReal(0, 0.));
    Vec2dReal stencil_outer(0, Vec1dReal(0, 0));
    Vec2dReal stencil_inner(0, Vec1dReal(0, 0));

    d1f_half( stencil_outer, coeffs_outer );
    d1f_half( stencil_inner, coeffs_inner );

    if (stencil_inner[0].size() <= 0)
    {
      std::cout << "***(Error) in Operators_FD::Del_A_Del3_f_ex***\n";
      std::cout << "Stencil(s) allocated improperly.\n";
      exit(1);
    }

    Vec1dUInt iphalf(stencil_inner[0].size(), 0);
    Vec1dUInt imhalf(stencil_inner[0].size(), 0);

    for (int i=0; i<stencil_outer.size(); i++)
    {
      // get df
      tmpvec.setflag_inrealspace(false);
      tmpvec = FieldType(0.);

      df.setflag_inrealspace(false);
      df = FieldType(0.);
      for (int j=0; j<stencil_inner.size(); j++)
      {
        iphalf = add_stencils(stencil_outer[i], stencil_inner[j]);
        tmpvec.subscript( f, iphalf );
        df.xpby_inplace( tmpvec, coeffs_inner[j] );
      }

      // get A_half
      A_half.setflag_inrealspace(false);
      A_half.zero();
      iphalf = add_stencils( stencil_outer[i], 0.5 );
      imhalf = add_stencils( stencil_outer[i], -0.5 );
      A_half.subscript( A, iphalf );
      tmpmat.subscript( A, imhalf );
      A_half += tmpmat;
      A_half *= FieldType(0.5);

      // take dot product in real space
      df.ifft();
      tmpvec.setflag_inrealspace(true);
      tmpvec = FieldType(0.);
      tmpvec.dot(A_half, df);
      tmpvec.fft();

      // use coeffs to get outer derivative
      dAdf.xpby_inplace( tmpvec, coeffs_outer[i] );

    }
  } // }}}

  // del_yz .( A del_yz f) -- PS
  { // {{{

    SmartFieldMat Del_f(NIcomp, _CurrGrid->PSDim() );
    Del_f.setflag_inrealspace(false);
    Del_f.outer(f, _Del); // outer product over Ndim

    SmartFieldMat A_del_f( NIcomp, _CurrGrid->PSDim() );
    A_del_f.setflag_inrealspace(true);

    Del_f.ifft();
    A_del_f.dot(A, Del_f); // inner product over Ndim
    A_del_f.fft();

    SmartFieldVec del_A_del_f(NIcomp);
    del_A_del_f.setflag_inrealspace(false);
    del_A_del_f.dot(A_del_f, _Del);

    dAdf += del_A_del_f;

  } // }}}

} // }}}

void Operators_FD::Del_A_Del_f_ex ( const SmartFieldMat &Del_f, 
                                const SmartFieldMat &A, 
                                      SmartFieldVec & dAdf )
{ // {{{

  // f: Rows = _NIcomp
  // A: Rows = _NIcomp, Cols = _NIcomp
  // dAdf: Rows = _NIcomp

  std::cout << "*** Caveat emptor: Operators_FD::Del_A_Del_f_ex() was called ***\n";
  std::cout << "***    It is not tested and the output may be garbage    ***\n";

  if ( Del_f.getflag_inrealspace() )
  {
    std::cout << "***(Error) Operators_FD::Del_A_F_ex***\n";
    std::cout << "Cannot perform gradient operations on a field in real space.\n";
    exit(1);
  }
 
  if ( ! A.getflag_inrealspace() )
  {
    std::cout << "***(Error) Operators_FD::Del_A_F_ex***\n";
    std::cout << "A matrix must be in real space\n";
    exit(1);
  }

  UInt NIcomp = dAdf.get_Nelem();

  // del.(A del f) = d/dx( A df/dx) + del_yz . (A del_yz f)

  dAdf.setflag_inrealspace(false);
  dAdf = FieldType(0.);

  // d/dx( A df/dx) -- FD
  { // {{{

    SmartFieldMat A_half( NIcomp, NIcomp ); 
    SmartFieldVec df( NIcomp );
    SmartFieldVec tmpvec( NIcomp );
    SmartFieldMat tmpmat( NIcomp, NIcomp );

    // For FD coeffs: (stencil width x Nx)
    Vec2dReal coeffs_outer(0, Vec1dReal(0, 0.));
    Vec2dReal coeffs_inner(0, Vec1dReal(0, 0.));
    Vec2dReal stencil_outer(0, Vec1dReal(0, 0));
    Vec2dReal stencil_inner(0, Vec1dReal(0, 0));

    d1f_half( stencil_outer, coeffs_outer );
    d1f_half( stencil_inner, coeffs_inner );

    if (stencil_inner[0].size() <= 0)
    {
      std::cout << "***(Error) in Operators_FD::Del_A_F_ex***\n";
      std::cout << "Stencil(s) allocated improperly.\n";
      exit(1);
    }

    Vec1dUInt iphalf(stencil_inner[0].size(), 0);
    Vec1dUInt imhalf(stencil_inner[0].size(), 0);

    for (int i=0; i<stencil_outer.size(); i++)
    {
      // get df
      tmpvec.setflag_inrealspace(false);
      for (int m=0; m<NIcomp; m++)
        tmpvec[m] = Del_f(m,0);

      df.setflag_inrealspace(false);
      for (int j=0; j<stencil_inner.size(); j++)
      {
        iphalf = add_stencils(stencil_outer[i], stencil_inner[j]);
        df.subscript( tmpvec, iphalf );
      }

      // get A_half
      A_half.setflag_inrealspace(false);
      iphalf = add_stencils( stencil_outer[i], 0.5 );
      imhalf = add_stencils( stencil_outer[i], -0.5 );
      A_half.subscript( A, iphalf );
      tmpmat.subscript( A, imhalf );
      A_half += tmpmat;
      A_half *= FieldType(0.5);

      // take dot product in real space
      df.ifft();
      tmpvec.setflag_inrealspace(true);
      tmpvec = FieldType(0.);
      tmpvec.dot(A_half, df);
      tmpvec.fft();

      // use coeffs to get outer derivative
      dAdf.xpby_inplace( tmpvec, coeffs_outer[i] );

    }
  } // }}}

  // del_yz .( A del_yz f) -- PS
  { // {{{

    SmartFieldMat del_f(NIcomp, _CurrGrid->PSDim() );
    del_f.setflag_inrealspace(false);

    for (int i=0; i<NIcomp; i++)
    {
      for (int j=0; j<_CurrGrid->PSDim(); j++)
      {
        del_f(i,j) = Del_f(i,j+1);
      }
    }

    SmartFieldMat A_del_f( NIcomp, _CurrGrid->PSDim() );
    A_del_f.setflag_inrealspace(true);

    del_f.ifft();
    A_del_f.dot(A, del_f); // inner product over Ndim
    A_del_f.fft();

    SmartFieldVec del_A_del_f(NIcomp);
    del_A_del_f.setflag_inrealspace(false);
    del_A_del_f.dot(A_del_f, _Del);

    dAdf += del_A_del_f;

  } // }}}

}
void Operators_FD::Del_A_Del3_f_ex( const SmartFieldVec &f, 
                                    const SmartFieldMat &A, 
                                          SmartFieldVec & dAd3f)
{ // {{{

  // f: Rows = _NIcomp
  // A: Rows = _NIcomp, Cols = _NIcomp
  // dAd3f: Rows = _NIcomp

  if ( f.getflag_inrealspace() )
  {
    std::cout << "***(Error) Operators_FD::Del_A_Del3_f_ex***\n";
    std::cout << "Cannot perform gradient operations on a field in real space.\n";
    exit(1);
  }
 
  if ( ! A.getflag_inrealspace() )
  {
    std::cout << "***(Error) Operators_FD::Del_A_Del3_f_ex***\n";
    std::cout << "A matrix must be in real space\n";
    exit(1);
  }

  UInt NIcomp = f.get_Nelem();

  // del.(A del del2 f) = 
  //    d/dx ( A d^3 f/dx^3 ) + 
  //    d/dx ( A d/dx del_yz^2 f) +
  //    del_yz . (A del_yz d^2 f /dx^2) +
  //    del_yz . (A del_yz del_yz^2 f)

  dAd3f.setflag_inrealspace(false);
  dAd3f = FieldType(0.);

  // d/dx ( A d^3 f/dx^3) -- Pure FD
  { // {{{

    SmartFieldMat A_half( NIcomp, NIcomp );
    SmartFieldVec d3f( NIcomp );
    SmartFieldMat tmpmat( NIcomp, NIcomp );
    SmartFieldVec tmpvec( NIcomp );

    // For FD coeffs: (stencil width x Nx)
    Vec2dReal coeffs_outer(0, Vec1dReal(0, 0.));
    Vec2dReal coeffs_inner(0, Vec1dReal(0, 0.));
    Vec2dReal stencil_outer(0, Vec1dReal(0, 0));
    Vec2dReal stencil_inner(0, Vec1dReal(0, 0));

    d1f_half( stencil_outer, coeffs_outer );
    d3f_half( stencil_inner, coeffs_inner );

    if (stencil_inner[0].size() <= 0)
    {
      std::cout << "***(Error) in Operators_FD::Del_A_Del3_f_ex***\n";
      std::cout << "Stencil(s) allocated improperly.\n";
      exit(1);
    }

    Vec1dUInt iphalf(stencil_inner[0].size(), 0);
    Vec1dUInt imhalf(stencil_inner[0].size(), 0);

    for (int i=0; i<stencil_outer.size(); i++)
    {
      // get d3f (inner FD derivative)
      tmpvec.setflag_inrealspace(false);
      tmpvec = FieldType(0.);

      d3f.setflag_inrealspace(false);
      d3f = FieldType(0.);

      for (int j=0; j<stencil_inner.size(); j++)
      {
        iphalf = add_stencils(stencil_outer[i], stencil_inner[j]);  
        tmpvec.subscript( f, iphalf );
        d3f.xpby_inplace( tmpvec, coeffs_inner[j] );
      }

      // get A_half
      A_half.setflag_inrealspace(false);
      iphalf = add_stencils( stencil_outer[i], 0.5 );
      A_half.subscript( A, iphalf ); 
      imhalf = add_stencils( stencil_outer[i], -0.5 );
      tmpmat.subscript( A, imhalf ); 
      A_half += tmpmat;
      A_half *= FieldType(0.5);

      // take dot product in real space
      d3f.ifft();
      tmpvec.setflag_inrealspace(true);
      tmpvec.dot(A_half, d3f);
      tmpvec.fft();

      // get outer derivative
      dAd3f.xpby_inplace( tmpvec, coeffs_outer[i] );

    }

  } // }}}

  // d/dx ( A d/dx del_yz^2 f) FD outer, PS inner
  { // {{{

    SmartFieldMat A_half( NIcomp, NIcomp );
    SmartFieldVec d3f( NIcomp );
    SmartFieldMat tmpmat( NIcomp, NIcomp );
    SmartFieldVec tmpvec( NIcomp );

    Vec2dReal coeffs_outer(0, Vec1dReal(0, 0.));
    Vec2dReal coeffs_inner(0, Vec1dReal(0, 0.));
    Vec2dReal stencil_outer(0, Vec1dReal(0, 0));
    Vec2dReal stencil_inner(0, Vec1dReal(0, 0));

    // get 1st derivative stencils
    d1f_half( stencil_outer, coeffs_outer );
    d1f_half( stencil_inner, coeffs_inner );

    if (stencil_inner[0].size() <= 0)
    {
      std::cout << "***(Error) in Operators_FD::Del_A_Del3_f_ex***\n";
      std::cout << "Stencil(s) allocated improperly.\n";
      exit(1);
    }
    
    Vec1dUInt iphalf(stencil_inner[0].size(), 0);
    Vec1dUInt imhalf(stencil_inner[0].size(), 0);

    // get del2 for inner part
    SmartFieldVec del2_f( f );
    del2_f *= _Del2;

    for (int i=0; i<stencil_outer.size(); i++)
    {
      // get d3f
      tmpvec.setflag_inrealspace(false);
      tmpvec = FieldType(0.);
      d3f.setflag_inrealspace(false);
      d3f = FieldType(0.);

      for (int j=0; j<stencil_inner.size(); j++)
      {
        iphalf = add_stencils(stencil_outer[i], stencil_inner[j]);  
        tmpvec.subscript( del2_f, iphalf );
        d3f.xpby_inplace( tmpvec, coeffs_inner[j] );
      }

      // get A_half
      A_half.setflag_inrealspace(false);
      iphalf = add_stencils( stencil_outer[i], 0.5 );
      imhalf = add_stencils( stencil_outer[i], -0.5 );
      A_half.subscript( A, iphalf );
      tmpmat.subscript( A, imhalf );
      A_half += tmpmat;
      A_half *= FieldType(0.5);

      // take dot product in real space
      d3f.ifft();
      tmpvec.setflag_inrealspace(true);
      tmpvec.dot(A_half, d3f);
      tmpvec.fft();

      // use coeffs to get outer derivative
      dAd3f.xpby_inplace( tmpvec, coeffs_outer[i] );

    }
  } // }}}

  // del_yz .( A del_yz d^2 f /dx^2) -- PS outer, FD inner
  { // {{{

    SmartFieldVec tmpvec( NIcomp );

    // get d^2/dx^2 f
    Vec2dReal coeffs(0, Vec1dReal(0, 0.));
    Vec2dUInt stencil(0, Vec1dUInt(0, 0));
    d2f( stencil, coeffs );

    SmartFieldVec d2f_dx2( NIcomp );
    d2f_dx2.setflag_inrealspace( false );
    d2f_dx2 = FieldType(0.);
    
    tmpvec.setflag_inrealspace( false );
    for (int j = 0; j < stencil.size(); j++)
    {
      tmpvec.subscript( f, stencil[j] ); // tmp = f[i][stencil[j]];
      d2f_dx2.xpby_inplace( tmpvec, coeffs[j] );
    }

    SmartFieldMat Del3_f(NIcomp, _CurrGrid->PSDim() );
    Del3_f.setflag_inrealspace(false);
    Del3_f.outer( d2f_dx2, _Del ); // outer product over PSDim

    SmartFieldMat A_del3_f( NIcomp, _CurrGrid->PSDim() );
    Del3_f.ifft();
    A_del3_f.setflag_inrealspace(true);
    A_del3_f.dot(A, Del3_f);
    A_del3_f.fft();

    tmpvec.setflag_inrealspace(false);
    tmpvec.dot( A_del3_f, _Del );
    dAd3f += tmpvec;

  } // }}}

  // del_yz .( A del_yz del_yz^2 f) -- Pure PS
  { // {{{

    SmartFieldVec tmpvec( NIcomp );
    tmpvec.setflag_inrealspace(false);

    tmpvec = f;
    tmpvec *= _Del2;

    SmartFieldMat Del3_f(NIcomp, _CurrGrid->PSDim() );
    Del3_f.setflag_inrealspace(false);

    Del3_f.outer( tmpvec, _Del ); // outer product over PSDim

    SmartFieldMat A_del3_f( NIcomp, _CurrGrid->PSDim() );

    Del3_f.ifft();
    A_del3_f.setflag_inrealspace(true);
    A_del3_f.dot(A, Del3_f);
    A_del3_f.fft();

    tmpvec.setflag_inrealspace(false);
    tmpvec.dot(A_del3_f, _Del); // inner product over PSDim
    
    dAd3f += tmpvec;

  } // }}}

} // }}}

void Operators_FD::transverse_f_ex ( const SmartFieldVec &f, 
                                           SmartFieldVec & T_f )
{ // {{{

  // f: Rows = _Dim
  // T_f: Rows = _Dim
  //
  // transverse = del^2 - del del
  // ** this may help with computing the velocity
  //    without so much round off error in Model H **
  //
  // ( del^2 - del del ) f = 
  //   [  dyy + dzz   -dxy      -dxz    ] [ fx ] 
  //   [    -dyx    dxx + dzz   -dyz    ] [ fy ] =
  //   [    -dzx      -dzy    dyy + dzz ] [ fz ]
  //
  //   [  (dyy + dzz) fx -    dxy fy      -    dxz fz      ]
  //   [    -dyx fx      + (dxx + dzz) fy -    dyz fz      ]
  //   [    -dzx fx      -    dzy fy      + (dxx + dyy) fz ]
  //

  if ( f.getflag_inrealspace() )
  {
    std::cout << "***(Error) Operators_FD::transverse_f_ex***\n";
    std::cout << "Cannot perform gradient operations ";
    std::cout << "on a field in real space.\n";
    exit(1);
  }

  if ( f.get_Nelem() != _CurrGrid->Dim() )
  {
    std::cout << "***(Error) Operators_FD::transverse_f_ex***\n";
    std::cout << "Input vector must have Dim = ";
    std::cout << _CurrGrid->Dim() << "dimensions\n";
    exit(1);
  }
  
  T_f.setflag_inrealspace(false);

  // For FD coeffs: (stencil width x Nx)
  Vec2dReal d1_coeffs(0, Vec1dReal(0, 0.));
  Vec2dUInt d1_stencil(0, Vec1dUInt(0, 0));
  d1f( d1_stencil, d1_coeffs );

  Vec2dReal d2_coeffs(0, Vec1dReal(0, 0.));
  Vec2dUInt d2_stencil(0, Vec1dUInt(0, 0));
  d2f( d2_stencil, d2_coeffs );

  // for PS derivatives: fyy, fyz, fzy, fzz
  SmartFieldMat DelDel( _Del.get_Nelem(), _Del.get_Nelem() );
  DelDel.outer( _Del, _Del );

  SmartField tmp;
  tmp.setflag_inrealspace(false);

  // --- x-component ---
  // (T_f)_x = (dyy + dzz) fx - dxy fy - dxz fz

  T_f[0].zero();

  tmp = _Del2; // dyy + dzz
  tmp *= f[0]; // (dyy + dzz) fx
  T_f[0] += tmp;

  FD_der( d1_stencil, d1_coeffs, f[1], tmp ); // dx fy 
  tmp *= _Del[0]; // dy (dx fy)
  T_f[0] -= tmp;
  
  if (_CurrGrid->Dim()>2)
  {
    FD_der( d1_stencil, d1_coeffs, f[2], tmp ); // dx fz 
    tmp *= _Del[1]; // dz (dx fy)
    T_f[0] -= tmp;
  }
  
  // --- y-component ---
  // (T_f)_y = - dyx fx + (dxx + dzz) fy - dyz fz

  T_f[1].zero();

  FD_der( d1_stencil, d1_coeffs, f[0], tmp ); // dx fx
  tmp *= _Del[0]; // dy (dx fx)
  T_f[1] -= tmp;

  FD_der( d2_stencil, d2_coeffs, f[1], tmp ); // dxx fy
  T_f[1] += tmp;

  if (_CurrGrid->Dim() > 2)
  {
    tmp = DelDel(1,1); // dzz
    tmp *= f[1];  // dzz fy
    T_f[1] += tmp;
    
    tmp = DelDel(0,1); // dyz
    tmp *= f[2];  // - dyz fz
    T_f[1] -= tmp;
  }

  // --- z-component ---
  // (T_f)_z = - dzx fx - dzy fy + (dxx + dyy) fz

  if (_CurrGrid->Dim() > 2)
  {
    T_f[2].zero();

    FD_der( d1_stencil, d1_coeffs, f[0], tmp ); // dx fx
    tmp *= _Del[1]; // dz (dx fx)
    T_f[2] -= tmp;

    tmp = DelDel(1,0); // dzy
    tmp *= f[1];  // dzy fy
    T_f[2] -= tmp;

    FD_der( d2_stencil, d2_coeffs, f[2], tmp ); // dxx fz
    T_f[2] += tmp;

    tmp = DelDel(0,0); // dyy
    tmp *= f[2];  // dyy fz
    T_f[2] += tmp;
  }

} // }}}

// ----------------- Implicit Operators -----------------

void Operators_FD::Eye_im( SmartFieldOpMat & T )
{ // {{{

  if (T.getflag_inrealspace())
  {
    std::cout << "***(Error) in Operators_FD::Eye_im() ***\n";
    std::cout << "The implicit operator must be Fourier space.\n";
    exit(1);
  }

  SmartFieldMat Eye( T.Nrow(), T.Ncol() );
  Eye.setflag_inrealspace(false);
  Eye.eye( FieldType(1.) );
  T.AddBand(0, Eye);

} // }}}

void Operators_FD::Del2_im( SmartFieldOp & T )
{ // {{{

  // del^2 = d^2/dx^2 + del_yz^2

  if (T.getflag_inrealspace())
  {
    std::cout << "***(Error) in Operators_FD::Eye_im() ***\n";
    std::cout << "The implicit operator must be Fourier space.\n";
    exit(1);
  }

  // A . (d^2/dx^2)
  Vec1dFieldType coeffs(0, 0.);
  Vec1dInt stencil(0, 0);
  d2f_per( stencil, coeffs );

  SmartField band;

  // d^2/dx^2
  for (int i=0; i<stencil.size(); i++)
  {
    band.setflag_inrealspace(false);
    band = coeffs[i];
    T.AddBand(stencil[i], band);
  }

  // del_yz^2
  T.AddBand(0, _Del2);

} // }}}

void Operators_FD::Coulomb_im( SmartFieldOp & T )
{ // {{{

  // del^2 = d^2/dx^2 + del_yz^2

  std::cout << "*** Caveat emptor: Operators_FD::Coulomb_im() was called ***\n";
  std::cout << "***    It is not tested and the output may be garbage    ***\n";

  if (T.getflag_inrealspace())
  {
    std::cout << "***(Error) in Operators_FD::Coulomb_im() ***\n";
    std::cout << "The implicit operator must be Fourier space.\n";
    exit(1);
  }

  SmartField Coul_Op;
  Coul_Op.setflag_inrealspace(false);

  SmartField impulse;
  impulse = FieldType(-1.);
  impulse.setflag_inrealspace(false);

  // Del^2 operator
  SmartFieldOp Del2;
  Del2_im( Del2 );

  // Solve: Del2 Coul_f = -1 in Fourier space
  Solve_im_matrix( Del2, impulse, Coul_Op, 0.);

  // del_yz^2
  T.AddBand(0, Coul_Op);

} // }}}

void Operators_FD::A_Del2_im( const Vec2dFieldType &A, SmartFieldOpMat & T )
{ // {{{

  // del^2 = d^2/dx^2 + del_yz^2

  if (T.getflag_inrealspace())
  {
    std::cout << "***(Error) in Operators_FD::Eye_im() ***\n";
    std::cout << "The implicit operator must be Fourier space.\n";
    exit(1);
  }

  // A . (d^2/dx^2)
  Vec1dFieldType coeffs(0, 0.);
  Vec1dInt stencil(0, 0);
  d2f_per( stencil, coeffs );

  SmartFieldMat tmpmat( T.Nrow(), T.Ncol() );
  SmartFieldMat band( T.Nrow(), T.Ncol() );

  for (int i=0; i<stencil.size(); i++)
  {
    // make diagonal field matrix with FD coeffs
    tmpmat.setflag_inrealspace(false);
    tmpmat.eye( coeffs[i] );

    // dot product with A, to get band
    band.setflag_inrealspace(false);
    band.dot( A, tmpmat );
    T.AddBand(stencil[i], band);
  }

  // A . del_yz^2
  band.setflag_inrealspace(false);
  band = A;
  band *= _Del2;
  T.AddBand(0, band);

} // }}}

void Operators_FD::A_Coulomb_im( const Vec2dFieldType &A, SmartFieldOpMat & T )
{ // {{{

  // del^2 = d^2/dx^2 + del_yz^2

  std::cout << "*** Caveat emptor: Operators_FD::A_Coulomb_im() was called ***\n";
  std::cout << "***     It is not tested and the output may be garbage     ***\n";

  if (T.getflag_inrealspace())
  {
    std::cout << "***(Error) in Operators_FD::A_Coulomb_im() ***\n";
    std::cout << "The implicit operator must be Fourier space.\n";
    exit(1);
  }

  SmartField Coul_Op;
  Coul_Op.setflag_inrealspace(false);

  SmartField impulse;
  impulse = FieldType(-1.);
  impulse.setflag_inrealspace(false);

  // Del^2 operator
  SmartFieldOp Del2;
  Del2_im( Del2 );

  // Solve: Del2 Coul_f = -1 in Fourier space
  Solve_im_matrix( Del2, impulse, Coul_Op, 0.);

  // Get A_Coul = A * Coul_Op 
  SmartFieldMat A_Coul( T.Nrow(), T.Ncol() );

  A_Coul.setflag_inrealspace(false);
  A_Coul = A;
  A_Coul *= Coul_Op;
  T.AddBand(0, A_Coul);

} // }}}

void Operators_FD::A_Del4_im( const Vec2dFieldType &A, SmartFieldOpMat & T )
{ // {{{

  // del^2 del^2 f =  d^4/dx^4 f + 2*d^2/dx^2 del_yz^2 f + del_yz^4 f]

  if (T.getflag_inrealspace())
  {
    std::cout << "***(Error) in Operators_FD::Eye_im() ***\n";
    std::cout << "The implicit operator must be Fourier space.\n";
    exit(1);
  }

  // d^4/dx^4 f
  Vec1dFieldType coeffs (0, 0.);
  Vec1dInt stencil(0, 0);
  d4f_per( stencil, coeffs );

  SmartFieldMat tmpmat( T.Nrow(), T.Ncol() );
  SmartFieldMat band( T.Nrow(), T.Ncol() );

  for (int i=0; i<stencil.size(); i++)
  {
    // make diagonal field matrix with FD coeffs
    tmpmat.setflag_inrealspace(false);
    tmpmat.eye( coeffs[i] );

    // dot product with A, to get band
    band.setflag_inrealspace(false);
    band.dot( A, tmpmat );
    T.AddBand(stencil[i], band);
  }

  // 2*d^2/dx^2 del_yz^2 f
  d2f_per( stencil, coeffs );

  for (int i=0; i<stencil.size(); i++)
  {
    // make diagonal field matrix with FD coeffs
    tmpmat.setflag_inrealspace(false);
    tmpmat.eye( coeffs[i] );

    // dot product with A, to get band
    band.setflag_inrealspace(false);
    band.dot( A, tmpmat );
    band *= _Del2;
    band *= FieldType(2.);
    T.AddBand(stencil[i], band);
  }

  // del_yz^4 f]
  band.setflag_inrealspace(false);
  band = A;
  band *= _Del2;
  band *= _Del2;
  T.AddBand(0, band);

} // }}}

void Operators_FD::A_Del2_F_im( const Vec2dFieldType &A, 
                                const SmartFieldOpMat &F, 
                                      SmartFieldOpMat & T )
{ // {{{

  SmartFieldMat old_band( T.Nrow(), T.Ncol() );
  SmartFieldMat new_band( T.Nrow(), T.Ncol() );
  SmartFieldMat A_Del2( T.Nrow(), T.Ncol() );

  Vec2dFieldType tmpmat( T.Nrow(), Vec1dFieldType(T.Ncol(), FieldType(0.)) );

  std::vector< int > band_stencil = F.GetBandStencil();

  if (T.getflag_inrealspace())
  {
    std::cout << "***(Error) in Operators_FD::Eye_im() ***\n";
    std::cout << "The implicit operator must be Fourier space.\n";
    exit(1);
  }

  // A . (d^2/dx^2)
  Vec1dFieldType d2_coeffs(0, 0.);
  Vec1dInt d2_stencil(0, 0);
  d2f_per( d2_stencil, d2_coeffs );

  for ( int i=0; i<F.Nband(); i++ )
  {

    old_band.setflag_inrealspace(false);
    F.GetBand( band_stencil[i], old_band );

    for (int j=0; j<d2_stencil.size(); j++)
    {

      // make diagonal matrix with FD coeffs
      tmpmat = eye( T.Nrow(), T.Ncol(), d2_coeffs[j] );

      // dot product with A, to get band
      A_Del2.setflag_inrealspace(false);
      A_Del2 = matmul( A, tmpmat );
      new_band.setflag_inrealspace(false);
      new_band.dot( A_Del2, old_band );
      T.AddBand(band_stencil[i] + d2_stencil[j], new_band);
    }
  }

  // A . del_yz^2 F
  for ( int i=0; i<F.Nband(); i++ )
  {
    old_band.setflag_inrealspace(false);
    F.GetBand( band_stencil[i], old_band );
    old_band *= _Del2;

    new_band.setflag_inrealspace(false);
    new_band.dot( A, old_band );
    T.AddBand( band_stencil[i], new_band );
  }

} // }}}

// ----------------- Solvers -----------------

void Operators_FD::Solve_im_matrix( SmartFieldOpMat & A, 
                                    SmartFieldVec & rhs, 
                                    SmartFieldVec & x,
                                    RealType t )
{ // {{{

  // Solve A.x = rhs
  // t is the simulation time, needed for time-dependent BCs

  _CurrBCs->solve_BCs( A, rhs, x, t );

} // }}}

void Operators_FD::Solve_im_matrix( SmartFieldOp & A, 
                                    SmartField & rhs, 
                                    SmartField & x,
                                    RealType t )
{ // {{{

  // Create dummy matrices and vectors to pass to the other
  // solver function. This avoids having to re-write the 
  // whole BCs structure to accomodate single component 
  // solvers.
  SmartFieldOpMat A_mat(1, 1);
  SmartFieldVec rhs_vec(1);
  SmartFieldVec x_vec(1);

  A_mat.setflag_inrealspace( A.getflag_inrealspace() );
  rhs_vec.setflag_inrealspace( rhs.getflag_inrealspace() );
  x_vec.setflag_inrealspace( x.getflag_inrealspace() );

  A_mat[0][0] = A;
  rhs_vec[0] = rhs;
  x_vec[0] = x;

  Solve_im_matrix( A_mat, rhs_vec, x_vec, t );

  x = x_vec[0];

} // }}}

// ----------------- Finite Difference Formulas -----------------

// safe-end FD schemes
void Operators_FD::d1f( Vec2dUInt & stencil, Vec2dReal & coeffs  )
{ // {{{

  // df/dx: second-order finite-difference
  // df/dx_0 = (-3/2 f_0     + 2 f_1     - 1/2 f_2     )/dx
  // df/dx_i = (-1/2 f_{i-1}             + 1/2 f_{i+1} )/dx
  // df/dx_n = ( 1/2 f_{n-2} - 2 f_{n-1} + 3/2 f_{n}   )/dx

  RealType dx = _CurrGrid->Dx();
  int Nx = _CurrGrid->FDSize();

  stencil.resize( 3, Vec1dUInt(Nx, 0) );
  coeffs.resize( 3, Vec1dReal(Nx, 0.));

  // --- Periodic ---

  for (int m=0; m<Nx; m++)
  {
    stencil[0][m] = (m-1+Nx)%Nx;
    coeffs[0][m] = -1./(2*dx);

    stencil[1][m] = m;
    coeffs[1][m] =  0.;
  
    stencil[2][m] = (m+1)%Nx;
    coeffs[2][m] =  1./(2*dx);
  }

  // --- One-sided ends ---

//  // left node
//  stencil[0][0] = 0; // m = 0
//  coeffs[0][0] = -3./(2.*dx);
//  for (int m=1; m<Nx-1; m++)
//  {
//    stencil[0][m] = m-1;
//    coeffs[0][m] = -1./(2*dx);
//  }
//  stencil[0][Nx-1] = (Nx-1)-2; // m = Nx-1
//  coeffs[0][Nx-1] =  1./(2.*dx);
//
//  // middle node
//  stencil[1][0] = 1;
//  coeffs[1][0] =  2./dx;
//  for (int m = 1; m<Nx-1; m++)
//  {
//    stencil[1][m] = m;
//    coeffs[1][m] =  0.;
//  }
//  stencil[1][Nx-1] = (Nx-1)-1;
//  coeffs[1][Nx-1] = -2./dx;
//  
//  // right node
//  stencil[2][0] = 2;
//  coeffs[2][0] = -1./(2.*dx);
//  for (int m = 1; m<Nx-1; m++)
//  {
//    stencil[2][m] = m+1;
//    coeffs[2][m] =  1./(2*dx);
//  }
//  stencil[2][Nx-1] = (Nx-1)-0;
//  coeffs[2][Nx-1] =  3./(2.*dx);

} // }}}

void Operators_FD::d2f( Vec2dUInt & stencil, Vec2dReal & coeffs )
{ // {{{

  // d2f/dx2: second-order finite-difference
  // d2f/dx2_0 = (2 f_0 - 5 f_1 + 4 f_2 - f_3)/dx2
  // d2f/dx2_i = (- f_{i-1} + 2 f_{i} - f_{i+1})/dx2
  // d2f/dx2_n = (- f_{n-3} + 4 f_{n-2} - 5 f_{n-1} + 2 f_{n})/dx2
  
  RealType dx2 = _CurrGrid->Dx() * _CurrGrid->Dx();
  int Nx = _CurrGrid->FDSize();

  stencil.resize( 4, Vec1dUInt(Nx, 0) );
  coeffs.resize( 4, Vec1dReal(Nx, 0.));

  // --- Periodic ---

  for (int m=0; m<Nx; m++)
  {
    stencil[0][m] = (m-1+Nx)%Nx;
    coeffs[0][m] = 1./dx2;

    stencil[1][m] = m;
    coeffs[1][m] = -2./dx2;

    stencil[2][m] = (m+1)%Nx;
    coeffs[2][m] = 1./dx2;
  }

  // --- One-sided ends ---

//  // node: 0 (left-most)
//  stencil[0][0] = 0; // m = 0
//  coeffs[0][0] = 2./dx2;
//  for (int m=1; m<Nx-1; m++)
//  {
//    stencil[0][m] = m-1;
//    coeffs[0][m] = 1./dx2;
//  }
//  stencil[0][Nx-1] = (Nx-1)-3; // m = Nx-1
//  coeffs[0][Nx-1] =  -1./dx2;
//
//  // node 1: 
//  stencil[1][0] = 1; // m = 0
//  coeffs[1][0] = -5./dx2;
//  for (int m=1; m<Nx-1; m++)
//  {
//    stencil[1][m] = m;
//    coeffs[1][m] = -2./dx2;
//  }
//  stencil[1][Nx-1] = (Nx-1)-2; // m = Nx-1
//  coeffs[1][Nx-1] =  4./dx2;
//
//  // node 2:
//  stencil[2][0] = 2; // m = 0
//  coeffs[2][0] = 4./dx2;
//  for (int m=1; m<Nx-1; m++)
//  {
//    stencil[2][m] = m+1;
//    coeffs[2][m] = 1./dx2;
//  }
//  stencil[2][Nx-1] = (Nx-1)-1; // m = Nx-1
//  coeffs[2][Nx-1] =  -5./dx2;
//
//  // node 3: right most
//  stencil[3][0] = 3; // m = 0
//  coeffs[3][0] = -1./dx2;
//  for (int m=1; m<Nx-1; m++)
//  {
//    stencil[3][m] = m;
//    coeffs[3][m] = 0.;
//  }
//  stencil[3][Nx-1] = (Nx-1); // m = Nx-1
//  coeffs[3][Nx-1] =  2./dx2;

} // }}}

void Operators_FD::d3f( Vec2dUInt & stencil, Vec2dReal & coeffs )
{ // {{{

  // d3f/dx3: second-order finite-difference
  // d3f/dx3_0     = (-5/2 f_0     + 9 f_1     -12 f_2      + 7 f_3     - 3/2 f_4     )/dx3
  // d3f/dx3_1     = (-3/2 f_0     + 5 f_1     - 6 f_2      + 3 f_3     - 1/2 f_4     )/dx3
  // d3f/dx3_{i}   = (-1/2 f_{i-2} + f_{i-1}                - f_{i+1}   + 1/2 f_{i+2} )/dx3
  // d3f/dx3_{n-1} = ( 1/2 f_{n-4} - 3 f_{n-3} + 6 f_{n-2}  - 5 f_{n-1} + 3/2 f_{n}   )/dx3
  // d3f/dx3_{n}   = ( 3/2 f_{n-4} - 7 f_{n-3} + 12 f_{n-2} - 9 f_{n-1} + 5/2 f_{n}   )/dx3
  
  RealType dx3 = pow(_CurrGrid->Dx(), 3);
  int Nx = _CurrGrid->FDSize();

  stencil.resize( 5, Vec1dUInt(Nx, 0) );
  coeffs.resize( 5, Vec1dReal(Nx, 0.));

  // --- Periodic ---

  for (int m=0; m<Nx; m++)
  {
    stencil[0][m] = (m-2+Nx)%Nx;
    coeffs[0][m] = -1./(2*dx3);

    stencil[1][m] = (m-1+Nx)%Nx;
    coeffs[1][m] = 1./dx3;

    stencil[2][m] = m;
    coeffs[2][m] = 0.;

    stencil[3][m] = (m+1)%Nx;
    coeffs[3][m] = -1./dx3;

    stencil[4][m] = (m+2)%Nx;
    coeffs[4][m] = 1./(2*dx3);
  }

  // --- One-sided ends ---

//  // node: 0 (left-most)
//  stencil[0][0] = 0; // m = 0
//  coeffs[0][0] = -5./(2*dx3);
//  stencil[0][1] = 0; // m = 1
//  coeffs[0][1] = -3./(2*dx3);
//  for (int m=2; m<Nx-2; m++)
//  {
//    stencil[0][m] = m-2;
//    coeffs[0][m] = -1./(2*dx3);
//  }
//  stencil[0][Nx-2] = (Nx-2)-3; // m = Nx-2
//  coeffs[0][Nx-2] = 1./(2*dx3);
//  stencil[0][Nx-1] = (Nx-1)-4; // m = Nx-1
//  coeffs[0][Nx-1] =  3./(2*dx3);
//
//  // node: 1
//  stencil[1][0] = 1; // m = 0
//  coeffs[1][0] = 9./dx3;
//  stencil[1][1] = 1; // m = 1
//  coeffs[1][1] = 5./dx3;
//  for (int m=2; m<Nx-2; m++)
//  {
//    stencil[1][m] = m-1;
//    coeffs[1][m] = 1./dx3;
//  }
//  stencil[1][Nx-2] = (Nx-2)-2; // m = Nx-2
//  coeffs[1][Nx-2] = -3./dx3;
//  stencil[1][Nx-1] = (Nx-1)-3; // m = Nx-1
//  coeffs[1][Nx-1] =  -7./dx3;
//
//  // node: 2 (center node)
//  stencil[2][0] = 2; // m = 0
//  coeffs[2][0] = -12./dx3;
//  stencil[2][1] = 2; // m = 1
//  coeffs[2][1] = -6./dx3;
//  for (int m=2; m<Nx-2; m++)
//  {
//    stencil[2][m] = m;
//    coeffs[2][m] = 0.;
//  }
//  stencil[2][Nx-2] = (Nx-2)-1; // m = Nx-2
//  coeffs[2][Nx-2] = 6./dx3;
//  stencil[2][Nx-1] = (Nx-1)-2; // m = Nx-1
//  coeffs[2][Nx-1] = 12./dx3;
//
//  // node: 3
//  stencil[3][0] = 3; // m = 0
//  coeffs[3][0] = 7./dx3;
//  stencil[3][1] = 3; // m = 1
//  coeffs[3][1] = 3./dx3;
//  for (int m=2; m<Nx-2; m++)
//  {
//    stencil[3][m] = m+1;
//    coeffs[3][m] = -1./dx3;
//  }
//  stencil[3][Nx-2] = (Nx-2); // m = Nx-2
//  coeffs[3][Nx-2] = -5./dx3;
//  stencil[3][Nx-1] = (Nx-1)-1; // m = Nx-1
//  coeffs[3][Nx-1] = -9./dx3;
//
//  // node: 4 (right most)
//  stencil[4][0] = 4; // m = 0
//  coeffs[4][0] = -3./(2.*dx3);
//  stencil[4][1] = 4; // m = 1
//  coeffs[4][1] = -1./(2*dx3);
//  for (int m=2; m<Nx-2; m++)
//  {
//    stencil[4][m] = m+2;
//    coeffs[4][m] = 1./(2*dx3);
//  }
//  stencil[4][Nx-2] = (Nx-2)+1; // m = Nx-2
//  coeffs[4][Nx-2] = 3./(2*dx3);
//  stencil[4][Nx-1] = (Nx-1); // m = Nx-1
//  coeffs[4][Nx-1] = 5./(2*dx3);

} // }}}

void Operators_FD::d4f( Vec2dUInt & stencil, Vec2dReal & coeffs )
{ // {{{

  // d4f/dx4: second-order finite-difference
  // d4f/dx4_0     = (  3 f_0     - 14 f_1     + 26 f_2      - 24 f_3     + 11 f_4     - 2 f_5   )/dx4
  // d4f/dx4_1     = (  2 f_0     -  9 f_1     + 16 f_2      - 14 f_3     +  6 f_4     - f_5     )/dx4
  // d4f/dx4_{i}   = (    f_{i-2} -  4 f_{i-1} +  6 f_{i}    -  4 f_{i+1} +    f_{i+2}           )/dx4
  // d4f/dx4_{n-1} = ( -1 f_{n-5} +  6 f_{n-4} - 14 f_{n-3}  + 16 f_{n-2} -  9 f_{n-1} + 2 f_{n} )/dx4
  // d4f/dx4_{n}   = ( -2 f_{n-5} + 11 f_{n-4} - 24 f_{n-3}  + 26 f_{n-2} - 14 f_{n-1} + 3 f_{n} )/dx4
  
  RealType dx4 = pow(_CurrGrid->Dx(), 4);
  int Nx = _CurrGrid->FDSize();

  stencil.resize( 6, Vec1dUInt(Nx, 0) );
  coeffs.resize( 6, Vec1dReal(Nx, 0.));

  // --- Periodic ---

  for (int m=0; m<Nx; m++)
  {

    stencil[0][m] = (m-2+Nx)%Nx;
    coeffs[0][m] = 1./dx4;

    stencil[1][m] = (m-1+Nx)%Nx;
    coeffs[1][m] = -4./dx4;

    stencil[2][m] = m;
    coeffs[2][m] = 6./dx4;

    stencil[3][m] = (m+1)%Nx;
    coeffs[3][m] = -4./dx4;

    stencil[4][m] = (m+2)%Nx;
    coeffs[4][m] = 1./dx4;

    stencil[5][m] = (m+3)%Nx;
    coeffs[5][m] = 0.;

  }

  // --- One-sided ends ---

//  // node: 0 (left-most)
//  stencil[0][0] = 0; // m = 0
//  coeffs[0][0] = 3./dx4;
//  stencil[0][1] = 0; // m = 1
//  coeffs[0][1] = 2./dx4;
//  for (int m=2; m<Nx-2; m++)
//  {
//    stencil[0][m] = m-2;
//    coeffs[0][m] = 1./dx4;
//  }
//  stencil[0][Nx-2] = Nx-6; // m = Nx-2
//  coeffs[0][Nx-2] = -1./dx4;
//  stencil[0][Nx-1] = Nx-6; // m = Nx-1
//  coeffs[0][Nx-1] = -2./dx4;
//
//  // node: 1 
//  stencil[1][0] = 1; // m = 0
//  coeffs[1][0] = -14./dx4;
//  stencil[1][1] = 1; // m = 1
//  coeffs[1][1] = -9./dx4;
//  for (int m=2; m<Nx-2; m++)
//  {
//    stencil[1][m] = m-1;
//    coeffs[1][m] = -4./dx4;
//  }
//  stencil[1][Nx-2] = Nx-5; // m = Nx-2
//  coeffs[1][Nx-2] = 6./dx4;
//  stencil[1][Nx-1] = Nx-5; // m = Nx-1
//  coeffs[1][Nx-1] = 11./dx4;
//
//  // node: 2 (center-left node)
//  stencil[2][0] = 2; // m = 0
//  coeffs[2][0] = 26./dx4;
//  stencil[2][1] = 2; // m = 1
//  coeffs[2][1] = 16./dx4;
//  for (int m=2; m<Nx-2; m++)
//  {
//    stencil[2][m] = m;
//    coeffs[2][m] = 6./dx4;
//  }
//  stencil[2][Nx-2] = Nx-4; // m = Nx-2
//  coeffs[2][Nx-2] = -14./dx4;
//  stencil[2][Nx-1] = Nx-4; // m = Nx-1
//  coeffs[2][Nx-1] = -24./dx4;
//
//  // node: 3 (center-right node)
//  stencil[3][0] = 3; // m = 0
//  coeffs[3][0] = -24./dx4;
//  stencil[3][1] = 3; // m = 1
//  coeffs[3][1] = -14./dx4;
//  for (int m=2; m<Nx-2; m++)
//  {
//    stencil[3][m] = m+1;
//    coeffs[3][m] = -4./dx4;
//  }
//  stencil[3][Nx-2] = Nx-3; // m = Nx-2
//  coeffs[3][Nx-2] = 16./dx4;
//  stencil[3][Nx-1] = Nx-3; // m = Nx-1
//  coeffs[3][Nx-1] = 26./dx4;
//
//  // node: 4
//  stencil[4][0] = 4; // m = 0
//  coeffs[4][0] = 11./dx4;
//  stencil[4][1] = 4; // m = 1
//  coeffs[4][1] = 6./dx4;
//  for (int m=2; m<Nx-2; m++)
//  {
//    stencil[4][m] = m+2;
//    coeffs[4][m] = 1./dx4;
//  }
//  stencil[4][Nx-2] = Nx-2; // m = Nx-2
//  coeffs[4][Nx-2] = -9./dx4;
//  stencil[4][Nx-1] = Nx-2; // m = Nx-1
//  coeffs[4][Nx-1] = -14./dx4;
//
//  // node: 5
//  stencil[5][0] = 5; // m = 0
//  coeffs[5][0] = -2./dx4;
//  stencil[5][1] = 5; // m = 1
//  coeffs[5][1] = -1./dx4;
//  for (int m=2; m<Nx-2; m++)
//  {
//    stencil[5][m] = m;
//    coeffs[5][m] = 0.;
//  }
//  stencil[5][Nx-2] = Nx-1; // m = Nx-2
//  coeffs[5][Nx-2] = 2./dx4;
//  stencil[5][Nx-1] = Nx-1; // m = Nx-1
//  coeffs[5][Nx-1] = 3./dx4;

} // }}}

// half-nodes schemes for nested derivatives
void Operators_FD::d1f_half( Vec2dReal & stencil, Vec2dReal & coeffs  )
{ // {{{

  // df/dx: second-order finite-difference
  // df/dx_0 = (-3/2 f_0       + 2 f_1     - 1/2 f_2       )/dx
  // df/dx_i = (-1   f_{i-1/2}             + 1   f_{i+1/2} )/dx
  // df/dx_n = ( 1/2 f_{n-2}   - 2 f_{n-1} + 3/2 f_{n}     )/dx
  //
  // Note: let edges be on exact nodes, 
  // since we cannot simultaneously 
  // let them be on half nodes and use
  // a self adjoint operator
  //

  RealType dx = _CurrGrid->Dx();
  UInt Nx = _CurrGrid->FDSize();

  stencil.resize( 3, Vec1dReal(Nx, 0.));
  coeffs.resize( 3, Vec1dReal(Nx, 0.));

  // --- Periodic ---

  for (int m=0; m<Nx; m++)
  {
    stencil[0][m] = m-1./2;
    coeffs[0][m] = -1./dx;

    stencil[1][m] = m-1./2;
    coeffs[1][m] =  0.;

    stencil[2][m] = m+1./2;
    coeffs[2][m] =  1./dx;
  }

  // --- One-sided ends ---

//  // left node
//  stencil[0][0] = 0; // m = 0
//  coeffs[0][0] = -3./(2.*dx);
//  for (int m=1; m<Nx-1; m++)
//  {
//    stencil[0][m] = m-1./2;
//    coeffs[0][m] = -1./dx;
//  }
//  stencil[0][Nx-1] = (Nx-1)-2; // m = Nx-1
//  coeffs[0][Nx-1] =  1./(2.*dx);
//
//  // middle node
//  stencil[1][0] = 1;
//  coeffs[1][0] =  2./dx;
//  for (int m = 1; m<Nx-1; m++)
//  {
//    stencil[1][m] = m-1./2;
//    coeffs[1][m] =  0.;
//  }
//  stencil[1][Nx-1] = (Nx-1)-1;
//  coeffs[1][Nx-1] = -2./dx;
//  
//  // right node
//  stencil[2][0] = 2;
//  coeffs[2][0] = -1./(2.*dx);
//  for (int m = 1; m<Nx-1; m++)
//  {
//    stencil[2][m] = m+1./2;
//    coeffs[2][m] =  1./dx;
//  }
//  stencil[2][Nx-1] = (Nx-1)-0;
//  coeffs[2][Nx-1] =  3./(2.*dx);

} // }}}

void Operators_FD::d3f_half( Vec2dReal & stencil, Vec2dReal & coeffs )
{ // {{{

  // d3f/dx3: second-order finite-difference, offset by 1/2
  // d^3f/dx^3_0     = (-5/2 f_0       + 9 f_1       - 12 f_2      + 7 f_3       - 3/2 f_4       )/dx^3
  // d^3f/dx^3_1     = (-2   f_{1/2}   + 7 f_{3/2}   - 9 f_{5/2}   + 5 f_{7/2}   - 1   f_{9/2}   )/dx^3
  // d^3f/dx^3_i     = (-1   f_{i-3/2} + 3 f_{i-1/2}               - 3 f_{i+1/2} + 1   f_{i+3/2} )/dx^3
  // d^3f/dx^3_{n-1} = ( 1   f_{n-9/2} - 5 f_{n-7/2} + 9 f_{n-5/2} - 7 f_{n-3/2} + 2   f_{n-1/2} )/dx^3
  // d^3f/dx^3_n     = ( 3/2 f_{n-4}   - 7 f_{n-3}   + 12 f_{n-2}  - 9 f{n-1}    + 5/2 f_{n}     )/dx^3
  // *** In our case n = Nx-1
  //
  // Note: let edges be on exact nodes, 
  // since we cannot simultaneously 
  // let them be on half nodes and use
  // a self adjoint operator
  //

  RealType dx3 = pow(_CurrGrid->Dx(), 3);

  UInt Nx = _CurrGrid->FDSize();

  stencil.resize( 5, Vec1dReal(Nx, 0) );
  coeffs.resize( 5, Vec1dReal(Nx, 0.));

  // --- Periodic ---

  for (int m=0; m<Nx; m++)
  {

    stencil[0][m] = m-3./2;
    coeffs[0][m] = -1./dx3;

    stencil[1][m] = m-1./2;
    coeffs[1][m] = 3./dx3;

    stencil[2][m] = m+1./2;
    coeffs[2][m] = 0.;

    stencil[3][m] = m+1./2;
    coeffs[3][m] = -3./dx3;

    stencil[4][m] = m+3./2;
    coeffs[4][m] = 1./dx3;

  }

  // --- One-sided ends ---

//  // node: 0 (left-most)
//  stencil[0][0] = 0; // m = 0
//  coeffs[0][0] = -5./(2*dx3);
//  stencil[0][1] = 1./2; // m = 1
//  coeffs[0][1] = -2./dx3;
//  for (int m=2; m<Nx-2; m++)
//  {
//    stencil[0][m] = m-3./2;
//    coeffs[0][m] = -1./dx3;
//  }
//  stencil[0][Nx-2] = (Nx-1)-9./2; // m = Nx-2
//  coeffs[0][Nx-2] = 1./dx3;
//  stencil[0][Nx-1] = (Nx-1)-4; // m = Nx-1
//  coeffs[0][Nx-1] =  3./(2*dx3);
//
//  // node: 1
//  stencil[1][0] = 1; // m = 0
//  coeffs[1][0] = 9./dx3;
//  stencil[1][1] = 3./2; // m = 1
//  coeffs[1][1] = 7./dx3;
//  for (int m=2; m<Nx-2; m++)
//  {
//    stencil[1][m] = m-1./2;
//    coeffs[1][m] = 3./dx3;
//  }
//  stencil[1][Nx-2] = (Nx-1)-7./2; // m = Nx-2
//  coeffs[1][Nx-2] = -5./dx3;
//  stencil[1][Nx-1] = (Nx-1)-3; // m = Nx-1
//  coeffs[1][Nx-1] =  -7./dx3;
//
//  // node: 2 (center node)
//  stencil[2][0] = 2; // m = 0
//  coeffs[2][0] = -12./dx3;
//  stencil[2][1] = 5./2; // m = 1
//  coeffs[2][1] = -9./dx3;
//  for (int m=2; m<Nx-2; m++)
//  {
//    stencil[2][m] = m+1./2;
//    coeffs[2][m] = 0.;
//  }
//  stencil[2][Nx-2] = (Nx-1)-5./2; // m = Nx-2
//  coeffs[2][Nx-2] = 9./dx3;
//  stencil[2][Nx-1] = (Nx-1)-2; // m = Nx-1
//  coeffs[2][Nx-1] = 12./dx3;
//
//  // node: 3
//  stencil[3][0] = 3; // m = 0
//  coeffs[3][0] = 7./dx3;
//  stencil[3][1] = 7./2; // m = 1
//  coeffs[3][1] = 5./dx3;
//  for (int m=2; m<Nx-2; m++)
//  {
//    stencil[3][m] = m+1./2;
//    coeffs[3][m] = -3./dx3;
//  }
//  stencil[3][Nx-2] = (Nx-1)-3./2; // m = Nx-2
//  coeffs[3][Nx-2] = -7./dx3;
//  stencil[3][Nx-1] = (Nx-1)-1; // m = Nx-1
//  coeffs[3][Nx-1] = -9./dx3;
//
//  // node: 4 (right most)
//  stencil[4][0] = 4; // m = 0
//  coeffs[4][0] = -3./(2.*dx3);
//  stencil[4][1] = 9./2; // m = 1
//  coeffs[4][1] = -1./dx3;
//  for (int m=2; m<Nx-2; m++)
//  {
//    stencil[4][m] = m+3./2;
//    coeffs[4][m] = 1./dx3;
//  }
//  stencil[4][Nx-2] = (Nx-1)-1./2; // m = Nx-2
//  coeffs[4][Nx-2] = 2./dx3;
//  stencil[4][Nx-1] = (Nx-1); // m = Nx-1
//  coeffs[4][Nx-1] = 5./(2*dx3);

} // }}}

// periodic schemes for implicit matrices
void Operators_FD::d1f_per( Vec1dInt & stencil, Vec1dFieldType & coeffs  )
{ // {{{

  // df/dx: second-order finite-difference
  // df/dx_i = (-1/2 f_{i-1} + 1/2 f_{i+1})/dx

  RealType dx = _CurrGrid->Dx();

  stencil.resize( 2, 0 );
  coeffs.resize( 2, 0.);

  stencil[0] = -1;
  stencil[1] = 1;

  coeffs[0] = -1./(2*dx);
  coeffs[1] =  1./(2*dx);

} // }}}

void Operators_FD::d2f_per( Vec1dInt & stencil, Vec1dFieldType & coeffs )
{ // {{{

  // d2f/dx2: second-order finite-difference
  // d2f/dx2_i = (f_{i-1} - 2 f_{i} + f_{i+1})/dx2
  
  RealType dx2 = _CurrGrid->Dx() * _CurrGrid->Dx();

  stencil.resize( 3, 0 );
  coeffs.resize( 3, 0. );

  stencil[0] = -1;
  stencil[1] = 0;
  stencil[2] = 1;

  coeffs[0] = 1./dx2;
  coeffs[1] = -2./dx2;
  coeffs[2] = 1./dx2;

} // }}}

void Operators_FD::d3f_per( Vec1dInt & stencil, Vec1dFieldType & coeffs )
{ // {{{

  // d3f/dx3: second-order finite-difference
  // d3f/dx3_{i} = (-1/2 f_{i-2} + f_{i-1} - f_{i+1} + 1/2 f_{i+2} )/dx3
  
  RealType dx3 = pow(_CurrGrid->Dx(), 3);

  stencil.resize( 4, 0 );
  coeffs.resize( 4, 0. );

  stencil[0] = -2;
  stencil[1] = -1;
  stencil[2] =  1;
  stencil[3] =  2;

  coeffs[0] = -1./(2*dx3);
  coeffs[1] = 1./dx3;
  coeffs[2] = -1./dx3;
  coeffs[3] = 1./(2*dx3);

} // }}}

void Operators_FD::d4f_per( Vec1dInt & stencil, Vec1dFieldType & coeffs )
{ // {{{

  // d4f/dx4: second-order finite-difference
  // d4f/dx4_{i} = ( f_{i-2} - 4 f_{i-1} + 6 f_{i} - 4 f_{i+1} + f_{i+2} )/dx4
  
  RealType dx4 = pow(_CurrGrid->Dx(), 4);

  stencil.resize( 5, 0 );
  coeffs.resize( 5, 0. );

  stencil[0] = -2;
  stencil[1] = -1;
  stencil[2] =  0;
  stencil[3] =  1;
  stencil[4] =  2;

  coeffs[0] =  1./dx4;
  coeffs[1] = -4./dx4;
  coeffs[2] =  6./dx4;
  coeffs[3] = -4./dx4;
  coeffs[4] =  1./dx4;

} // }}}

// useful stencil operations
Vec1dUInt Operators_FD::add_stencils( Vec1dReal & stencil_a, 
                                      Vec1dReal & stencil_b )
{ // {{{

  if ( stencil_a.size() != stencil_b.size() )
  {
    std::cout << "***(Error) Operators_FD::add_stencils***\n";
    std::cout << "Stencil sizes do not match\n";
    exit(1);
  }

  Vec1dUInt lhs( stencil_a.size(), 0 );
  RealType idx;

  for (UInt i=0; i<stencil_a.size(); i++)
  {
    idx = stencil_a[i] + stencil_b[i] - i;
    // (begin) added for periodic FD stencils
    if (idx < 0) idx += stencil_a.size();
    if (idx > stencil_a.size()-1) idx -= stencil_a.size();
    // (end) added for periodic FD stencils
    if (std::floor(idx) == idx)
    {
      lhs[i] = UInt(idx);
    }
    else
    {
      std::cout << "***(Error) Operators_FD::add_stencils***\n";
      std::cout << "Stencil sum is not an integer\n";
      exit(1);
    }
  }

  return lhs;

} // }}}

Vec1dUInt Operators_FD::add_stencils( Vec1dReal & stencil_a, RealType b )
{ // {{{

  Vec1dUInt lhs( stencil_a.size(), 0 );
  RealType idx;

  // skip beginning and ending nodes
  //for (UInt i=1; i<stencil_a.size()-1; i++)
  for (UInt i=0; i<stencil_a.size(); i++)
  {
    idx = stencil_a[i] + b;
    // (begin) added for periodic FD stencils
    if (idx < 0) idx += stencil_a.size();
    if (idx > stencil_a.size()-1) idx -= stencil_a.size();
    // (end) added for periodic FD stencils
    if (std::floor(idx) == idx)
    {
      lhs[i] = UInt(idx);
    }
    else
    {
      std::cout << "***(Error) Operators_FD::add_stencils***\n";
      std::cout << "Stencil sum is not an integer\n";
      exit(1);
    }
  }
  
  return lhs;

} // }}}

void Operators_FD::FD_der( const Vec2dUInt & stencil, const Vec2dReal & coeffs,
                           const SmartField & f, SmartField & df )
{ // {{{

  // (in) stencil: array of stencil arrays
  // (in) coeffs: array of FD coeff arrays for each stencil
  // (in) f: incoming function
  // (out) df: outgoing FD derivative of some order

  SmartField tmp;
  tmp.setflag_inrealspace(false);
  tmp.zero();

  df.zero();
  for (int j = 0; j < stencil.size(); j++)
  {
    tmp.subscript( f, stencil[j] ); // tmp = f[ stencil[j] ];
    df.xpby_inplace( tmp, coeffs[j] );
  }

} // }}}

