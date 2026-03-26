/***************************************************************
*
* Const Vel boundary condition solver
* 
* DRT -- Mon, 17 Apr 2017
*
****************************************************************/

#include "BCs_Const_Vel.h"

// ----------------- Constructor/Destructor -----------------

BCs_Const_Vel :: BCs_Const_Vel(std::string filename, int nicomp, 
                       Grid * gridarg ) : 
BCs_Base( nicomp, gridarg),
_Dim( gridarg->Dim() )
{ // {{{

  // Read in Boundary conditions
  std::string tmp_str;

  if (jsonFile(filename.c_str()))
  {
    std::ifstream f2(filename.c_str());
    if ( f2.is_open() )
    {
      rapidjson::IStreamWrapper isw(f2);
      rapidjson::Document d;
      d.ParseStream(isw);
      _BCFlag = d["BCFlag"].GetInt();

      _Vel0.resize(_Dim, 0.);
      rapidjson::Value& Vel0_in = d["vel_0"];
      for ( UInt i=0; i<_Dim; i++)
      {
        _Vel0[i] = Vel0_in[i].GetDouble();
      }

      _VelL.resize(_Dim, 0.);
      rapidjson::Value& VelL_in = d["vel_L"];
      for ( UInt i=0; i<_Dim; i++)
      {
        _VelL[i] = VelL_in[i].GetDouble();
      }

    }
    else
    {
      std::cout << "cannot open " << filename << std::endl;
      exit(1);
    }
  }
  else
  {//legacy file format
    std::ifstream f2(filename.c_str());
    if ( f2.is_open() )
    {
      std::getline(f2, tmp_str, '#');
      _BCFlag = atoi(tmp_str.c_str());
      std::getline(f2, tmp_str, '\n');

      _Vel0.resize(_Dim, 0.);
      for ( UInt i=0; i<_Dim; i++)
      {
        std::getline(f2, tmp_str, '#');
        _Vel0[i] = atof(tmp_str.c_str());
        std::getline(f2, tmp_str, '\n');
      }

      _VelL.resize(_Dim, 0.);
      for ( UInt i=0; i<_Dim; i++)
      {
        std::getline(f2, tmp_str, '#');
        _VelL[i] = atof(tmp_str.c_str());
        std::getline(f2, tmp_str, '\n');
      }
    }
    else
    {
      std::cout << "cannot open " << filename << std::endl;
      exit(1);
    }
  }

  std::cout << "\n" << std::scientific;
  std::cout << " * Vel Boundary Conditions Initiated\n" << std::scientific;
  std::cout << "   - BC Parameters:\n" << std::scientific;
  std::cout << "      BCFlag             = " << _BCFlag << std::endl;
  for (int i=0; i<_Dim; ++i)
  {
    std::cout << "        Vel(x=0)[" << i << "]";
    std::cout <<"      = "<< _Vel0[i] << "\n";
  }
  for (int i=0; i<_Dim; ++i)
  {
    std::cout << "        Vel(x=L)[" << i << "]";
    std::cout <<"      = "<< _VelL[i] << "\n";
  }

} // }}}

BCs_Const_Vel :: ~BCs_Const_Vel()
{ // {{{
} // }}}

// ----------------- Constant BCs --------------------

void BCs_Const_Vel::solve_BCs( SmartFieldOpMat & T, 
                               SmartFieldVec & rhs, 
                               SmartFieldVec & x,
                               RealType t )
{ // {{{

  UInt Dim = T.Nrow();
  UInt BW = T.BW();
  UInt Nx = T.FDSize();
  UInt Nbands = 2*BW+1;

  MatInterFieldMat T_int( Nbands, Nx, Dim, Dim );
  VecInterFieldVec rhs_int( Nx, Dim );
  VecInterFieldVec x_int( Nx, Dim );

  T_int.setflag_inrealspace(false);
  rhs_int.setflag_inrealspace(false);
  x_int.setflag_inrealspace(false);

  T_int.smart_to_inter( T );
  rhs_int.smart_to_inter( rhs );

  // update both T and rhs using boundary conditions
  calc_BCs( T_int, rhs_int );

  // do a banded LU decomposition
  T_int.band_LU();

  // Use the LU matrix to solve the linear sytem
  T_int.band_solve( rhs_int, x_int );
  x_int.inter_to_smart( x );

} // }}}

void BCs_Const_Vel::calc_BCs( MatInterFieldMat & T, VecInterFieldVec & b )
{ // {{{

  // T: incoming intercalated block-banded matrix 
  //    (Nbands, Nx, Dim, Dim, Ny/Nz)
  //
  // [ B C        ...   A ]
  // [ A B C      ...     ]
  // [   A B C    ...     ]
  // [         ...        ]
  // [     ...    A B C   ]
  // [     ...      A B C ]
  // [ C   ...        A B ]
  //
  // b: rhs vector
  //    (Nx, Dim, Ny/Nz)

  int Nbands = T.Nrow();
  int Nx = T.Ncol();
  int BW = (Nbands-1)/2;
  int center = (Nbands-1)/2;
  RealType Dx = _CurrGrid->Dx();
  RealType Dx2 = Dx*Dx;

  InterFieldVec tmpvec( _Dim );

  // --- Modify the rows of T and b using FD equations ---
  // - T is banded, need to access (j-i+BW) like in LU

  // -------
  // x == 0
  // -------

  // Row 0, vel = BC on boundary
  int i = 0; // Row
  int max_j = std::min( i+BW+1, Nbands );

  for( int j=0; j<max_j; j++)
  {
    if (j-i+BW == center)
      T[j-i+BW][j].eye(1.);
    else
      T[j-i+BW][j].zero();
  }

  //b[i].zero();
  tmpvec.setflag_inrealspace(true);
  tmpvec.zero();
  for (int k=0; k<_Dim; k++) 
  {
    tmpvec[k] = _Vel0[k];
  }
  tmpvec.fft();
  b[i] = tmpvec;

  // -------
  // x == L
  // -------

  // Row Nx-1, vel = BC on boundary
 
  i = Nx-1;
  int min_j = std::max( i-BW, i-Nbands );

  for( int j=min_j; j<Nx; j++)
  {
    if (j-i+BW == center)
      T[j-i+BW][j].eye(1.);
    else
      T[j-i+BW][j].zero();
  }

  //b[i].zero();
  tmpvec.setflag_inrealspace(true);
  tmpvec.zero();
  for (int k=0; k<_Dim; k++) 
  {
    tmpvec[k] = _VelL[k];
  }
  tmpvec.fft();
  b[i] = tmpvec;

} // }}}

