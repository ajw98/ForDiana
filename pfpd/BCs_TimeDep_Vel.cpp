/***************************************************************
*
* Const Phi boundary condition solver
* 
* DRT -- Sat 15 Apr 2017
*
****************************************************************/

#include "BCs_TimeDep_Vel.h"

// ----------------- Constructor/Destructor -----------------

BCs_TimeDep_Vel :: BCs_TimeDep_Vel(std::string filename, int nicomp, 
                       Grid * gridarg ) : 
BCs_Base( nicomp, gridarg), 
_Dim( gridarg->Dim() )
{ // {{{

  // Read in Boundary conditions
  std::string tmp_str;

  if (jsonFile(filename.c_str()))
  {
    std::ifstream f1(filename.c_str());
    if ( f1.is_open() )
    {
      rapidjson::IStreamWrapper isw(f1);
      rapidjson::Document d;
      d.ParseStream(isw);

      _BCFlag = d["BCFlag"].GetInt();

      _Nt = d["Nt"].GetInt();

    }
    else
    {
      std::cout << "*** Error in BCs_TimeDep_Vel::BCs_TimeDep_Vel ***\n";
      std::cout << "*** Cannot open " << filename << "***" << std::endl;
      exit(1);
    }
  }
  else
  {//Legacy file format
    std::ifstream f1(filename.c_str());
    if ( f1.is_open() )
    {

      std::getline(f1, tmp_str, '#');
      _BCFlag = atoi(tmp_str.c_str());
      std::getline(f1, tmp_str, '\n');

      std::getline(f1, tmp_str, '#');
      _Nt = atoi(tmp_str.c_str());
      std::getline(f1, tmp_str, '\n');

    }
    else
    {
      std::cout << "*** Error in BCs_TimeDep_Vel::BCs_TimeDep_Vel ***\n";
      std::cout << "*** Cannot open " << filename << "***" << std::endl;
      exit(1);
    }
  }

  // read in time dependent BCs
  std::string timedep_file("params_BCs_vel_time.in");
  std::ifstream f2(timedep_file.c_str());
  if ( f2.is_open() )
  {

    _BC_times.resize(_Nt, 0.);
    std::getline(f2, tmp_str); // header
    for ( UInt i=0; i<_Nt; i++)
    {
      std::getline(f2, tmp_str);
      _BC_times[i] = atof(tmp_str.c_str());
    }
    std::getline(f2, tmp_str); // extra space

    _Vel0.resize(_Dim, Vec1dReal(_Nt, 0.));
    for ( UInt i=0; i<_Dim; i++)
    {
      std::getline(f2, tmp_str); // header
      for ( UInt j=0; j<_Nt; j++)
      {
        std::getline(f2, tmp_str);
        _Vel0[i][j] = atof(tmp_str.c_str());
      }
      std::getline(f2, tmp_str); // extra eol
    }

    _VelL.resize(_Dim, Vec1dReal(_Nt, 0.));
    for ( UInt i=0; i<_Dim; i++)
    {
      std::getline(f2, tmp_str); // header
      for ( UInt j=0; j<_Nt; j++)
      {
        std::getline(f2, tmp_str);
        _VelL[i][j] = atof(tmp_str.c_str());
      }
      std::getline(f2, tmp_str); // extra eol
    }

  }
  else
  {
    std::cout << "*** Error in BCs_TimeDep_Vel::BCs_TimeDep_Vel ***\n";
    std::cout << "*** Cannot open " << timedep_file << "***" << std::endl;
    exit(1);
  }

  std::cout << "\n" << std::scientific;
  std::cout << " * Vel Boundary Conditions Initiated\n" << std::scientific;
  std::cout << "   - BC Parameters:\n" << std::scientific;
  std::cout << "      BCFlag             = " << _BCFlag << std::endl;
  std::cout << "      Nt                 = " << _Nt << std::endl;
  for (int i=0; i<_Dim; ++i)
  {
    std::cout << "      Vel(x=0)[" << i << "][t=0]";
    std::cout <<"   = "<< _Vel0[i][0] << "\n";
  }
  for (int i=0; i<_Dim; ++i)
  {
    std::cout << "      Vel(x=L)[" << i << "][t=0]";
    std::cout <<"   = "<< _VelL[i][0] << "\n";
  }

} // }}}

BCs_TimeDep_Vel :: ~BCs_TimeDep_Vel()
{ // {{{
} // }}}

// ----------------- Constant BCs --------------------

void BCs_TimeDep_Vel::solve_BCs( SmartFieldOpMat & T, 
                              SmartFieldVec & rhs, 
                              SmartFieldVec & x,
                              RealType t )
{ // {{{

  UInt NIcomp = T.Nrow();
  UInt BW = T.BW();
  UInt Nx = T.FDSize();
  UInt Nbands = 2*BW+1;

  MatInterFieldMat T_int( Nbands, Nx, NIcomp, NIcomp );
  VecInterFieldVec rhs_int( Nx, NIcomp );
  VecInterFieldVec x_int( Nx, NIcomp );

  T_int.setflag_inrealspace(false);
  rhs_int.setflag_inrealspace(false);
  x_int.setflag_inrealspace(false);

  T_int.smart_to_inter( T );
  rhs_int.smart_to_inter( rhs );

  // update both T and rhs using boundary conditions
  calc_BCs( T_int, rhs_int, t );

  // do a banded LU decomposition
  T_int.band_LU();

  // Use the LU matrix to solve the linear sytem
  T_int.band_solve( rhs_int, x_int );
  x_int.inter_to_smart( x );

} // }}}

void BCs_TimeDep_Vel::calc_BCs( MatInterFieldMat & T, VecInterFieldVec & b, RealType t )
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
    tmpvec[k] = linterp(t, _BC_times, _Vel0[k]);
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
    tmpvec[k] = linterp(t, _BC_times, _VelL[k]);
  }
  tmpvec.fft();
  b[i] = tmpvec;

} // }}}

