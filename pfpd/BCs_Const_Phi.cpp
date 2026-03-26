/***************************************************************
*
* Const Phi boundary condition solver
* 
* DRT -- Sat 15 Apr 2017
*
****************************************************************/

#include "BCs_Const_Phi.h"

// ----------------- Constructor/Destructor -----------------

BCs_Const_Phi :: BCs_Const_Phi(std::string filename, int nicomp, 
                       Grid * gridarg ) : 
BCs_Base( nicomp, gridarg )
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

      _NBCs = d["NBCs"].GetInt();

      _GhostNodeFlag = d["GhostNodeFlag"].GetInt();

      _Order_xeq0.resize(_NBCs, Vec1dInt(_NIcomp, 0));
      rapidjson::Value& xeq0_in = d["Order_xeq0"];
      for ( UInt i=0; i<_NBCs; i++)
      {
        rapidjson::Value& xeq0_in_2 = xeq0_in[i];
        for ( UInt j=0; j<_NIcomp; j++)
        {
          _Order_xeq0[i][j] = xeq0_in_2[j].GetInt();
        }
      }

      _Order_xeqL.resize(_NBCs, Vec1dInt(_NIcomp, 0));
      rapidjson::Value& xeqL_in = d["Order_xeqL"];
      for ( UInt i=0; i<_NBCs; i++)
      {
        rapidjson::Value& xeqL_in_2 = xeqL_in[i];
        for ( UInt j=0; j<_NIcomp; j++)
        {
          _Order_xeqL[i][j] = xeqL_in_2[j].GetInt();
        }
      }

      _BC_xeq0.resize(_NBCs, Vec1dReal(_NIcomp, 0.));
      rapidjson::Value& BC_xeq0_in = d["BC_xeq0"];
      for ( UInt i=0; i<_NBCs; i++)
      {
        rapidjson::Value& BC_xeq0_in_2 = BC_xeq0_in[i];
        for ( UInt j=0; j<_NIcomp; j++)
        {
          _BC_xeq0[i][j] = BC_xeq0_in_2[j].GetDouble();
        }
      }

      _BC_xeqL.resize(_NBCs, Vec1dReal(_NIcomp, 0.));
      rapidjson::Value& BC_xeqL_in = d["BC_xeqL"];
      for ( UInt i=0; i<_NBCs; i++)
      {
        rapidjson::Value& BC_xeqL_in_2 = BC_xeqL_in[i];
        for ( UInt j=0; j<_NIcomp; j++)
        {
          _BC_xeqL[i][j] = BC_xeqL_in_2[j].GetDouble();
        }
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

      std::getline(f2, tmp_str, '#');
      _NBCs = atoi(tmp_str.c_str());
      std::getline(f2, tmp_str, '\n');

      std::getline(f2, tmp_str, '#');
      _GhostNodeFlag = atoi(tmp_str.c_str());
      std::getline(f2, tmp_str, '\n');

      _Order_xeq0.resize(_NBCs, Vec1dInt(_NIcomp, 0));
      for ( UInt i=0; i<_NBCs; i++)
      {
        for ( UInt j=0; j<_NIcomp; j++)
        {
          std::getline(f2, tmp_str, '#');
          _Order_xeq0[i][j] = atoi(tmp_str.c_str());
          std::getline(f2, tmp_str, '\n');
        }
      }

      _Order_xeqL.resize(_NBCs, Vec1dInt(_NIcomp, 0));
      for ( UInt i=0; i<_NBCs; i++)
      {
        for ( UInt j=0; j<_NIcomp; j++)
        {
          std::getline(f2, tmp_str, '#');
          _Order_xeqL[i][j] = atoi(tmp_str.c_str());
          std::getline(f2, tmp_str, '\n');
        }
      }

      _BC_xeq0.resize(_NBCs, Vec1dReal(_NIcomp, 0.));
      for ( UInt i=0; i<_NBCs; i++)
      {
        for ( UInt j=0; j<_NIcomp; j++)
        {
          std::getline(f2, tmp_str, '#');
          _BC_xeq0[i][j] = atof(tmp_str.c_str());
          std::getline(f2, tmp_str, '\n');
        }
      }

      _BC_xeqL.resize(_NBCs, Vec1dReal(_NIcomp, 0.));
      for ( UInt i=0; i<_NBCs; i++)
      {
        for ( UInt j=0; j<_NIcomp; j++)
        {
          std::getline(f2, tmp_str, '#');
          _BC_xeqL[i][j] = atof(tmp_str.c_str());
          std::getline(f2, tmp_str, '\n');
        }
      }

    }
    else
    {
      std::cout << "cannot open " << filename << std::endl;
      exit(1);
    }
  }

  std::cout << "\n" << std::scientific;
  std::cout << " * Phi Boundary Conditions Initiated\n" << std::scientific;
  std::cout << "   - BC Parameters:\n" << std::scientific;
  std::cout << "      BCFlag             = " << _BCFlag << std::endl;
  std::cout << "      NBCs               = " << _NBCs << std::endl;
  std::cout << "      GhostNodeFlag      = " << _GhostNodeFlag << std::endl;
  for (int i=0; i<_NBCs; ++i)
  {
    for (int j=0; j<_NIcomp; ++j)
    {
      std::cout << "      Order_xeq0[" << i << "]" << "[" << j <<"]";
      std::cout <<"   = "<< _Order_xeq0[i][j] << "\n";
    }
  }
  for (int i=0; i<_NBCs; ++i)
  {
    for (int j=0; j<_NIcomp; ++j)
    {
      std::cout << "      Order_xeqL[" << i << "]" << "[" << j <<"]";
      std::cout <<"   = "<< _Order_xeqL[i][j] << "\n";
    }
  }
  for (int i=0; i<_NBCs; ++i)
  {
    for (int j=0; j<_NIcomp; ++j)
    {
      std::cout << "      BC_xeq0[" << i << "]" << "[" << j <<"]";
      std::cout <<"      = "<< _BC_xeq0[i][j] << "\n";
    }
  }
  for (int i=0; i<_NBCs; ++i)
  {
    for (int j=0; j<_NIcomp; ++j)
    {
      std::cout << "      BC_xeqL[" << i << "]" << "[" << j <<"]";
      std::cout <<"      = "<< _BC_xeqL[i][j] << "\n";
    }
  }

} // }}}

BCs_Const_Phi :: ~BCs_Const_Phi()
{ // {{{
} // }}}

// ----------------- Constant BCs --------------------

void BCs_Const_Phi::solve_BCs( SmartFieldOpMat & T, 
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
  calc_BCs( T_int, rhs_int );

  // do a banded LU decomposition
  T_int.band_LU();

  // Use the LU matrix to solve the linear sytem
  T_int.band_solve( rhs_int, x_int );
  x_int.inter_to_smart( x );

} // }}}

void BCs_Const_Phi::calc_BCs( MatInterFieldMat & T, VecInterFieldVec & b )
{ // {{{

  // T: incoming intercalated block-banded matrix 
  //    (Nbands, Nx, NIcomp, NIcomp, Ny/Nz)
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
  //    (Nx, NIcomp, Ny/Nz)

  int Nbands = T.Nrow();
  int Nx = T.Ncol();
  int BW = (Nbands-1)/2;
  int center = (Nbands-1)/2;
  RealType Dx = _CurrGrid->Dx();
  RealType Dx2 = Dx*Dx;
  RealType Dx3 = Dx2*Dx;

  InterFieldVec tmpvec( _NIcomp );

  // --- set finite difference equations for each component ---

  Vec3dReal FD_xeq0( _NBCs, Vec2dReal(Nbands, Vec1dReal(_NIcomp, 0.)) );
  Vec3dReal FD_xeqL( _NBCs, Vec2dReal(Nbands, Vec1dReal(_NIcomp, 0.)) );

  for (int i=0; i<_NBCs; i++)
  {
    for (int k=0; k<_NIcomp; k++)
    {
 
      // ------------- 
      //     x = 0
      // ------------- 

      // f(x=0) = _BC
      if (_Order_xeq0[i][k] == 0)
      {
        for (int j=0; j<Nbands; j++)
        {
          if (j == center-i)
            FD_xeq0[i][j][k] = 1.;
          else
            FD_xeq0[i][j][k] = 0.;
        }
      }

      // f'(x=0) = _BC
      if (_Order_xeq0[i][k] == 1)
      {
        for (int j=0; j<Nbands; j++)
        {
          if (j == center-i)
            FD_xeq0[i][j][k] = -3./(2*Dx);
          else if (j == center+1-i)
            FD_xeq0[i][j][k] = 4./(2*Dx);
          else if (j == center+2-i)
            FD_xeq0[i][j][k] = -1./(2*Dx);
          else
            FD_xeq0[i][j][k] = 0.;
        }
      }

      // f''(x=0) = _BC
      if (_Order_xeq0[i][k] == 2 and _GhostNodeFlag == 0)
      {
        for (int j=0; j<Nbands; j++)
        {
          if (j == center-i)
            FD_xeq0[i][j][k] = 1./(Dx2);
          else if (j == center+1-i)
            FD_xeq0[i][j][k] = -2./(Dx2);
          else if (j == center+2-i)
            FD_xeq0[i][j][k] = 1./(Dx2);
          else
            FD_xeq0[i][j][k] = 0.;
        }
      }

      if (_Order_xeq0[i][k] == 2 and _GhostNodeFlag == 1)
      {
        for (int j=0; j<Nbands; j++)
        {
          if (j == center-i)
            FD_xeq0[i][j][k] = 2./(Dx2);
          else if (j == center+1-i)
            FD_xeq0[i][j][k] = -5./(Dx2);
          else if (j == center+2-i)
            FD_xeq0[i][j][k] = 4./(Dx2);
          else if (j == center+3-i)
            FD_xeq0[i][j][k] = -1./(Dx2);
          else
            FD_xeq0[i][j][k] = 0.;
        }
      }

      // f'''(x=0) = _BC
      if (_Order_xeq0[i][k] == 3)
      {
        for (int j=0; j<Nbands; j++)
        {
          if (j == center-i)
            FD_xeq0[i][j][k] = -1./(Dx3);
          else if (j == center+1-i)
            FD_xeq0[i][j][k] = 3./(Dx3);
          else if (j == center+2-i)
            FD_xeq0[i][j][k] = -3./(Dx3);
          else if (j == center+3-i)
            FD_xeq0[i][j][k] = 1./(Dx3);
          else
            FD_xeq0[i][j][k] = 0.;
        }
      }

      // ------------- 
      //     x = L
      // ------------- 

      // f(x=L) = _BC
      if (_Order_xeqL[i][k] == 0)
      {
        for (int j=0; j<Nbands; j++)
        {
          if (j == center+i)
            FD_xeqL[i][j][k] = 1.;
          else
            FD_xeqL[i][j][k] = 0.;
        }
      }

      // f'(x=L) = _BC
      if (_Order_xeqL[i][k] == 1)
      {
        for (int j=0; j<Nbands; j++)
        {
          if (j == center-2+i)
            FD_xeqL[i][j][k] = 1./(2*Dx);
          else if (j == center-1+i)
            FD_xeqL[i][j][k] = -4./(2*Dx);
          else if (j == center+i)
            FD_xeqL[i][j][k] = 3./(2*Dx);
          else
            FD_xeqL[i][j][k] = 0.;
        }
      }

      // f''(x=L) = _BC
      if (_Order_xeqL[i][k] == 2 and _GhostNodeFlag == 0)
      {
        for (int j=0; j<Nbands; j++)
        {
          if (j == center-2+i)
            FD_xeqL[i][j][k] = 1./Dx2;
          else if (j == center-1+i)
            FD_xeqL[i][j][k] = -2./Dx2;
          else if (j == center+i)
            FD_xeqL[i][j][k] = 1./Dx2;
          else
            FD_xeqL[i][j][k] = 0.;
        }
      }

      if (_Order_xeqL[i][k] == 2 and _GhostNodeFlag == 1)
      {
        for (int j=0; j<Nbands; j++)
        {
          if (j == center-2+i)
            FD_xeqL[i][j][k] = -1./Dx2;
          if (j == center-2+i)
            FD_xeqL[i][j][k] = 4./Dx2;
          else if (j == center-1+i)
            FD_xeqL[i][j][k] = -5./Dx2;
          else if (j == center+i)
            FD_xeqL[i][j][k] = 2./Dx2;
          else
            FD_xeqL[i][j][k] = 0.;
        }
      }


      // f'''(x=L) = _BC
      // (Warning: only first order derivative)
      if (_Order_xeq0[i][k] == 3)
      {
        for (int j=0; j<Nbands; j++)
        {
          if (j == center-3+i)
            FD_xeqL[i][j][k] = -1./(Dx3);
          else if (j == center-2+i)
            FD_xeqL[i][j][k] = 3./(Dx3);
          else if (j == center-1+i)
            FD_xeqL[i][j][k] = -3./(Dx3);
          else if (j == center+i)
            FD_xeqL[i][j][k] = 1./(Dx3);
          else
            FD_xeqL[i][j][k] = 0.;
        }
      }

    }
  }

  // --- Modify the rows of T and b using FD equations from above ---

  // x == 0
  // - loop over rows for # of BCs
  // - T is banded, need to access (j-i+BW) like in LU

  int max_j = 0;

  for (int i=0; i<_NBCs; i++) // row where we will install the equation
  {

    max_j = std::min( i+BW+1, Nbands );

    for( int j=0; j<max_j; j++)
    {
      
      T[j-i+BW][j].zero();
      for( int k=0; k<_NIcomp; k++)
      {
        T[j-i+BW][j][k][k] = FD_xeq0[i][j-i+BW][k];
      }

    }

    tmpvec.setflag_inrealspace(true);
    tmpvec.zero();
    for (int k=0; k<_NIcomp; k++) 
    {
      tmpvec[k] = _BC_xeq0[i][k];
    }
    tmpvec.fft();
    b[i] = tmpvec;

  }

  // x == L
  // - loop over rows for # of BCs
  // - T is banded, need to access (j-i+BW) like in LU

  int min_j = 0;

  for (int i=Nx-_NBCs; i<Nx; i++)
  {

    min_j = std::max( i-BW, i-Nbands );

    for( int j=min_j; j<Nx; j++)
    {

      T[j-i+BW][j].zero();
      for( int k=0; k<_NIcomp; k++)
      {
        T[j-i+BW][j][k][k] = FD_xeqL[Nx-_NBCs-i+1][j-i+BW][k];
      }

    }

    tmpvec.setflag_inrealspace(true);
    tmpvec.zero();
    for (int k=0; k<_NIcomp; k++) 
    {
      tmpvec[k] = _BC_xeqL[Nx-_NBCs-i+1][k];
    }
    tmpvec.fft();
    b[i] = tmpvec;

  }

} // }}}

