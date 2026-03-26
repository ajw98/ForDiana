/***************************************************************
*
* VIPS boundary condition solver
* 
* JUG -- Wed 28 Nov 2018
* JUG -- Edited 27 Aug 2019
*
****************************************************************/

#include "BCs_VIPS.h"

// ----------------- Constructor/Destructor -----------------

BCs_VIPS :: BCs_VIPS(std::string filename, int nicomp, 
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
      
      _DdxAcc = d["ddxAcc"].GetInt();

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
  {//Legacy file format
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
      
      std::getline(f2, tmp_str, '#');
      _DdxAcc = atoi(tmp_str.c_str());
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
      
      _BC_xeq0.resize(_NBCs, Vec1dReal (_NIcomp, 0.));
      for ( UInt i=0; i<_NBCs; i++)
      {
        for ( UInt j=0; j<_NIcomp; j++)
        {
          std::getline(f2, tmp_str, '#');
          _BC_xeq0[i][j] = atof(tmp_str.c_str());
          std::getline(f2, tmp_str, '\n');
        }
      }

      _BC_xeqL.resize(_NBCs, Vec1dReal (_NIcomp, 0.));
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
  std::cout << "      BCFlag                      = " << _BCFlag << std::endl;
  std::cout << "      NBCs                        = " << _NBCs << std::endl;
  std::cout << "      GhostNodeFlag               = " << _GhostNodeFlag << std::endl;
  std::cout << "      Order accuracy of BD at Lx  = " << _DdxAcc << std::endl;
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
  check_input_Lx();  

  // Read in Constant thermodynamic parameters from file.
  _Ncomp = _NIcomp;
  _Ncomp += 1;
  
  if(jsonFile("params_EnergyModel.in"))
  {
    std::ifstream f3("params_EnergyModel.in");
    if ( f3.is_open() )
    {
      rapidjson::IStreamWrapper isw(f3);
      rapidjson::Document d;
      d.ParseStream(isw);
      
      // std::cout << "went into the  ifstream" << std::endl;
      _EnergyModelFlag = d["EnergyModelFlag"].GetInt();
      
      // std::cout << "EModelFlag=" << _EnergyModelFlag << std::endl;
      _Nr = d["Nr"].GetDouble();
      

      //std::cout << "Nr=" << _Nr << std::endl;
      _CReg = d["C_reg"].GetDouble();
      
      //std::cout << "CReg=" << _CReg << std::endl;

      _Delta = d["delta"].GetDouble();
      //std::cout << "Delta=" << _Delta << std::endl;

      _N.resize(_Ncomp, 0.);
      _alpha.resize(_Ncomp, 0.);
      const rapidjson::Value& N_in = d["N"];
      for( int i=0; i<_Ncomp; i++ )
      {
        
        _N[i] = N_in[i].GetDouble();
        

        _alpha[i] = _N[i]/_Nr;
      }
      // std::cout << "Np=" << _N[0] << std::endl;
      _Chi.resize(_Ncomp, Vec1dReal(_Ncomp, 0.));
      _ChiN.resize(_Ncomp, Vec1dReal(_Ncomp, 0.));
      rapidjson::Value& Chi_in = d["Chi"];
      for( int i=0; i<_Ncomp; i++ )
      {
        rapidjson::Value& Chi_in_2 = Chi_in[i];
        for( int j=i+1; j<_Ncomp; j++ )
        {
          
          _Chi[i][j] = Chi_in_2[j].GetDouble();
          
        
          _Chi[j][i] = _Chi[i][j];
          _ChiN[i][j] = _Chi[i][j]*_Nr;
          _ChiN[j][i] = _ChiN[i][j];
        }
      }
      //std::cout << "Chi12=" << _Chi[0][1] << std::endl;

      _Kappa.resize(_NIcomp, 0.);
      const rapidjson::Value& Kappa_in = d["Kappa"];
      for( int i=0; i<_NIcomp; i++ )
      {
        _Kappa[i] = Kappa_in[i].GetDouble();
      }

    // _Hmax.resize(_NIcomp, Vec1dFieldType(_NIcomp, FieldType(0.)));
      // std::cout << "Got to the end of reading parameters" << std::endl;
    }
    else
    { 
      std::cout << "Error opening params_EnergyModel.in" << std::endl;
      exit(1);
    }
    f3.close();
  }
  else
  { //Legacy file format
    std::ifstream f3("params_EnergyModel.in");
    if ( f3.is_open() )
    {
      // std::cout << "went into the  ifstream" << std::endl;
      std::getline(f3, tmp_str, '#');
      _EnergyModelFlag = atoi(tmp_str.c_str());
      std::getline(f3, tmp_str, '\n');
      
      // std::cout << "EModelFlag=" << _EnergyModelFlag << std::endl;
      std::getline(f3, tmp_str, '#');
      _Nr = atof(tmp_str.c_str());
      std::getline(f3, tmp_str, '\n');
      

      //std::cout << "Nr=" << _Nr << std::endl;
      std::getline(f3, tmp_str, '#');
      _CReg = atof(tmp_str.c_str());
      std::getline(f3, tmp_str, '\n');
      
      //std::cout << "CReg=" << _CReg << std::endl;
    
      std::getline(f3, tmp_str, '#');
      _Delta = atof(tmp_str.c_str());
      std::getline(f3, tmp_str, '\n');
      //std::cout << "Delta=" << _Delta << std::endl;

      _N.resize(_Ncomp, 0.);
      _alpha.resize(_Ncomp, 0.);
      for( int i=0; i<_Ncomp; i++ )
      {
        std::getline(f3, tmp_str, '#');
        _N[i] = atof(tmp_str.c_str());
        std::getline(f3, tmp_str, '\n');

        _alpha[i] = _N[i]/_Nr;
      }
      // std::cout << "Np=" << _N[0] << std::endl;
      _Chi.resize(_Ncomp, Vec1dReal(_Ncomp, 0.));
      _ChiN.resize(_Ncomp, Vec1dReal(_Ncomp, 0.));
      for( int i=0; i<_Ncomp; i++ )
      {
        for( int j=i+1; j<_Ncomp; j++ )
        {
          std::getline(f3, tmp_str, '#');
          _Chi[i][j] = atof(tmp_str.c_str());
          std::getline(f3, tmp_str, '\n');
        
          _Chi[j][i] = _Chi[i][j];
          _ChiN[i][j] = _Chi[i][j]*_Nr;
          _ChiN[j][i] = _ChiN[i][j];
        }
      }
      //std::cout << "Chi12=" << _Chi[0][1] << std::endl;

      _Kappa.resize(_NIcomp, 0.);
      for( int i=0; i < _NIcomp; i++ )
      {
        std::getline(f3, tmp_str, '#');
        _Kappa[i] = atof(tmp_str.c_str());
        std::getline(f3, tmp_str, '\n');
      }

    // _Hmax.resize(_NIcomp, Vec1dFieldType(_NIcomp, FieldType(0.)));
      // std::cout << "Got to the end of reading parameters" << std::endl;
    }
    else
    { 
      std::cout << "Error opening params_EnergyModel.in" << std::endl;
      exit(1);
    }
    f3.close();
  }
    

  // set up ChiN_ijM
  // e.g. ChiN_123 = _ChiN_12 - _ChiN_13 - _ChiN_23
  _ChiN_ijM.resize(_NIcomp, Vec1dReal(_NIcomp, 0.));
  for( int i=0; i<_NIcomp; i++ )
  {
    for( int j=0; j<_NIcomp; j++ )
    {
      if (i != j)
      {
        _ChiN_ijM[i][j] = _ChiN[i][j] - _ChiN[i][_Ncomp-1] - _ChiN[j][_Ncomp-1];
      }
      else
      {
        _ChiN_ijM[i][j] = 0.;
      }
    }
  }
  ////std::cout << "Just read in all thermo params." << std::endl; 
} // }}}

BCs_VIPS :: ~BCs_VIPS()
{ // {{{
} // }}}

// ----------------- Constant BCs --------------------

void BCs_VIPS::solve_BCs( SmartFieldOpMat & T, 
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
  x_int.smart_to_inter( x );
  
  
  // convert coefficients to VIPS
  convert_BCs2VIPS( x, Nx );
 
  // update both T and rhs using boundary conditions
  calc_BCs( T_int, rhs_int );

  // do a banded LU decomposition
  T_int.band_LU();

  // Use the LU matrix to solve the linear sytem
  T_int.band_solve( rhs_int, x_int );
  x_int.inter_to_smart( x );

} // }}}

void BCs_VIPS::calc_BCs( MatInterFieldMat & T, VecInterFieldVec & b )
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

void BCs_VIPS :: check_input_Lx()
{ // {{
  //check that the first BC for the nonsolvent is either zero or first order!
  if ( _Order_xeqL[0][1] > 1 )
  {
    std::cout << "***(Error)*** First BC for nonsolvent can only be set either as zero or first order" << std::endl;
    exit(1);
  }
  //tell user that 2nd BC for nonsolvent will be ignored.
  std::cout << "\n Second BC input for nonsolvent will be overwritten with VIPS formulation." << std::endl;
} // }}

void BCs_VIPS :: convert_BCs2VIPS ( SmartFieldVec & phi, int Nx )
{ // {{{ 
  //Converts the input BCs for phiN at x=L to appropriate VIPS BCs
  //i.e., changes _Order_xeqL[i][1] and _BC_xeqL[i][1] 
  // phi comes in k-space with index [component][x].

  //Step 1: Calculate averages at x=L
  //This is a shortcut, instead of having a locally varying f for each slice.
  //Note: the -1 offset is due to C++ counting from zero (verified index is correct).
  RealType phiP_xeqL = phi[0][Nx-1].getaverage().real();
  RealType phiN_xeqL = phi[1][Nx-1].getaverage().real();
  RealType phiN_xeqL_1 = phi[1][Nx-2].getaverage().real();
  RealType phiN_xeqL_2 = phi[1][Nx-3].getaverage().real();
  RealType Dx = _CurrGrid->Dx();
  RealType flux = 0; //calculate actual value in Step 3.
  
  //Step 2: Calculate g based from average phi's and thermodynamic model inputs.
  RealType g;
  calc_g( g, phiP_xeqL, phiN_xeqL );
  //g calculation is accurate.
  
  //Step 3: Calculate "flux"
  //We need dphi/dx for Step 4 but if user input is Dirichlet,
  //we need to calculate dphi/dx ourselves.
  if (_Order_xeqL[0][1] == 0) //0th order input
  {  
    if ( _DdxAcc == 6 ){
      //6th-order-accurate backward difference 
      //f' = ( 147phi(Lx) - 360phi(Lx-1) + 450phi(Lx-2) - 400phi(Lx-3)
      //       +225phi(Lx-4) -72phi(Lx-5) + 10phi(Lx-6)  )/(60Dx)
      RealType phiN_xeqL_3 = phi[1][Nx-4].getaverage().real();
      RealType phiN_xeqL_4 = phi[1][Nx-5].getaverage().real();
      RealType phiN_xeqL_5 = phi[1][Nx-6].getaverage().real();
      RealType phiN_xeqL_6 = phi[1][Nx-7].getaverage().real();
      flux = 147*phiN_xeqL;
      flux -= 360*phiN_xeqL_1;
      flux += 450*phiN_xeqL_2;
      flux -= 400*phiN_xeqL_3;
      flux += 225*phiN_xeqL_4;
      flux -= 72*phiN_xeqL_5;
      flux += 10*phiN_xeqL_6;
      flux /= (60*Dx);

    } else if ( _DdxAcc == 4 ){
      //4th-order-accurate backward difference 
      //f' = ( 25phi(Lx) - 48phi(Lx-1) + 36phi(Lx-2) - 16phi(Lx-3)
      //       +3phi(Lx-4) )/(12Dx)
      RealType phiN_xeqL_3 = phi[1][Nx-4].getaverage().real();
      RealType phiN_xeqL_4 = phi[1][Nx-5].getaverage().real();
      flux = 25*phiN_xeqL;
      flux -= 48*phiN_xeqL_1;
      flux += 36*phiN_xeqL_2;
      flux -= 16*phiN_xeqL_3;
      flux += 3*phiN_xeqL_4;
      flux /= (12*Dx);

    } else {
      //2nd-order-accurate backward difference---default 
      //f' = (3phi(Lx) - 4phi(Lx-1) + phi(Lx-2))/(2Dx)
      flux = 3*phiN_xeqL;
      flux -= 4*phiN_xeqL_1; 
      flux += phiN_xeqL_2; 
      flux /= (2*Dx);
    }
  }
  else //1st order input---use input by user.
  {
    flux = _BC_xeqL[0][1];
  } 

  //Step 4: d3f/dx = g*df/dx
  //Overwrites dummy inputs for 2nd BC for nonsolvent.
  _Order_xeqL[1][1] = 3;
  _BC_xeqL[1][1] = flux;
  _BC_xeqL[1][1] *= g;
  //std::cout << "_BC_xeqL[1][1]= " << _BC_xeqL[1][1] << std::endl; 
  
  //Step 5: check if polymer flux at x=L is actually zero.---only for debugging
  //print_jp (phiP_xeqL, phiN_xeqL,_BC_xeqL[0][0],_BC_xeqL[1][0], flux, _BC_xeqL[1][1]);

}// }}}

void BCs_VIPS :: print_jp ( RealType phiP0, RealType phiN0, RealType dp, RealType d3p,
		           RealType dn, RealType d3n )
{ // {{{ 
  // Verify that the polymer flux is actually zero.
  
  calc_Hess (phiP0, phiN0);

  RealType Mpp = phiP0*(1 - phiP0);
  RealType Mpn = -phiP0*phiN0;

  RealType jp = 0.;
  
  RealType t1 = Mpp;
  t1 *= _Hpp;
  t1 *= dp;  

  RealType t2 = Mpn;
  t2 *= _Hpn;
  t2 *= dp;  
  
  RealType t3 = Mpp;
  t3 *= _Kappa[0];
  t3 *= -d3p;  

  RealType t4 = Mpn;
  t4 *= 0.0;//Kpn
  t4 *= -d3p;  
  
  RealType t5 = Mpp;
  t5 *= _Hpn;
  t5 *= dn;  
  
  RealType t6 = Mpn;
  t6 *= _Hnn;
  t6 *= dn;  
  
  RealType t7 = Mpp;
  t7 *= 0.0;//Knp
  t7 *= -d3n;  
  
  RealType t8 = Mpn;
  t8 *= _Kappa[1];
  t8 *= -d3n; 

  jp = t1;
  jp += t2;
  jp += t3;
  jp += t4;
  jp += t5;
  jp += t6;
  jp += t7;
  jp += t8;

  std::cout << "***Polymer flux jp = " << jp << std::endl;
}// }}}


void BCs_VIPS :: calc_g ( RealType &g, RealType &phiP0, RealType &phiN0 )
{ // {{{ 
  //Calculate g from thermodynamic parameters and phi values
  
  //Step 1: Read in thermodynamic parameters
  calc_Hess( phiP0, phiN0 );
 
  //Step 2: Calculate g based from average phi's and thermodynamic model inputs.
  RealType num;
  
  g = phiP0;
  g *= phiN0;  
  g *= _Hnn;

  num = -phiP0;
  num += 1.;
  num *= -phiP0;
  num *= _Hpn;
  
  g += num;
  g /= _Kappa[1];
  g /= phiP0;
  g /= phiN0;

}// }}}

void BCs_VIPS :: calc_Hess( RealType & phiP0, RealType & phiN0 )
{ 
   //std::cout << "Just called calc_Hess" << std::endl;
   RealType phiM, tmp3;
   phiM = 1.;
   phiM -= phiN0;
   phiM -= phiP0;
  
   //term3 
   tmp3 = 1.;
   tmp3 /= phiM;
   tmp3 /= _alpha[_Ncomp-1];
  
   //Hpn
   _Hpn = tmp3;
   _Hpn += _ChiN_ijM[0][1];//verify  
  
   //Hpp 
   _Hpp = 1.;
   _Hpp /= phiP0;
   _Hpp /= _alpha[0];
   _Hpp += tmp3;
   _Hpp -= 2*_ChiN[0][2];

   //Hnn
   _Hnn = 1.;
   _Hnn /= phiN0;
   _Hnn /= _alpha[1];
   _Hnn += tmp3;
   _Hnn -= 2*_ChiN[1][2];
   //std::cout << "Finished calc_Hess calculation" << std::endl;
}
