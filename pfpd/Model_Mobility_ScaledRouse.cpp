/***************************************************************
*
* Class defining a Mobility model scaled by viscosity.
* Viscosity is sigmoidal model.
*
* JUG -- Thu, 28 Jun 2018
*
****************************************************************/

#include "Model_Mobility_ScaledRouse.h"

// ----------------- Constructor/Destructor -----------------
Model_Mobility_ScaledRouse :: Model_Mobility_ScaledRouse( int ncomp, int nicomp, Model_Visc_Base* visc_model ):
Model_Mobility_Base( ncomp, nicomp, visc_model )
{ // {{{
} // }}}

Model_Mobility_ScaledRouse :: ~Model_Mobility_ScaledRouse()
{ // {{{
} // }}}

// ----------------- Input/Output -----------------
void Model_Mobility_ScaledRouse :: set_params()
{ // {{{

  std::string tmp_str;

  if (jsonFile("params_MobilityModel.in"))
  {
    std::ifstream f1("params_MobilityModel.in");
    if ( f1.is_open() )
    {
      rapidjson::IStreamWrapper isw(f1);
      rapidjson::Document d;
      d.ParseStream(isw);
      _MobilityModelFlag = d["MobilityModelFlag"].GetInt();

      _Const.resize(_NIcomp, Vec1dReal(_NIcomp, 0.));
      rapidjson::Value& C_in = d["C"];
      for( int i=0; i<_NIcomp; i++ )
      {
        rapidjson::Value& C_in_2 = C_in[i];
        for( int j=0; j<_NIcomp; j++ )
        {
          _Const[i][j] = C_in_2[j].GetDouble();
         
        }
      }
    }
    else
    { 
      std::cout << "Error opening params_MobilityModel.in" << std::endl;
      exit(1);
    }
    f1.close();
  }
  else
  { // Legacy file format
    std::ifstream f3("params_MobilityModel.in");
    if ( f3.is_open() )
    {
      std::getline(f3, tmp_str, '#');
      _MobilityModelFlag = atoi(tmp_str.c_str());
      std::getline(f3, tmp_str, '\n');

      _Const.resize(_NIcomp, Vec1dReal(_NIcomp, 0.));
      for( int i=0; i<_NIcomp; i++ )
      {
        for( int j=0; j<_NIcomp; j++ )
        {
          std::getline(f3, tmp_str, '#');
          _Const[i][j] = atof(tmp_str.c_str());
          std::getline(f3, tmp_str, '\n');
        }
      }
    }
    else
    { 
      std::cout << "Error opening params_MobilityModel.in" << std::endl;
      exit(1);
    }
    f3.close();
  }
  
  std::cout << "   - Mobility Model Parameters:\n" << std::scientific;
  std::cout << "      MobilityModelFlag = " << _MobilityModelFlag << "\n";

} // }}}

// ----------------- Calculate Mobility -----------------
void Model_Mobility_ScaledRouse :: set_mobility( const SmartFieldVec &phi, SmartFieldMat &M )
{ // {{{

  //Calculate viscosity field 
  SmartField eta; 
  eta.setflag_inrealspace(true);
  eta.zero();

  _ViscModel->set_eta(phi, eta);  
  
  //Calculate mobility matrix
  M.setflag_inrealspace(phi.getflag_inrealspace()); // We need to do this because we start by setting elements of M with explicit constants
  for (int i=0; i<_NIcomp; i++)
  {
    for (int j=0; j<_NIcomp; j++)
    {
      // M_ii = phi_i * (1 - phi_i)
      if (i == j)
      {
        M(i,j) = 1.;
        M(i,j) -= phi[i];
        M(i,j) *= phi[i];
      }
      // M_ij = - phi_i * phi_j
      else // i != j
      {
        M(i,j).zero();
        M(i,j).accumulateproduct_inplace(phi[i], phi[j], -1.0);
      }
      M(i,j) /= eta; //mobility is scaled by viscosity
      M(i,j) *= _Const[i][j]; // also multiplied by prefactor
    }
  }

} // }}}

