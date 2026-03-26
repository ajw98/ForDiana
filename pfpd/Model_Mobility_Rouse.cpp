/***************************************************************
*
* Class defining a Flory Huggins free energy model
* 
* DRT -- Wed, 18 Aug 2015
*
****************************************************************/

#include "Model_Mobility_Rouse.h"

// ----------------- Constructor/Destructor -----------------
Model_Mobility_Rouse :: Model_Mobility_Rouse( int ncomp, int nicomp, Model_Visc_Base* visc_model ):
Model_Mobility_Base( ncomp, nicomp, visc_model )
{ // {{{
} // }}}

Model_Mobility_Rouse :: ~Model_Mobility_Rouse()
{ // {{{
} // }}}

// ----------------- Input/Output -----------------
void Model_Mobility_Rouse :: set_params()
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

    }
    else
    { 
      std::cout << "Error opening params_MobilityModel.in" << std::endl;
      exit(1);
    }
    f1.close();
  }
  else
  {//Legacy file format
    std::ifstream f3("params_MobilityModel.in");
    if ( f3.is_open() )
    {
      std::getline(f3, tmp_str, '#');
      _MobilityModelFlag = atoi(tmp_str.c_str());
      std::getline(f3, tmp_str, '\n');

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

void Model_Mobility_Rouse :: set_mobility( const SmartFieldVec &phi, SmartFieldMat &M )
{ // {{{

  M.setflag_inrealspace(phi.getflag_inrealspace()); // We need to do this because we start by setting elements of M with explicit constants
  for (int i=0; i<_NIcomp; i++)
  {
    // Diagonals
    // M_ii = C_ii * phi_i * (1 - phi_i)
    M(i,i) = 1.;
    M(i,i) -= phi[i];
    M(i,i) *= phi[i];

    // Off diagonals
    // M_ij = - C_ij * phi_i * phi_j
    for (int j=i+1; j<_NIcomp; j++)
    {
        M(i,j).zero();
        M(i,j).accumulateproduct_inplace(phi[i], phi[j], -1.0);
        M(j,i) = M(i,j); // Symmetric matrix
    }
  }

} // }}}

