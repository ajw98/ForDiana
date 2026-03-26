/***************************************************************
*
* Class defining a Flory Huggins free energy model
* 
* DRT -- Wed, 18 Aug 2015
*
****************************************************************/

#include "Model_Mobility_Const.h"

// ----------------- Constructor/Destructor -----------------
Model_Mobility_Const :: Model_Mobility_Const( int ncomp, int nicomp, Model_Visc_Base* visc_model ):
Model_Mobility_Base( ncomp, nicomp, visc_model )
{ // {{{
} // }}}

Model_Mobility_Const :: ~Model_Mobility_Const()
{ // {{{
} // }}}

// ----------------- Input/Output -----------------
void Model_Mobility_Const :: set_params()
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

      _Mobility.resize(_NIcomp, Vec1dReal(_NIcomp, 0.));
      rapidjson::Value& mob_in = d["Mob"];
      for( int i=0; i<_NIcomp; i++ )
      {
        rapidjson::Value& mob_in_2 = mob_in[i];
        for( int j=0; j<_NIcomp; j++ )
        {
          _Mobility[i][j] = mob_in_2[j].GetDouble();
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
  {//Legacy file format
    
    std::ifstream f1("params_MobilityModel.in");
    if ( f1.is_open() )
    {
      std::getline(f1, tmp_str, '#');
      _MobilityModelFlag = atoi(tmp_str.c_str());
      std::getline(f1, tmp_str, '\n');

      _Mobility.resize(_NIcomp, Vec1dReal(_NIcomp, 0.));
      for( int i=0; i<_NIcomp; i++ )
      {
        for( int j=0; j<_NIcomp; j++ )
        {
          std::getline(f1, tmp_str, '#');
          _Mobility[i][j] = atof(tmp_str.c_str());
          std::getline(f1, tmp_str, '\n');
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
    

  std::cout << "   - Mobility Model Parameters:\n" << std::scientific;
  std::cout << "      MobilityModelFlag = " << _MobilityModelFlag << "\n";
  for( int i=0; i<_NIcomp; i++ )
  {
    for( int j=0; j<_NIcomp; j++ )
    {
      std::cout << "      M" << i << j << "             = ";
      std::cout << _Mobility[i][j] << "\n";
    }
  }

} // }}}

void Model_Mobility_Const :: set_mobility( const SmartFieldVec &phi, SmartFieldMat &M )
{// {{{

  M.setflag_inrealspace(phi.getflag_inrealspace()); // We need to do this because we start by setting elements of M with explicit constants
  for( int i=0; i<_NIcomp; i++ )
  {
    for( int j=0; j<_NIcomp; j++ )
    {
      M(i,j) = _Mobility[i][j];
    }
  }

} // }}}

