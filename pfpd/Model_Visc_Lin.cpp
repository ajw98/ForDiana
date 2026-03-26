/***************************************************************
*
* Class defining a linear viscosity model
* 
* DRT -- Wed, 18 Aug 2015
*
****************************************************************/

#include "Model_Visc_Lin.h"

// ----------------- Constructor/Destructor -----------------
Model_Visc_Lin :: Model_Visc_Lin( int ncomp, int nicomp ):
Model_Visc_Base( ncomp, nicomp )
{ // {{{
} // }}}

Model_Visc_Lin :: ~Model_Visc_Lin()
{ // {{{
} // }}}

// ----------------- Input/Output -----------------
void Model_Visc_Lin :: set_params()
{ // {{{
  
  if (jsonFile("params_ViscModel.in"))
  {
    std::ifstream f1("params_ViscModel.in");
    if ( f1.is_open() )
    {
      rapidjson::IStreamWrapper isw(f1);
      rapidjson::Document d;
      d.ParseStream(isw);

      _ViscModelFlag = d["ViscModelFlag"].GetInt();

      _EtaR = d["eta_r"].GetDouble();

      _Eta.resize(_Ncomp, 0.);
      rapidjson::Value& eta = d["eta"];
      for( int i=0; i<_Ncomp; i++ )
      {
        _Eta[i] = eta[i].GetDouble();
      }
    }
    else
    {
      std::cout << "Error opening params_ViscModel.in" << std::endl;
      exit(0);
    }
    f1.close();
  }
  else
  {//Legacy file format
    std::string tmp_str;
    std::ifstream f1("params_ViscModel.in");
    if ( f1.is_open() )
    {

      std::getline(f1, tmp_str, '#');
      _ViscModelFlag = atoi(tmp_str.c_str());
      std::getline(f1, tmp_str, '\n');

      std::getline(f1, tmp_str, '#');
      _EtaR = atof(tmp_str.c_str());
      std::getline(f1, tmp_str, '\n');

      _Eta.resize(_Ncomp, 0.);
      for( int i=0; i<_Ncomp; i++ )
      {
        std::getline(f1, tmp_str, '#');
        _Eta[i] = atof(tmp_str.c_str());
        std::getline(f1, tmp_str, '\n');
      }

    }
    else
    {
      std::cout << "Error opening params_ViscModel.in" << std::endl;
      exit(0);
    }
    f1.close();
  }

  std::cout << "   - Viscosity Model Parameters:\n" << std::scientific;
  std::cout << "      ViscModelFlag   = " << _ViscModelFlag << "\n";
  std::cout << "      EtaR            = " << _EtaR << "\n";
  for( int i=0; i<_Ncomp; i++ )
  {
    std::cout << "      Eta" << i << "            = ";
    std::cout << _Eta[i] << "\n";
  }

} // }}}

void Model_Visc_Lin :: set_eta( const SmartFieldVec &phi, SmartField &eta )
{ // {{{

//  SmartField tmp;

  // e.g. for ternary system:
  // eta = 1. + phi1*(eta_1-eta_3)/(eta_r) + phi2*(eta_2-eta_3)/(eta_r)
//  eta = 1.;
//  for( int i=0; i<_NIcomp; i++ )
//  {
//    tmp = phi[i];
//    tmp *= (_Eta[i]-_Eta[_Ncomp-1])/_EtaR;
//    eta += tmp;
//  }

  eta = 1.;
  for(int i=0; i<_NIcomp; ++i)
    eta.xpby_inplace(phi[i], (_Eta[i]-_Eta[_Ncomp-1])/_EtaR);

} // }}}

