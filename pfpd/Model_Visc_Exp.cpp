/***************************************************************
*
* Class defining a sigmoidal viscosity model
* 
* DRT -- Wed, 18 Aug 2015
*
****************************************************************/

#include "Model_Visc_Exp.h"

// ----------------- Constructor/Destructor -----------------
Model_Visc_Exp :: Model_Visc_Exp( int ncomp, int nicomp ):
Model_Visc_Base( ncomp, nicomp )
{ // {{{
} // }}}

Model_Visc_Exp :: ~Model_Visc_Exp()
{ // {{{
} // }}}

// ----------------- Input/Output -----------------
void Model_Visc_Exp :: set_params()
{ // {{{

  if(jsonFile("params_ViscModel.in"))
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

      _PhiT.resize(_NIcomp, 0.);
      rapidjson::Value& phiT_in = d["phi_T"];
      for( int i=0; i<_NIcomp; i++ )
      {
        _PhiT[i] = phiT_in[i].GetDouble();
      }

      _W.resize(_NIcomp, 0.);
      rapidjson::Value& w_in = d["w"];
      for( int i=0; i<_NIcomp; i++ )
      {
        _W[i] = w_in[i].GetDouble();
      }

      _ExpOn.resize(_NIcomp, 0.);
      rapidjson::Value& expon_in = d["expon"];
      for( int i=0; i<_NIcomp; i++ )
      {
        _ExpOn[i] = expon_in[i].GetInt();
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

      _PhiT.resize(_NIcomp, 0.);
      for( int i=0; i<_NIcomp; i++ )
      {
        std::getline(f1, tmp_str, '#');
        _PhiT[i] = atof(tmp_str.c_str());
        std::getline(f1, tmp_str, '\n');
      }
     
      _W.resize(_NIcomp, 0.);
      for( int i=0; i<_NIcomp; i++ )
      {
        std::getline(f1, tmp_str, '#');
        _W[i] = atof(tmp_str.c_str());
        std::getline(f1, tmp_str, '\n');
      }

      _ExpOn.resize(_NIcomp, 0.);
      for( int i=0; i<_NIcomp; i++ )
      {
        std::getline(f1, tmp_str, '#');
        _ExpOn[i] = atoi(tmp_str.c_str());
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

  for( int i=0; i<_NIcomp; i++ )
  {
    std::cout << "     PhiT" << i << "            = ";
    std::cout << _PhiT[i] << "\n";
  }
  
  for( int i=0; i<_NIcomp; i++ )
  {
    std::cout << "        W" << i << "            = ";
    std::cout << _W[i] << "\n";
  }

  for( int i=0; i<_NIcomp; i++ )
  {
    std::cout << "    ExpOn" << i << "            = ";
    std::cout << _ExpOn[i] << "\n";
  }

} // }}}

void Model_Visc_Exp :: set_eta( const SmartFieldVec &phi, SmartField &eta )
{ // {{{

  SmartField tmp;

  // eta = 1. + sum_i f(phi_i)*(eta_iMR)
  //   where eta_iMR = (eta_i - eta_M)/eta_R
  //   with M = Ncomp
  eta = 1.;
  for (int i=0; i<_NIcomp; i++)
  {
    if (_ExpOn[i])
    {
      tmp = phi[i];
      f_exp(tmp, _PhiT[i], _W[i]); // tmp -> f_sigmoid(tmp)
      eta.xpby_inplace(tmp, (_Eta[i]-_Eta[_Ncomp-1])/_EtaR);
    } else {
      eta.xpby_inplace(phi[i], (_Eta[i]-_Eta[_Ncomp-1])/_EtaR);
    }
//    tmp *= (_Eta[i]-_Eta[_Ncomp-1])/_EtaR;
//    eta += tmp;
  }
  
} // }}}

void Model_Visc_Exp :: f_exp( SmartField &phi, RealType phiT, RealType w)
{ // {{{

  // (in) phi: incoming concentration
  // (in) phiT: threshold
  // (in) w: characteristic width of transition in phi
  // (out) phi: sigmoida
  //
  // f(phi) = exp((phi-phiT)/w)

  // phi will be replaced, so work in-place without tmp field
  phi += -phiT;
  phi.exponentiate(phi, 1./w);

} // }}}

