/***************************************************************
*
* Class defining an exponential viscosity model similar to the
* form by Vogel-Fulcher-Tamman-Hesse (VFTH)
* 
* JUG --- 8 November 2018
*
****************************************************************/

#include "Model_Visc_VFTH.h"

// ----------------- Constructor/Destructor -----------------
Model_Visc_VFTH :: Model_Visc_VFTH( int ncomp, int nicomp ):
Model_Visc_Base( ncomp, nicomp )
{ // {{{
} // }}}

Model_Visc_VFTH :: ~Model_Visc_VFTH()
{ // {{{
} // }}}

// ----------------- Input/Output -----------------
void Model_Visc_VFTH :: set_params()
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
      
      _ViscCapFlag = d["etaCapFlag"].GetInt();

      _EtaR = d["eta_r"].GetDouble();

      _EtaS = d["eta_s"].GetDouble();
      
      _Gamma = d["gamma"].GetDouble();

      _PhiDagger = d["phiDagger"].GetDouble();
      
      _EtaCap = d["etaCap"].GetDouble();
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
      _ViscCapFlag = atoi(tmp_str.c_str());
      std::getline(f1, tmp_str, '\n');

      std::getline(f1, tmp_str, '#');
      _EtaR = atof(tmp_str.c_str());
      std::getline(f1, tmp_str, '\n');

      std::getline(f1, tmp_str, '#');
      _EtaS = atof(tmp_str.c_str());
      std::getline(f1, tmp_str, '\n');
      
      std::getline(f1, tmp_str, '#');
      _Gamma = atof(tmp_str.c_str());
      std::getline(f1, tmp_str, '\n');

      std::getline(f1, tmp_str, '#');
      _PhiDagger = atof(tmp_str.c_str());
      std::getline(f1, tmp_str, '\n');
      
      std::getline(f1, tmp_str, '#');
      _EtaCap = atof(tmp_str.c_str());
      std::getline(f1, tmp_str, '\n');
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
  std::cout << "      ViscCapFlag   = " << _ViscCapFlag << "\n";
  std::cout << "      EtaR            = " << _EtaR << "\n";
  std::cout << "      EtaS            = " << _EtaS << "\n";
  std::cout << "      Gamma           = " << _Gamma << "\n";
  std::cout << "      PhiDagger       = " << _PhiDagger << "\n";
  std::cout << "      EtaCap          = " << _EtaCap << "\n";

} // }}}

void Model_Visc_VFTH :: set_eta( const SmartFieldVec &phi, SmartField &eta )
{ // {{{
  // (in) phi: incoming concentrations (but we're only using phi[0], the polymer).
  // (out) eta: viscosity 
  //
  // eta = eta_s/eta_r *exp{-gamma/phiDagger *(1+phiDagger/(phi_p - phiDagger))}}
  
  // cap polymer field below glass concentration threshold.
  SmartField tmp = phi[0];//operate on a copy, in case capping is required! 
  if ( _ViscCapFlag > 0 )
  {
    RealType phiCap = phi_of_Eta ( _EtaCap );

    capField ( tmp, phiCap );
    if ( tmp.maxsigned().real() > phiCap )  
    {
      std::cout << "* Warning: capField called, but FieldMax = " << tmp.maxsigned().real();    
      std::cout << " while FieldCap = " << phiCap << std::endl; 
    }    
  }    

  //eta = 1/(phi_p - phiDagger)
  eta = tmp;
  eta += -_PhiDagger;
  eta.log(eta);
  eta.exponentiate(eta, -1.);

  //eta = exp{-gamma/phiDagger *(1+phiDagger/(phi_p - phiDagger))}
  eta *= _PhiDagger;
  eta += 1.;
  eta.exponentiate(eta, -_Gamma/_PhiDagger);  
  
  // eta = eta_s/eta_r *exp{-gamma/phiDagger *(1+phiDagger/(phi_p - phiDagger))}}
  eta *= _EtaS/_EtaR;
  
  // Print out maximum viscosity, for testing. Commented out for normal use.
  // std::cout << "EtaMax = " << eta.maxsigned().real() << std::endl;
 
  // Check for violators. 
  if ( eta.maxsigned().real() > _EtaCap )  
  {
    std::cout << "* Warning: tried to cap Eta, but EtaMax = " << eta.maxsigned().real();    
    std::cout << " while EtaCap = " << _EtaCap << std::endl;
    std::cout << " Stopping simulation here---check phi fields for divergence!" << std::endl;
    exit(1); 
  }    
  
} // }}}

RealType Model_Visc_VFTH :: phi_of_Eta ( const RealType eta )
{ // {{{
  // Calculate the corrsponding phi to the given eta (inverse function of VFTH).
  // phi_of_Eta = phiDagger - gamma/log(eta/A);
  // A = etaS/etaR*exp(-gamma/phiDagger);

  RealType A = exp(-_Gamma/_PhiDagger);
  A *= _EtaS/_EtaR;
  
  RealType B = log(eta/A);
  RealType phiCap  = _Gamma/B;
  phiCap *= -1.;
  phiCap += _PhiDagger;
  return phiCap; 

} // }}}{

void Model_Visc_VFTH :: capField ( SmartField &F, const RealType k )
{ // {{{
  // Caps the field to a maximum value of k. 
   
  F -= k;
  F.zeropos();
  F += k;
  
} // }}}{
