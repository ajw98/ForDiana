/***************************************************************
*
* Class defining a Linear free energy model
* 
* DRT -- Wed, 18 Aug 2015
*
****************************************************************/
//
// mu_i = sum_j ( H_ij phi_j ) - K_i del^2 phi_i
//

#include "Model_Energy_Linear.h"

// ----------------- Constructor/Destructor -----------------
Model_Energy_Linear::Model_Energy_Linear( Operators_Base & OpObj, int ncomp, int nicomp ): 
Model_Energy_Base( OpObj, ncomp, nicomp )
{ // {{{
} // }}}

Model_Energy_Linear::~Model_Energy_Linear()
{ // {{{
} // }}}

// ----------------- Input/Output -----------------
void Model_Energy_Linear::set_params()
{ // {{{

  std::string tmp_str;
  
  if (jsonFile("params_EnergyModel.in"))
  {
    std::ifstream f2("params_EnergyModel.in");
    if ( f2.is_open() )
    {
      rapidjson::IStreamWrapper isw(f2);
      rapidjson::Document d;
      d.ParseStream(isw);
      
      _EnergyModelFlag = d["EnergyModelFlag"].GetInt();
      
      rapidjson::Value& Hessin = d["Hess"];
      _Hess.resize(_NIcomp, Vec1dFieldType(_NIcomp, 0.));
      for( int i=0; i<_NIcomp; i++ )
      {
        rapidjson::Value& dataIn = Hessin[i];
        for( int j=0; j<_NIcomp; j++ )
        {
          _Hess[i][j] = FieldType(dataIn[j].GetDouble(), 0.);
        }
      }

      _Kappa.resize(_NIcomp, 0.);
      const rapidjson::Value& KappaIn = d["Kappa"];
      for( int i=0; i<_NIcomp; i++ )
      {
        _Kappa[i] = KappaIn[i].GetDouble();
      }

    }
    else
    { 
      std::cout << "Error opening params_EnergyModel.in" << std::endl;
      exit(1);
    }
    f2.close();
  }
  
  else
  {//Legacy file format
    
    std::ifstream f2("params_EnergyModel.in");
    if ( f2.is_open() )
    {
      std::getline(f2, tmp_str, '#');
      _EnergyModelFlag = atoi(tmp_str.c_str());
      std::getline(f2, tmp_str, '\n');

      _Hess.resize(_NIcomp, Vec1dFieldType(_NIcomp, 0.));
      for( int i=0; i<_NIcomp; i++ )
      {
        for( int j=0; j<_NIcomp; j++ )
        {
          std::getline(f2, tmp_str, '#');
          _Hess[i][j] = FieldType(atof(tmp_str.c_str()), 0.);
          std::getline(f2, tmp_str, '\n');
        }
      }

      _Kappa.resize(_NIcomp, 0.);
      for( int i=0; i<_NIcomp; i++ )
      {
        std::getline(f2, tmp_str, '#');
        _Kappa[i] = atof(tmp_str.c_str());
        std::getline(f2, tmp_str, '\n');
      }

    }
    else
    { 
      std::cout << "Error opening params_EnergyModel.in" << std::endl;
      exit(1);
    }
    f2.close();
  
  }
 
  std::cout << "   - Energy Model Parameters:\n" << std::scientific;
  std::cout << "      EnergyModelFlag = " << _EnergyModelFlag << "\n";
  for( int i=0; i<_NIcomp; i++ )
  {
    for( int j=0; j<_NIcomp; j++ )
    {
      std::cout << "      H" << i << j << "             = ";
      std::cout << _Hess[i][j] << "\n";
    }
  }
  for( int i=0; i<_NIcomp; i++ )
  {
    std::cout << "      K" << i << "              = ";
    std::cout << _Kappa[i] << "\n";
  }

} // }}}

// Calculates the full chemical potential (in Fourier Space!)
void Model_Energy_Linear::calc_mu(const SmartFieldVec & phi, SmartFieldVec & mu)
{ // {{{ 

  // mu = sum_j ( df/dphi_j - K_ij del^2 phi_j )
  // mu = sum_j ( H_ij phi_j - K_ij del^2 phi_j )

  if (phi.getflag_inrealspace() != false)
  {
    std::cout << "*** Error in Model_Energy_Linear::calc_mu ***\n";
    std::cout << "*** phi must be in Fourier space ***\n";
  }

  mu.setflag_inrealspace(false);

  // get constant terms
  for (int i = 0; i < _NIcomp; i++)
  {
    mu[i] = phi[0];
    mu[i] *= _Hess[i][0];
    for (int j = 1; j < _NIcomp; j++)
    {
      mu[i].xpby_inplace(phi[j],_Hess[i][j]);
    }
  }

  //take Del2 of fields
  SmartFieldVec grad2_phi(2);
  grad2_phi.setflag_inrealspace(false);
  _OpObj->Del2_f_ex( phi, grad2_phi );

  // get graident terms
  for (int i = 0; i < _NIcomp; i++)
  {
    mu[i].xpby_inplace(grad2_phi[i], -_Kappa[i]);
  }

} // }}}

// Calculates the linear part of the chemical potential explicitly
void Model_Energy_Linear::calc_mu_lin_exp(const SmartFieldVec & phi, SmartFieldVec & mu_lin)
{ // {{{ 

  // mu = sum_j ( df/dphi_j - K_ij del^2 phi_j )
  // mu_lin = H_ij phi_j - K_ij del^2 phi_j
  // (in this case, both are the same)

  if (phi.getflag_inrealspace() != false)
  {
    std::cout << "*** Error in Model_Energy_Linear::calc_mu_lin_exp ***\n";
    std::cout << "*** phi must be in Fourier space ***\n";
  }

  calc_mu( phi, mu_lin ); // linear model, so same as above

} // }}}

// Calculates the linear part of the chemical potential implicitly
// Only returns the operator F, not mu, such that mu = F.phi
void Model_Energy_Linear::calc_mu_lin_imp(const SmartFieldVec & phi, SmartFieldOpMat & F)
{ // {{{ 

  // mu = sum_j ( df/dphi_j - K_ij del^2 phi_j )
  // mu_lin = H_ij phi_j - K_ij del^2 phi_j
  // (in this case, both are the same)

  F.zero(); // Must be initialized because Op functions ADD contributions
  F.setflag_inrealspace(false);

  // get gradient terms fist
  Vec2dFieldType K(_NIcomp, Vec1dFieldType(_NIcomp, FieldType(0.)));
  K = kroneckerdelta(_Kappa);
  _OpObj->A_Del2_im (K, F);
  F *= -1;

  // Add linear, diagonal terms
  F.AddBand(0, _Hess);

} // }}}

// Calculates the gradient of full chemical potential (in Fourier Space!)
void Model_Energy_Linear::calc_del_mu(const SmartFieldVec & phi, SmartFieldMat & del_mu)
{ // {{{ 

  // mu = sum_j ( df/dphi_j - K_ij del^2 phi_j )
  // mu = sum_j ( H_ij phi_j - K_ij del^2 phi_j )

  if (phi.getflag_inrealspace() != false)
  {
    std::cout << "*** Error in Model_Energy_Linear::calc_mu ***\n";
    std::cout << "*** phi must be in Fourier space ***\n";
  }

  del_mu.setflag_inrealspace(false);
  SmartFieldVec mu( _NIcomp);
  mu.setflag_inrealspace(false);

  // get constant terms
  for (int i = 0; i < _NIcomp; i++)
  {
    mu[i] = phi[0];
    mu[i] *= _Hess[i][0];
    for (int j = 1; j < _NIcomp; j++)
    {
      mu[i].xpby_inplace(phi[j],_Hess[i][j]);
    }
  }

  //take Del2 of fields
  SmartFieldVec grad2_phi(2);
  grad2_phi.setflag_inrealspace(false);
  _OpObj->Del2_f_ex( phi, grad2_phi );

  // get graident terms
  for (int i = 0; i < _NIcomp; i++)
  {
    mu[i].xpby_inplace(grad2_phi[i], -_Kappa[i]);
  }

  _OpObj->Del_f_ex( mu, del_mu );

} // }}} 

