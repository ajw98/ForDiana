/***************************************************************
*
* Class defining a Flory Huggins free energy model
* 
* DRT -- Wed, 18 Aug 2015
* DRT -- Wed, 05 Apr 2017
*
****************************************************************/
//
// mu_i = 1/alpha_i ( 1 + ln( phi_i ) ) ) - 1/alpha_M ( 1 + ln( phi_M ) )
//        + (1-2 phi_i) ChiN_iM + sum_{j \ne i}^{1, M-1} ( ChiN_ijM  phi_j )
//        - K_i del^2 phi_i
//
// where M = number of components, Ncomp
//

#include "Model_Energy_FHG.h"

// ----------------- Constructor/Destructor -----------------
Model_Energy_FHG :: Model_Energy_FHG( Operators_Base & OpObj, int ncomp, int nicomp ):
Model_Energy_Base( OpObj, ncomp, nicomp ),
_PrecalcMuExp(false),
_PrecalcMuImp(false),
_KappaGrad2Phi( nicomp )
{ // {{{
} // }}}

Model_Energy_FHG :: ~Model_Energy_FHG()
{ // {{{
} // }}}

// ----------------- Input/Output -----------------
void Model_Energy_FHG :: set_params()
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
      
      _Nr = d["Nr"].GetDouble();
      _CReg = d["C_reg"].GetDouble();
      _Delta = d["delta"].GetDouble();
      
      _N.resize(_Ncomp, 0.);
      _alpha.resize(_Ncomp, 0.);
      const rapidjson::Value& N_in = d["N"];
      for( int i=0; i<_Ncomp; i++ )
      {
        
        _N[i] = N_in[i].GetDouble();
        

        _alpha[i] = _N[i]/_Nr;
      }

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
      
      _Kappa.resize(_NIcomp, 0.);
      const rapidjson::Value& Kappa_in = d["Kappa"];
      for( int i=0; i<_NIcomp; i++ )
      {
        _Kappa[i] = Kappa_in[i].GetDouble();
      }


      _Hmax.resize(_NIcomp, Vec1dFieldType(_NIcomp, FieldType(0.)));

      _Dim = _OpObj->GetGrid()->Dim();

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

      std::getline(f2, tmp_str, '#');
      _Nr = atof(tmp_str.c_str());
      std::getline(f2, tmp_str, '\n');

      std::getline(f2, tmp_str, '#');
      _CReg = atof(tmp_str.c_str());
      std::getline(f2, tmp_str, '\n');

      std::getline(f2, tmp_str, '#');
      _Delta = atof(tmp_str.c_str());
      std::getline(f2, tmp_str, '\n');

      _N.resize(_Ncomp, 0.);
      _alpha.resize(_Ncomp, 0.);
      for( int i=0; i<_Ncomp; i++ )
      {
        std::getline(f2, tmp_str, '#');
        _N[i] = atof(tmp_str.c_str());
        std::getline(f2, tmp_str, '\n');

        _alpha[i] = _N[i]/_Nr;
      }
    
      _Chi.resize(_Ncomp, Vec1dReal(_Ncomp, 0.));
      _ChiN.resize(_Ncomp, Vec1dReal(_Ncomp, 0.));
      for( int i=0; i<_Ncomp; i++ )
      {
        for( int j=i+1; j<_Ncomp; j++ )
        {
          std::getline(f2, tmp_str, '#');
          _Chi[i][j] = atof(tmp_str.c_str());
          std::getline(f2, tmp_str, '\n');
        
          _Chi[j][i] = _Chi[i][j];
          _ChiN[i][j] = _Chi[i][j]*_Nr;
          _ChiN[j][i] = _ChiN[i][j];
        }
      }

      _Kappa.resize(_NIcomp, 0.);
      for( int i=0; i < _NIcomp; i++ )
      {
        std::getline(f2, tmp_str, '#');
        _Kappa[i] = atof(tmp_str.c_str());
        std::getline(f2, tmp_str, '\n');
      }

      _Hmax.resize(_NIcomp, Vec1dFieldType(_NIcomp, FieldType(0.)));

      _Dim = _OpObj->GetGrid()->Dim();

    }
    else
    { 
      std::cout << "Error opening params_EnergyModel.in" << std::endl;
      exit(1);
    }
    f2.close();
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

  std::cout << "   - Energy Model Parameters:\n" << std::scientific;
  std::cout << "      EnergyModelFlag = " << _EnergyModelFlag << "\n";
  std::cout << "      Nr              = " << _Nr << "\n";
  std::cout << "      C_reg           = " << _CReg  << "\n";
  std::cout << "      delta           = " << _Delta  << "\n";

  for( int i=0; i<_Ncomp; i++ )
  {
    std::cout << "      N" << i << "              = " << _N[i] << "\n";
  }

  for( int i=0; i<_Ncomp; i++ )
  {
    for( int j=i+1; j<_Ncomp; j++ )
    {
      std::cout << "      Chi" << i << j << "           = " << _Chi[i][j] << "\n";
    }
  }

  for( int i=0; i<_NIcomp; i++ )
  {
    std::cout << "      K" << i << "              = " << _Kappa[i] << "\n";
  }

} // }}}

// Calculates the full chemical potential
void Model_Energy_FHG::calc_mu(const SmartFieldVec & phi, SmartFieldVec & mu)
{ // {{{ 

  // (in) phi -- in fourier space
  // (out) mu -- in fourier space

  if (phi.getflag_inrealspace() != false)
  {
    std::cout << "*** Error in Model_Energy_FHG::calc_mu ***\n";
    std::cout << "*** phi must be in Fourier space ***\n";
  }

  // real space copy of phi
  SmartFieldVec phi_r(phi);
  phi_r.ifft();

  mu.setflag_inrealspace(true);

  // phi_M = 1 - phi_1 - phi_2 - ... phi_{M-1}
  SmartField logphi_M;
  logphi_M.setflag_inrealspace(true);
  logphi_M = 1.;
  for (int i=0; i<_NIcomp; i++)
    logphi_M -= phi_r[i];
  logphi_M.logreal(logphi_M);
  logphi_M += FieldType(1.);
  logphi_M *= (1./_alpha[_Ncomp-1]);

  // entropy terms
  // mu = 1/alpha_i (1 + log(phi_i)) - 1/alpha_M (1 + log(phi_M))
  for (int i = 0; i < _NIcomp; i++)
  {
    mu[i].logreal(phi_r[i]);
    mu[i] += FieldType(1.);
    mu[i] *= (1./_alpha[i]);

    mu[i] -= logphi_M;
  }

  // enthalpy terms
  // mu += (1-2 phi_i) ChiN_iM + sum_{j \ne i}^{1, M-1} ( ChiN_ijM  phi_j )
  for (int i = 0; i < _NIcomp; i++)
  {
    mu[i] += _ChiN[i][_Ncomp-1];
    mu[i].xpby_inplace(phi_r[i], -2.*_ChiN[i][_Ncomp-1]);

    for (int j = 0; j < _NIcomp; j++)
    {
      if (j != i)
        mu[i].xpby_inplace(phi_r[j], _ChiN_ijM[i][j]);
    }
  }

  mu.fft();

  // take Del2 of fields
  _KappaGrad2Phi.setflag_inrealspace(false);
  _OpObj->Del2_f_ex( phi, _KappaGrad2Phi ); // Up to now, KappaGrad2Phi contains only Grad2Phi. Multiply Kappa below.

  // get gradient terms
  for (int i = 0; i < _NIcomp; i++)
  {
    _KappaGrad2Phi[i] *= _Kappa[i];
    mu[i] -= _KappaGrad2Phi[i];
    //mu[i].xpby_inplace(_KappaGrad2Phi[i], -Kappa[i]); // Could use this, but _KappaGrad2Phi needs to be updated for linearized mu
  }

  // get maximum of Hessian (for implicit/explicit parts)
  _Hmax = calc_Hmax( phi_r );

  _PrecalcMuExp=true;
  _PrecalcMuImp=true;

} // }}}

// Calculates the linear part of the chemical potential explicitly
void Model_Energy_FHG::calc_mu_lin_exp(const SmartFieldVec & phi, SmartFieldVec & mu_lin)
{ // {{{ 

  // mu = sum_j ( df/dphi_j - K_ij del^2 phi_j )
  // mu_lin = H_ij phi_j - K_ij del^2 phi_j

  if (phi.getflag_inrealspace() != false)
  {
    std::cout << "*** Error in Model_Energy_FHG::calc_mu_lin_exp ***\n";
    std::cout << "*** phi must be in Fourier space ***\n";
  }

  if (_PrecalcMuExp==false)
  {
    std::cout << "*** Error in Model_Energy_FHG::calc_mu_lin_exp ***\n";
    std::cout << "*** Must call calc_mu before calc_mu_lin_exp to precalculate some quantities. ***\n";
  }

  // Find mu_lin using Hmax
  mu_lin.setflag_inrealspace(false);

  for (int i = 0; i < _NIcomp; i++)
  {
    mu_lin[i] = phi[0];
    mu_lin[i] *= _Hmax[i][0];
    for (int j = 1; j < _NIcomp; j++)
    {
      mu_lin[i].xpby_inplace(phi[j],_Hmax[i][j]);
    }
  }

  // subtract gradient terms
  for (int i = 0; i < _NIcomp; i++)
    mu_lin[i] -= _KappaGrad2Phi[i];

  _PrecalcMuExp=false;

} // }}}

// Calculates the linear part of the chemical potential implicitly
// Only returns the operator F, not mu, such that mu = F.phi
void Model_Energy_FHG::calc_mu_lin_imp(const SmartFieldVec & phi, SmartFieldOpMat & F)
{ // {{{ 

  if (_PrecalcMuImp==false)
  {
    std::cout << "*** Error in Model_Energy_FHG::calc_mu_lin_imp ***\n";
    std::cout << "*** Must call calc_mu before calc_mu_lin_imp to precalculate some quantities. ***\n";
  }

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

  //SmartFieldMat tmpmat(_NIcomp, _NIcomp);
  //tmpmat.setflag_inrealspace(false);
  //tmpmat = _Hmax;
  //F.AddBand(0, tmpmat);
  F.AddBand(0, _Hmax);

  _PrecalcMuImp=false;

} // }}}

void Model_Energy_FHG::calc_Hessian( const SmartFieldVec & phi, SmartFieldMat & Hess )
{ // {{{

  if (phi.getflag_inrealspace() != true)
  {
    std::cout << "*** Error in Model_Energy_FHG::calc_Hessian ***\n";
    std::cout << "*** phi must be in Real space ***\n";
  }

  SmartField tmp1;
  tmp1.setflag_inrealspace(true);
  Hess.setflag_inrealspace(true);

  // phi_M = 1 - phi_1 - phi_2 - ... phi_{M-1}
  SmartField phi_M;
  phi_M.setflag_inrealspace(true);
  phi_M = 1.;
  for (int i=0; i<_NIcomp; i++)
    phi_M -= phi[i];
  tmp1 = 1./_alpha[_Ncomp-1];
  tmp1 /= phi_M;

  for (int i = 0; i < _NIcomp; i++)
  {
    // Diagonal elements
    // H_ii = 1/(alpha_i phi_i) + 1/(alpha_M phi_M) - 2 ChiN_iM
    Hess(i,i) = 1./_alpha[i];
    Hess(i,i) /= phi[i];

    Hess(i,i) += tmp1;

    Hess(i,i) += -2*_ChiN[i][_Ncomp-1];

    // Off diagonals
    // H_ij = 1/(alpha_M phi_M) + 2 ChiN_ijM
    for (int j = i+1; j < _NIcomp; j++)
    {
      Hess(i,j) = tmp1;
      Hess(i,j) += _ChiN_ijM[i][j];
      Hess(j,i) = Hess(i,j); // Symmetric matrix
    }
  }

} // }}}

Vec2dFieldType Model_Energy_FHG::calc_Hmax( const SmartFieldVec & phi )
{ // {{{

  if (phi.getflag_inrealspace() != true)
  {
    std::cout << "*** Error in Model_Energy_FHG::calc_Hmax ***\n";
    std::cout << "*** phi must be in Real space ***\n";
  }

  Vec2dFieldType Hmax( _NIcomp, Vec1dFieldType( _NIcomp, 0.) );

  // find max of Hessian
  SmartFieldMat Hess( _NIcomp, _NIcomp );
  Hess.setflag_inrealspace(true);
  calc_Hessian( phi, Hess );

  SmartField H_norm;
  UInt max_indx;

//  Hess.norm(H_norm, true);
  Hess.norm2(H_norm, true); // Using norm^2, which is faster and does not change max c.f. norm()
  max_indx = H_norm.maxidx();
  Hmax = Hess.getelement(max_indx);

  // ensure that Hmax is positive definite
  make_posdef_ReSymm( Hmax );
  return Hmax;

} // }}}

// Calculates the gradient of full chemical potential
void Model_Energy_FHG::calc_del_mu(const SmartFieldVec & phi, SmartFieldMat & del_mu)
{ // {{{ 

  // (in) phi -- in fourier space
  // (out) del_mu -- in fourier space

  if (phi.getflag_inrealspace() != false)
  {
    std::cout << "*** Error in Model_Energy_FHG::calc_mu ***\n";
    std::cout << "*** phi must be in Fourier space ***\n";
  }

  // real space copy of phi
  SmartFieldVec phi_r(phi);
  phi_r.ifft();

  SmartFieldMat del_phi( _NIcomp, _Dim );
  del_phi.setflag_inrealspace(false);
  _OpObj->Del_f_ex( phi, del_phi );
  del_phi.ifft();

  del_mu.setflag_inrealspace(true);
  SmartFieldMat Hess( _NIcomp, _NIcomp );
  Hess.setflag_inrealspace(true);
  calc_Hessian( phi_r, Hess );

  del_mu.dot( Hess, del_phi );
  del_mu.fft();

  // take Del2 of fields
  _KappaGrad2Phi.setflag_inrealspace(false);
  _OpObj->Del2_f_ex( phi, _KappaGrad2Phi ); // Up to now, KappaGrad2Phi contains only Grad2Phi. Multiply Kappa below.

  SmartFieldVec tmpvec(_Dim);
  tmpvec.setflag_inrealspace(false); 

  // get gradient terms
  for (int i = 0; i < _NIcomp; i++)
  {
    _KappaGrad2Phi[i] *= _Kappa[i];
    _OpObj->Del_f_ex( _KappaGrad2Phi[i], tmpvec );
    for (int j = 0; j < _Dim; j++)
      del_mu(i,j) -= tmpvec[j];
  }

  del_mu.ifft();

  // get maximum of Hessian (for implicit/explicit parts)
  SmartField H_norm;
  UInt max_indx;
  Hess.norm2(H_norm, true); // Using norm^2, which is faster and does not change max c.f. norm()
  max_indx = H_norm.maxidx();
  _Hmax = Hess.getelement(max_indx);

  // ensure that Hmax is positive definite
  make_posdef_ReSymm( _Hmax );

  _PrecalcMuExp=true;
  _PrecalcMuImp=true;

} // }}} 

