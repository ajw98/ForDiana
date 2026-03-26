/***************************************************************
*
* Class defining a Ohta-Kawasaki-Doi-Uneyama free energy model
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

#include "Model_Energy_A_B_ABlinear.h"

// ----------------- Constructor/Destructor -----------------
Model_Energy_A_B_ABlinear :: Model_Energy_A_B_ABlinear( Operators_Base & OpObj, int ncomp, int nicomp ):
Model_Energy_Base( OpObj, ncomp, nicomp ),
_KappaGrad2Phi( nicomp ),
_PrecalcMuExp(false),
_PrecalcMuImp(false),
_RegFlag(false)
{ // {{{
} // }}}

Model_Energy_A_B_ABlinear :: ~Model_Energy_A_B_ABlinear()
{ // {{{
} // }}}

// ----------------- Input/Output -----------------
void Model_Energy_A_B_ABlinear :: set_params()
{ // {{{

  std::string tmp_str;

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
	
    if (_CReg != 0.)
      _RegFlag=true;      // determine if regularization is needed

    _N.resize(_Ncomp, 0.);
    _alpha.resize(_Ncomp, 0.);
    for( int i=0; i<_Ncomp; i++ )
    {
      std::getline(f2, tmp_str, '#');
      _N[i] = atof(tmp_str.c_str());
      std::getline(f2, tmp_str, '\n');

      _alpha[i] = _N[i]/_Nr;
    }

    _f = _N[1]/(_N[1]+_N[2]);
    _sf.resize(2, 0.);
    _sf = get_sf(_f);

    _alphaD = _alpha[1] + _alpha[2];

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

    std::getline(f2, tmp_str, '#');
    _KappaM = atof(tmp_str.c_str());
    std::getline(f2, tmp_str, '\n'); 
    
    _C.resize(_Ncomp, Vec1dReal(_Ncomp, 0.));
    for( int i=0; i<_Ncomp; i++ )
    {
      for( int j=i; j<_Ncomp; j++ )
      {
        std::getline(f2, tmp_str, '#');
        _C[i][j] = atof(tmp_str.c_str());
        std::getline(f2, tmp_str, '\n');

        _C[j][i] = _C[i][j];
      }
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

  std::cout << "      KM              = " << _KappaM  << "\n";

  for( int i=0; i<_Ncomp; i++ )
  {
    for( int j=i; j<_Ncomp; j++ )
    {
      std::cout << "      C" << i << j << "             = " << _C[i][j] << "\n";
    }
  }
 

} // }}}


// Calculates the full chemical potential
void Model_Energy_A_B_ABlinear::calc_mu(const SmartFieldVec & phi, SmartFieldVec & mu)
{ // {{{ 

  // (in) phi -- in fourier space
  // (out) mu -- in fourier space

  if (phi.getflag_inrealspace() != false)
  {
    std::cout << "*** Error in Model_Energy_A_B_ABlinear::calc_mu ***\n";
    std::cout << "*** phi must be in Fourier space ***\n";
  }

  // real space copy of phi
  SmartFieldVec phi_r(phi);
  phi_r.ifft();

  SmartField tmp1, tmp2;

  mu.setflag_inrealspace(true);
  tmp1.setflag_inrealspace(true);
  tmp2.setflag_inrealspace(true);

  // phi_M = 1 - phi_1 - phi_2 - ... phi_{M-1}
  SmartField logphi_M;
  logphi_M.setflag_inrealspace(true);
  SmartField phiM; // Fourier space copy of phi_M for nonlocal terms
  logphi_M = 1.;
  for (int i=0; i<_NIcomp; i++)
    logphi_M -= phi_r[i];
  phiM = logphi_M; // Extract phi_M to get Fourier space version
  logphi_M.logreal(logphi_M);
  logphi_M += FieldType(1.);
  logphi_M *= (1./_alpha[_Ncomp-1]);
  phiM.fft();
  phiM.setflag_inrealspace(false);

  // entropy terms
  // mu = 1/alpha_i (1 + log(phi_i)) - 1/alpha_M (1 + log(phi_M))
  for (int i = 0; i < _NIcomp; i++)
  {
    mu[i].logreal(phi_r[i]);
    mu[i] += FieldType(1.);

    if (i == 1)
      mu[i] *= (_sf[0]/_alpha[i]);
    else if (i == 2)
      mu[i] *= (_sf[1]/_alpha[i]);
    else
      mu[i] *= (1./_alpha[i]);

    mu[i] -= logphi_M;
  }
  
  // Diblock's surfactant terms
  // Set 1 = AD and 2 = BD, alpha[1] = f*Nd/Nr, alpha[2] = (1-f)*Nd/Nr
  // mu_1 -= sqrt(phi_2/phi_1) / (2 alpha_1)
  // mu_2 -= sqrt(phi_1/phi_2) / (2 alpha_2)
  // Use 0.4015 as multiplying factor for phi_1/2 log phi_1/2 with this
  
  tmp1.sqrt(phi_r[1]);
  tmp2.sqrt(phi_r[2]);
  tmp2 /= tmp1;
  tmp2 /= (2*_alphaD*std::sqrt(_f*(1.-_f)));
  mu[1] -= tmp2;
  
  tmp1.sqrt(phi_r[2]);
  tmp2.sqrt(phi_r[1]);
  tmp2 /= tmp1;
  tmp2 /= (2*_alphaD*std::sqrt(_f*(1.-_f)));
  mu[2] -= tmp2;  

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

  // regularization
  if (_RegFlag)
  {
    SmartField phi_M(phiM);
    phi_M.ifft();
    tmp1.exponentiate(phi_r[0], -1./_Delta);
    mu[0].xpby_inplace(tmp1, -1*(_CReg/_Delta));
    tmp1.exponentiate(phi_M, -1./_Delta);
    for (int i = 0; i < _NIcomp; i++)
    {
      mu[i].xpby_inplace(tmp1, (_CReg/_Delta));
    }
  }

  mu.fft();

  //take Del2 of fields
  _KappaGrad2Phi.setflag_inrealspace(false);
  _OpObj->Del2_f_ex( phi, _KappaGrad2Phi );
    
  SmartField KappaGrad2phiM;
  KappaGrad2phiM.setflag_inrealspace(false);
  _OpObj->Del2_f_ex( phiM, KappaGrad2phiM );
  FieldType phiM_mean = phiM.getaverage();
  KappaGrad2phiM *= (2./phiM_mean*_KappaM);

  // get gradient terms
  for (int i = 0; i < _NIcomp; i++)
  {
    FieldType phii_mean = phi[i].getaverage();
    _KappaGrad2Phi[i] *= (2./phii_mean*_Kappa[i]);
    _KappaGrad2Phi[i] -= KappaGrad2phiM;
    mu[i] -= _KappaGrad2Phi[i];
  }

  // Diblock's Coulomb terms
  // mu_i = integral_dr' ( G(r-r') * (1 - sqrt(phi_i (r')/phi_j(r')) ) ) where i != j
  
  FieldType phiD_mean = phi[1].getaverage()/_f;
  SmartField tau(phi[1]);
  tau.setflag_inrealspace( false );
  _Coul_tau.setflag_inrealspace( false );
  tau.axpby_inplace( phi[2], (1./_f), (-1./(1.-_f)) );
  _OpObj->Coulomb_f_ex( tau, _Coul_tau );
  _Coul_tau *= (_C[1][1]/(phiD_mean*_alphaD*_alphaD*2.));
  mu[1].xpby_inplace( _Coul_tau, (1./_f) );
  mu[2].xpby_inplace( _Coul_tau, (-1./(1.-_f)) ); 

  // get maximum of Hessian (for implicit/explicit parts)
  _Hmax = calc_Hmax( phi_r );
  _PrecalcMuExp=true;
  _PrecalcMuImp=true;

} // }}}


// Calculates the linear part of the chemical potential explicitly
void Model_Energy_A_B_ABlinear::calc_mu_lin_exp(const SmartFieldVec & phi, SmartFieldVec & mu_lin)
{ // {{{ 

  // mu = sum_j ( df/dphi_j - K_ij del^2 phi_j )
  // mu_lin = H_ij phi_j - K_ij del^2 phi_j

  if (phi.getflag_inrealspace() != false)
  {
    std::cout << "*** Error in Model_Energy_A_B_ABlinear::calc_mu_lin_exp ***\n";
    std::cout << "*** phi must be in Fourier space ***\n";
  }

  if(_PrecalcMuExp==false)
  {
    std::cout << "*** Error in Model_Energy_A_B_ABlinear::calc_mu_lin_exp ***\n";
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
  {
    mu_lin[i] -= _KappaGrad2Phi[i];
  }

  // get Coulomb terms
  mu_lin[1].xpby_inplace( _Coul_tau, (1./_f) );
  mu_lin[2].xpby_inplace( _Coul_tau, (-1./(1.-_f)) );

  _PrecalcMuExp = false;

} // }}}


// Calculates the linear part of the chemical potential implicitly
// Only returns the operator F, not mu, such that mu = F.phi
void Model_Energy_A_B_ABlinear::calc_mu_lin_imp(const SmartFieldVec & phi, SmartFieldOpMat & F)
{ // {{{ 

  if (_PrecalcMuImp==false)
  {
    std::cout << "*** Error in Model_Energy_A_B_ABlinear::calc_mu_lin_imp ***\n";
    std::cout << "*** Must call calc_mu before calc_mu_lin_imp to precalculate some quantitities. ***\n";
  }

  // mu = sum_j ( df/dphi_j - K_ij del^2 phi_j )
  // mu_lin = H_ij phi_j - K_ij del^2 phi_j
  // (in this case, both are the same)
  F.zero();
  F.setflag_inrealspace(false);

  // phi_M = 1 - phi_1 - phi_2 - ... phi_{M-1}
  FieldType phiM_mean;
  phiM_mean = 1.;
  for (int i=0; i<_NIcomp; i++)
    phiM_mean -= phi[i].getaverage();
  
  // get gradient terms first
  Vec2dFieldType K(_NIcomp, Vec1dFieldType(_NIcomp, FieldType(0.)));
  K = kroneckerdelta(_Kappa);
  for (int i = 0; i < _NIcomp; i++)
  {
    FieldType phii_mean = phi[i].getaverage();
    for (int j = 0; j < _NIcomp; j++)
    {
      K[i][j] *= (2./phii_mean);
      K[i][j] += (2.*_KappaM/phiM_mean);
    }
  }

  _OpObj->A_Del2_im (K, F);
  F *= -1;
  
  // get Coulomb terms next
  Vec2dFieldType CI(_NIcomp, Vec1dFieldType(_NIcomp, FieldType(0.)));
  FieldType phiD_mean = phi[1].getaverage()/_f;
  
  for (int i = 0; i < _NIcomp; i++)
  {
    for (int j = 0; j < _NIcomp; j++)
    {
      CI[i][j] = (_C[i][j] - _C[i][_NIcomp] - _C[_NIcomp][j])/(phiD_mean*_alpha[i]*_alpha[j]*2.);
    }
  }
  _OpObj->A_Coulomb_im (CI, F);   
  
  // Add linear, diagonal terms = maximum of Hessian
  F.AddBand(0, _Hmax);

  _PrecalcMuImp=false;

} // }}}

void Model_Energy_A_B_ABlinear::calc_Hessian( const SmartFieldVec & phi, SmartFieldMat & Hess )
{ // {{{

  if (phi.getflag_inrealspace() != true)
  {
    std::cout << "*** Error in Model_Energy_A_B_ABlinear::calc_Hmax ***\n";
    std::cout << "*** phi must be in Real space ***\n";
  }

  SmartField tmp1, tmp2, tmp3;
  tmp1.setflag_inrealspace(true);
  tmp2.setflag_inrealspace(true);
  tmp3.setflag_inrealspace(true);
  Hess.setflag_inrealspace(true);

  // Precalculate 1/(alpha_M phi_M) and 1/sqrt(phi_1 phi_2)/(4 sqrt(f*(1-f)) alpha_D)
  // since they appear repeatedly
  // phi_M = 1 - phi_1 - phi_2 - ... phi_{M-1}
  SmartField phi_M;
  phi_M.setflag_inrealspace(true);
  phi_M = 1.;
  for (int i=0; i<_NIcomp; i++)
    phi_M -= phi[i];
  tmp1 = 1./_alpha[_Ncomp-1];
  tmp1 /= phi_M;  

  tmp2.sqrt(phi[1]);
  tmp3.sqrt(phi[2]);
  tmp3 *= tmp2;
  tmp2 = FieldType( 1./(4.*std::sqrt(_alpha[1]*_alpha[2]) ) );
  tmp2 /= tmp3;

  for (int i = 0; i < _NIcomp; i++)
  {
    // Diagonal terms
    // H_ii = 1/(alpha_i phi_i) + 1/(alpha_M phi_M) - 2 ChiN_iM
    // H_11 += 1/sqrt(phi_1 phi_2)/(4 sqrt(f*(1-f)) alpha_D)  phi_2/phi_1
    // H_22 += 1/sqrt(phi_1 phi_2)/(4 sqrt(f*(1-f)) alpha_D)  phi_1/phi_2

    if (i == 1)
      Hess(i,i) = (_sf[0]/_alpha[i]);
    else if (i == 2)
      Hess(i,i) = (_sf[1]/_alpha[i]);
    else
      Hess(i,i) = (1./_alpha[i]);

    Hess(i,i) /= phi[i];
    Hess(i,i) += tmp1;
    Hess(i,i) -= (2.*_ChiN[i][_Ncomp-1]);

    if (i == 1)
    {
      tmp3 = tmp2;
      tmp3 *= phi[2];
      tmp3 /= phi[1];
      Hess(i,i) += tmp3;   
    }

    if (i == 2)
    {
      tmp3 = tmp2;
      tmp3 *= phi[1];
      tmp3 /= phi[2];
      Hess(i,i) += tmp3;   
    }

    for (int j=i+1; j < _NIcomp; j++)
    {
      // Off diagonal terms
      // H_ij = 1/(alpha_M phi_M) + ChiN_ijM
      // H_12, H_21 -= 1/sqrt(phi_1 phi_2)/(4 sqrt(f*(1-f)) alpha_D)
      Hess(i,j) = tmp1;
      Hess(i,j) += _ChiN_ijM[i][j];
      if (i == 1)
        if (j == 2)
          Hess(i,j) -= tmp2;

      Hess(j,i) = Hess(i,j); // Symmetric matrix
    }
  }

  // regularization
  if(_RegFlag)
  {
    tmp1.exponentiate(phi[0], -1./_Delta);
    Hess(0,0).xpby_inplace(tmp1, (_CReg/_Delta/_Delta));
    for (int i=0; i < _NIcomp; i++)
    {
      for (int j=0; j < _NIcomp; j++)
      { 
        tmp1.exponentiate(phi_M, -1./_Delta);
        Hess(i,j).xpby_inplace(tmp1, (_CReg/_Delta/_Delta));
      }
    }
  }

} // }}} 

Vec2dFieldType Model_Energy_A_B_ABlinear::calc_Hmax( const SmartFieldVec & phi )
{ // {{{

  if (phi.getflag_inrealspace() != true)
  {
    std::cout << "*** Error in Model_Energy_A_B_ABlinear::calc_Hmax ***\n";
    std::cout << "*** phi must be in Real space ***\n";
  }

  Vec2dFieldType Hmax( _NIcomp, Vec1dFieldType( _NIcomp, 0.) );

  // find max of Hessian
  SmartFieldMat Hess( _NIcomp, _NIcomp );
  Hess.setflag_inrealspace(true);
  calc_Hessian( phi, Hess );

  SmartField H_norm;
  UInt max_indx;

  Hess.norm2(H_norm, true); // Using norm^2, which is faster and does not change max c.f. norm()
  max_indx = H_norm.maxidx();
  Hmax = Hess.getelement(max_indx);

  // ensure that Hmax is positive definite
  make_posdef_ReSymm( Hmax );
   
  return Hmax;

} // }}}


// Calculates the full chemical potential gradient
void Model_Energy_A_B_ABlinear::calc_del_mu(const SmartFieldVec & phi, SmartFieldMat & del_mu)
{ // {{{ 

  // (in) phi -- in fourier space
  // (out) del_mu -- in real space

  if (phi.getflag_inrealspace() != false)
  {
    std::cout << "*** Error in Model_Energy_A_B_ABlinear::calc_del_mu ***\n";
    std::cout << "*** phi must be in Fourier space ***\n";
  }

  // real space copy of phi
  SmartFieldVec phi_r(phi);
  phi_r.ifft();

  SmartField phiM; // Fourier space copy of phi_M for nonlocal terms
  phiM.setflag_inrealspace(true);
  phiM = 1.;
  for (int i=0; i<_NIcomp; i++)
    phiM -= phi_r[i];
  phiM.fft();

  SmartFieldMat del_phi( _NIcomp, _Dim ); 
  del_phi.setflag_inrealspace(false);
  _OpObj->Del_f_ex( phi, del_phi ); 
  del_phi.ifft();

  del_mu.setflag_inrealspace(true);
  SmartFieldMat Hess( _NIcomp, _NIcomp );
  Hess.setflag_inrealspace(true);
  calc_Hessian( phi_r, Hess );

  del_mu.dot( Hess, del_phi );

  // take Del2 of fields
  // we could take Del3 of fields but having some nonlocal terms 
  // precomputed helps with explicit/implicit parts
  _KappaGrad2Phi.setflag_inrealspace(false);
  _OpObj->Del2_f_ex( phi, _KappaGrad2Phi );
    
  SmartField KappaGrad2phiM;
  KappaGrad2phiM.setflag_inrealspace(false);
  _OpObj->Del2_f_ex( phiM, KappaGrad2phiM );
  FieldType phiM_mean = phiM.getaverage();
  KappaGrad2phiM *= (-2./phiM_mean*_KappaM);

  SmartFieldVec tmpvec(_Dim);
  // get gradient terms
  for (int i = 0; i < _NIcomp; i++)
  {
    tmpvec.setflag_inrealspace(false);
    FieldType phii_mean = phi[i].getaverage();
    _KappaGrad2Phi[i] *= (2./phii_mean*_Kappa[i]);
    _KappaGrad2Phi[i] += KappaGrad2phiM;
    _OpObj->Del_f_ex( _KappaGrad2Phi[i], tmpvec );
    tmpvec.ifft();
    for (int j = 0; j < _Dim; j++)
      del_mu(i,j) -= tmpvec[j];
  }

  // Diblock's Coulomb terms
  // mu_i = integral_dr' ( G(r-r') * (1 - sqrt(phi_i (r')/phi_j(r')) ) ) where i != j
  tmpvec.setflag_inrealspace(false);
  FieldType phiD_mean = phi[1].getaverage()/_f;
  _Coul_tau.setflag_inrealspace( false );
  SmartField tau(phi[1]);
  tau.setflag_inrealspace(false);
  tau.axpby_inplace( phi[2], (1./_f), (-1./(1.-_f)) );
  _OpObj->Coulomb_f_ex( tau, _Coul_tau );
  _Coul_tau *= (_C[1][1]/(phiD_mean*_alphaD*_alphaD*2.));
  _OpObj->Del_f_ex( _Coul_tau, tmpvec);
  tmpvec.ifft();

  for (int j=0; j < _Dim; j++)
  {
    del_mu(1,j).xpby_inplace( tmpvec[j], (1./_f) );
    del_mu(2,j).xpby_inplace( tmpvec[j], (-1./(1.-_f)) ); 
  } 

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


// Get the coefficients of the entropy terms for diblock
// This correlation is obtained using Python for C(f) = 2/f*(s(f)+0.25)
// C*f^2 = 0.36321*f^4 - 0.7217*f^3 + 0.301*f^2 - 1.81178*f + 1.02972 -1.5*f*log(f))
// mean error 0.1%, max error 2.23% at f = 0.74

Vec1dReal Model_Energy_A_B_ABlinear::get_sf( const RealType &f )
{ // {{{

  Vec1dReal sf(2) ;

  if (f == 0.5)
  {
    sf[1] = 0.4015;
    sf[0] = 0.4015;
  }
  else
  {
    sf[0] = (0.36321*f*f*f - 0.7217*f*f + 0.301*f - 1.81178 + 1.02972/f -1.5*std::log(f))/2 - 0.25;
    sf[1] = (0.36321*(1.-f)*(1.-f)*(1.-f)-0.7217*(1.-f)*(1.-f) + 0.301*(1.-f)-1.81178+1.02972/(1.-f) -1.5*std::log(1.-f))/2-0.25;
  }

  return sf;

} // }}}

