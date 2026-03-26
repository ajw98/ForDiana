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

#include "Model_Reactions_Binary.h"

// ----------------- Constructor/Destructor -----------------
Model_Reactions_Binary :: Model_Reactions_Binary( int ncomp, int nicomp ):
Model_Reactions_Base( ncomp, nicomp )
{ // {{{
} // }}}

Model_Reactions_Binary :: ~Model_Reactions_Binary()
{ // {{{
} // }}}

// ----------------- Input/Output -----------------
void Model_Reactions_Binary :: set_params()
{ // {{{

  std::string tmp_str;

  if (jsonFile("params_Reactions.in"))
  {
    std::ifstream f2("params_Reactions.in");
    if ( f2.is_open() )
    {
      rapidjson::IStreamWrapper isw(f2);
      rapidjson::Document d;
      d.ParseStream(isw);

      _ReactionsModelFlag = d["ReactionsModelFlag"].GetInt();

      _kf = d["kf"].GetDouble();

      _kb = d["kb"].GetDouble();

      _Nr = d["Nr"].GetDouble();

      _N.resize(_Ncomp, 0.);
      rapidjson::Value& N_in = d["N"];
      for( int i=0; i<_Ncomp; i++ )
      {
        _N[i] = N_in[i].GetDouble();
      }

      _n = d["n"].GetInt();

      _Linratemax.resize( _NIcomp, Vec1dFieldType(_NIcomp, 0.) ); 

    }
    else
    { 
      std::cout << "Error opening params_Reactions.in" << std::endl;
      exit(1);
    }
    f2.close();
  }
  
  else
  {//Legacy file format
    std::ifstream f2("params_Reactions.in");
    if ( f2.is_open() )
    {
      std::getline(f2, tmp_str, '#');
      _ReactionsModelFlag = atoi(tmp_str.c_str());
      std::getline(f2, tmp_str, '\n');

      std::getline(f2, tmp_str, '#');
      _kf = atof(tmp_str.c_str());
      std::getline(f2, tmp_str, '\n');

      std::getline(f2, tmp_str, '#');
      _kb = atof(tmp_str.c_str());
      std::getline(f2, tmp_str, '\n');

      std::getline(f2, tmp_str, '#');
      _Nr = atof(tmp_str.c_str());
      std::getline(f2, tmp_str, '\n');

      _N.resize(_Ncomp, 0.);
      for( int i=0; i<_Ncomp; i++ )
      {
        std::getline(f2, tmp_str, '#');
        _N[i] = atof(tmp_str.c_str());
        std::getline(f2, tmp_str, '\n');
      }

      std::getline(f2, tmp_str, '#');
      _n = atoi(tmp_str.c_str());
      std::getline(f2, tmp_str, '\n');

      _Linratemax.resize( _NIcomp, Vec1dFieldType(_NIcomp, 0.) ); 
    }
    else
    { 
      std::cout << "Error opening params_Reactions.in" << std::endl;
      exit(1);
    }
    f2.close();
  }

  std::cout << "   - Reaction Parameters:\n" << std::scientific;
  std::cout << "      ReactionsModelFlag        = " << _ReactionsModelFlag << "\n";
  std::cout << "      Forward rate kf           = " << _kf << "\n";
  std::cout << "      Backward rate kb          = " << _kb << "\n";
  std::cout << "      Homopolymer/diblock ratio = " << _n  << "\n";
   
} // }}}

// Calculates the full reactions
void Model_Reactions_Binary::calc_reactions(const SmartFieldVec & phi, SmartFieldVec & rxn)
{ // {{{ 

  // (in) phi -- in fourier space
  // (out) d/dt( phi ) -- in fourier space

  if (phi.getflag_inrealspace() != false)
  {
    std::cout << "*** Error in Model_Reactions::calc_reaction_lin_exp ***\n";
    std::cout << "*** phi must be in Fourier space ***\n";
  }

  // real space copy of phi
  SmartFieldVec phi_r(phi);
  phi_r.ifft();

  rxn.zero();
  rxn.setflag_inrealspace(true);
  
  // phi_M = 1 - phi_1 - phi_2 - ... phi_{M-1}
  SmartField phi_M;
  phi_M.setflag_inrealspace(true);
  phi_M = 1.;
  for (int i=0; i<_NIcomp; i++)
    phi_M -= phi_r[i];
    
  // Stoichiometric reaction: A_h + B_h <---> n * A_d-B_d
  //
  // n = 1: rxn[1] = - rxn[0] = ( (kf/N_3) phi_0 phi_M - kb phi_1 ) Nr
  //        rxn[2] = - rxn[3] = ( (kf/N_0) phi_0 phi_M - kb phi_2 ) Nr
  //
  // n = 2: rxn[1] = - rxn[0] = ( (kf/N_3) phi_0 phi_M - kb*n/(N_1+N_2) (phi_1+phi_2) phi_1 ) Nr
  //        rxn[2] = - rxn[3] = ( (kf/N_0) phi_0 phi_M - kb*n/(N_1+N_2) (phi_1+phi_2) phi_2 ) Nr
  //
  // rxn[3] is not returned since component 3 is implicit
  
  // First get the forward rates if nonzero
  if ( _kf != 0. )
  {
    rxn[1].accumulateproduct_inplace(phi_r[0], phi_M, (_kf/_N[3]*_Nr));
    rxn[2].xpby_inplace(rxn[1], (_N[3]/_N[0]));
  }

  // Next get the backward rates if nonzero
  if (_kb != 0.)
  {
    if ( _n == 1 )
    {
      rxn[1].xpby_inplace(phi_r[1], (-1.*_kb*_Nr));
      rxn[2].xpby_inplace(phi_r[2], (-1.*_kb*_Nr));
    }
    else if ( _n == 2 )
    {
      SmartField tmp;
      tmp.setflag_inrealspace(true);
      tmp = phi_r[1];
      tmp += phi_r[2];
      RealType Kb;
      Kb = 2.*_kb*_Nr/(_N[1]+_N[2]);

      rxn[1].accumulateproduct_inplace(phi_r[1], tmp, (-1.*Kb));
      rxn[2].accumulateproduct_inplace(phi_r[2], tmp, (-1.*Kb));
    }
    else
    {
      std::cout << "Error. Check reaction parameters again" << std::endl;
    }
  }

  // Ensure conservation of mass between Homopolymer A and Diblock's A block
  rxn[0] -= rxn[1];
  rxn.fft();

  // get maximum of the linear rate for implicit/explicit parts
  _Linratemax = calc_linearratemax(phi_r);
  _PrecalcLinExp = true;
  _PrecalcLinImp = true;

} // }}}

// Calculates the linear part of the reaction explicitly
void Model_Reactions_Binary::calc_reactions_lin_exp(const SmartFieldVec & phi, SmartFieldVec & rxn_lin)
{ // {{{ 

  // rxn_lin = linratemax_ij phi_j

  if (phi.getflag_inrealspace() != false)
  {
    std::cout << "*** Error in Model_Reactions::calc_reaction_lin_exp ***\n";
    std::cout << "*** phi must be in Fourier space ***\n";
  }

  if (_PrecalcLinExp == false)
  {
    std::cout << "*** Error in Model_Reactions::calc_reactions_lin_exp ***\n";
    std::cout << "*** Must call calc_reactions before calc_reactions_lin_exp to precalculate some quantities. ***\n";
  }

  // Find mu_lin using Hmax
  rxn_lin.setflag_inrealspace(false);

  for (int i = 0; i < _NIcomp; i++)
  {
    rxn_lin[i] = phi[0];
    rxn_lin[i] *= _Linratemax[i][0];
    for (int j = 1; j < _NIcomp; j++)
    {
      rxn_lin[i].xpby_inplace(phi[j], _Linratemax[i][j]);
    }
  }
  
  _PrecalcLinExp = false;
     
} // }}}

// Calculates the linear part of the reaction rate implicitly
// Only returns the operator Rxn, not d/dt (phi), such that d/dt (phi) = Rxn.phi
void Model_Reactions_Binary::calc_reactions_lin_imp(const SmartFieldVec & phi, SmartFieldOpMat & Rxn)
{ // {{{ 

  // Rxn_ij = linratemax_ij

  if (phi.getflag_inrealspace() != false)
  {
    std::cout << "*** Error in Model_Reactions::calc_reaction_lin_imp ***\n";
    std::cout << "*** phi must be in Fourier space ***\n";
  }

  if (_PrecalcLinImp == false)
  {
    std::cout << "*** Error in Model_Reactions::calc_reactions_lin_imp ***\n";
    std::cout << "*** Must call calc_reactions before calc_reactions_lin_imp to precalculate some quantities. ***\n";
  }

  Rxn.zero();
  Rxn.setflag_inrealspace(false);
  Rxn.AddBand(0, _Linratemax);

  _PrecalcLinImp = false;

} // }}}

void Model_Reactions_Binary::calc_linearrate( const SmartFieldVec & phi, SmartFieldMat & linearrate )
{ // {{{

  if (phi.getflag_inrealspace() != true)
  {
    std::cout << "*** Error in Model_Reactions::calc_linearrate ***\n";
    std::cout << "*** phi must be in Real space ***\n";
  }

  linearrate.zero();
  linearrate.setflag_inrealspace(true);
 
  // phi_M = 1 - phi_1 - phi_2 - ... phi_{M-1}
  SmartField phi_M;
  phi_M.setflag_inrealspace(true);
  phi_M = 1.;
  for (int i=0; i<_NIcomp; i++)
    phi_M -= phi[i];
  
  // Forward reaction:
  // linearrate(1,0) = kf/N_3*Nr*(phi_M-phi_0)
  // linearrate(1,1) = linearrate(1,2) = -kf/N_3*Nr*phi_0
  // linearrate(0,i) = -linearrate(1,i)
  // linearrate(2,i) = N_3/N_0 linearrate(1,i)

  if (_kf != 0.)  
  {
    linearrate(1,0) = phi_M;
    linearrate(1,0).axpby_inplace( phi[0], (_kf*_Nr/_N[3]) , (-1.*_kf*_Nr/_N[3]) );
    linearrate(1,1).xpby_inplace( phi[0], (-1.*_kf*_Nr/_N[3]) );
    linearrate(1,2) = linearrate(1,1);
    for (int i=0; i<_NIcomp; i++)
    {
      linearrate(2,i).xpby_inplace(linearrate(1,i), (_N[3]/_N[0]) );
    }  
  }

  // Backward reaction:
  // n = 1: linearrate(1,1) = -kb*Nr  
  //        linearrate(2,2) = -kb*Nr  
  //        linearrate(0,1) = -linearrate(1,1)
  //
  // n = 2: linearrate(1,1) = -2*kb*Nr/(N_1+N_2) (2 phi_1 + phi_2)  
  //        linearrate(1,2) = -2*kb*Nr/(N_1+N_2) phi_1
  //        linearrate(2,2) = -2*kb*Nr/(N_1+N_2) (phi_1 + 2 phi_2)
  //        linearrate(2,1) = -2*kb*Nr/(N_1+N_2) phi_2
  //        linearrate(0,1) = 2*kb*Nr/(N_1+N_2) (2 phi_1 + phi_2)
  //        linearrate(0,2) = 2*kb*Nr/(N_1+N_2) phi_1
 
  if (_kb != 0.)
  {
    if ( _n == 1 )
    {
      linearrate(1,1) -= FieldType(_kb*_Nr);
      linearrate(2,2) -= FieldType(_kb*_Nr);
    }
    else if ( _n == 2 )
    {
      RealType Kb;
      Kb = 2.*_kb*_Nr/(_N[1]+_N[2]);
      linearrate(1,1).xpby_inplace(phi[1], (-2.*Kb));
      linearrate(1,1).xpby_inplace(phi[2], (-1.*Kb));
      linearrate(1,2).xpby_inplace(phi[1], (-1.*Kb));
      linearrate(2,1).xpby_inplace(phi[2], (-1.*Kb));
      linearrate(2,2).xpby_inplace(phi[1], (-1.*Kb));
      linearrate(2,2).xpby_inplace(phi[2], (-2.*Kb));
    } 
    else
    {
      std::cout << "Error. Check reaction parameters again" << std::endl;
    }
  } 

  // Enforce mass conservation: linearrate(0,i) = -linearrate(1,i)
  for (int i=0; i<_NIcomp; i++)
  {
    linearrate(0,i) -= linearrate(1,i);
  } 
   
    
} // }}} 

Vec2dFieldType Model_Reactions_Binary::calc_linearratemax( const SmartFieldVec & phi )
{ // {{{

  if (phi.getflag_inrealspace() != true)
  {
    std::cout << "*** Error in Model_Reactions::calc_linratemax ***\n";
    std::cout << "*** phi must be in Real space ***\n";
  }

  Vec2dFieldType linratemax( _NIcomp, Vec1dFieldType( _NIcomp, 0.) );

  // find max of Hessian
  SmartFieldMat linrate( _NIcomp, _NIcomp );
  linrate.setflag_inrealspace(true);
  calc_linearrate( phi, linrate );

  SmartField linrate_norm;
  UInt max_indx;

  linrate.norm(linrate_norm);
  max_indx = linrate_norm.maxidx();
  linratemax = linrate.getelement(max_indx);
  
  return linratemax;

} // }}}

