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

#include "Model_Reactions_None.h"

// ----------------- Constructor/Destructor -----------------
Model_Reactions_None :: Model_Reactions_None( int ncomp, int nicomp ):
Model_Reactions_Base( ncomp, nicomp )
{ // {{{
} // }}}

Model_Reactions_None :: ~Model_Reactions_None()
{ // {{{
} // }}}

// ----------------- Input/Output -----------------
void Model_Reactions_None :: set_params()
{ // {{{
  
  if (jsonFile("params_Reactions.in"))
  {
    std::ifstream f2("params_Reactions.in");
    if ( f2.is_open() )
    {
      rapidjson::IStreamWrapper isw(f2);
      rapidjson::Document d;
      d.ParseStream(isw);
      _ReactionsModelFlag = d["ReactionsModelFlag"].GetInt();
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
    std::string tmp_str;
    std::ifstream f2("params_Reactions.in");
    if ( f2.is_open() )
    {
      std::getline(f2, tmp_str, '#');
      _ReactionsModelFlag = atoi(tmp_str.c_str());
      std::getline(f2, tmp_str, '\n');
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
   
} // }}}

// Calculates the full reactions
void Model_Reactions_None::calc_reactions(const SmartFieldVec & phi, SmartFieldVec & rxn)
{ // {{{ 

  // should never get called, but it case it is, zero it.
  rxn.setflag_inrealspace(false);
  rxn.zero();

} // }}}

// Calculates the linear part of the reaction explicitly
void Model_Reactions_None::calc_reactions_lin_exp(const SmartFieldVec & phi, SmartFieldVec & rxn_lin)
{ // {{{ 

  // should never get called, but it case it is, zero it.
  rxn_lin.setflag_inrealspace(false);
  rxn_lin.zero();

} // }}}

// Calculates the linear part of the reaction rate implicitly
// Only returns the operator Rxn, not d/dt (phi), such that d/dt (phi) = Rxn.phi
void Model_Reactions_None::calc_reactions_lin_imp(const SmartFieldVec & phi, SmartFieldOpMat & Rxn)
{ // {{{ 

  // should never get called, but it case it is, zero it.
  Rxn.setflag_inrealspace(false);
  Rxn.zero();

} // }}}

