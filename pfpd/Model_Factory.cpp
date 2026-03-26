/***************************************************************
*
* Base class and factory for the energy model
* 
* DRT -- Wed, 18 Aug 2015
* RA  -- Fri, 20 Dec 2019
*
****************************************************************/

#include "Model_Factory.h"

// ----------------- Constructor/Destructor -----------------

Model_Factory :: Model_Factory()
{ // {{{
} // }}}

Model_Factory :: ~Model_Factory()
{ // {{{
} // }}}

// ----------------- Factory -----------------

Model_Energy_Base* Model_Factory :: make_energy_model( int key, Operators_Base & OpObj, int Ncomp, int NIcomp )
{ // {{{

  switch(key) 
  {

    case 0: // Constant free energy (diffusion only)
      return new Model_Energy_Linear( OpObj, Ncomp, NIcomp );
      break;
    
    case 1: // Ternary Flory Huggins free energy
      return new Model_Energy_FHG( OpObj, Ncomp, NIcomp );
      break;

    case 2: // A + B + A-B Doi Uneyama free energy
      return new Model_Energy_DUN( OpObj, Ncomp, NIcomp );
      break;
 
    case 3: // A + B + A-B Doi Uneyama free energy
      return new Model_Energy_A_B_ABlinear( OpObj, Ncomp, NIcomp );
      break;

    case 4: // Ternary Flory Huggins free energy with Kappa3
      return new Model_Energy_FHG_K3( OpObj, Ncomp, NIcomp );
      break;
   
    default: // error catch
      std::cout << std::endl << "Error: Bad Energy key. No Energy Model initialized" << std::endl;
      return NULL;
      break;
    
  }

  return NULL;

} // }}}

Model_Mobility_Base* Model_Factory :: make_mobility_model( int key, int Ncomp, int NIcomp, Model_Visc_Base* visc_model )
{ // {{{

  switch(key) 
  {

    case 0: // Constant mobility values
      return new Model_Mobility_Const( Ncomp, NIcomp, visc_model );
      break;
    
    case 1: // Rouse dynamics in the mobility
      return new Model_Mobility_Rouse( Ncomp, NIcomp, visc_model );
      break;
    
    case 2: // Rouse dynamics scaled by the viscosity 
      return new Model_Mobility_ScaledRouse( Ncomp, NIcomp, visc_model );
      break;

    default: // error catch
      std::cout << std::endl << "Error: Bad key. No mobility model initialized" << std::endl;
      return NULL;
      break;

  }

  return NULL;

} // }}}

Model_Visc_Base* Model_Factory :: make_visc_model( int key, int Ncomp, int NIcomp )
{ // {{{

  switch(key) 
  {

    case 0: // linear viscosity
      return new Model_Visc_Lin( Ncomp, NIcomp );
      break;

    case 1: // non-linear viscosity
      return new Model_Visc_Exp( Ncomp, NIcomp );
      break;
    
    case 2: // non-linear viscosity
      return new Model_Visc_Sig( Ncomp, NIcomp );
      break;
    
    case 3: // VFTH viscosity model
      return new Model_Visc_VFTH( Ncomp, NIcomp );
      break;

    default: // error catch
      std::cout << std::endl << "Error: Bad key. No viscosity model initialized" << std::endl;
      return NULL;
      break;
    
  }

  return NULL;

} // }}}

Model_Reactions_Base* Model_Factory :: make_reactions_model( int key, int Ncomp, int NIcomp )
{ // {{{

  switch(key) 
  {

    case 0: // linear viscosity
      return new Model_Reactions_None( Ncomp, NIcomp );
      break;

    case 1: // binary reactive blending model
      return new Model_Reactions_Binary( Ncomp, NIcomp );
      break;
    
    default: // error catch
      std::cout << std::endl << "Error: Bad key. No viscosity model initialized" << std::endl;
      return NULL;
      break;
    
  }

  return NULL;

} // }}}

