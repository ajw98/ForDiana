/***************************************************************
*
* Factory class for the time integration schemes
* 
* DRT -- Mon 16 May 2016
*
****************************************************************/

#include "TimeInt_Factory.h"

// ----------------- Constructor/Destructor -----------------

TimeInt_Factory :: TimeInt_Factory()
{ // {{{
} // }}}

TimeInt_Factory :: ~TimeInt_Factory()
{ // {{{
} // }}}

// ----------------- Factory -----------------

TimeInt_Base* TimeInt_Factory::make_timeint( int key, 
                                             int Ncomp, 
                                             int NIcomp,
                                             Grid* CurrGrid )
{ // {{{

  switch(key) 
  {

    case 0: // Model B
      return new TimeInt_ModelB( Ncomp, NIcomp, CurrGrid );
      break;

    case 1: // Model H
      return new TimeInt_ModelH( Ncomp, NIcomp, CurrGrid );
      break;

    case 2: // Model H Grid Deformation
      return new TimeInt_ModelH_GridDef( Ncomp, NIcomp, CurrGrid );
      break;

    case 3: // Model B Del_Mu formulation
      return new TimeInt_ModelB_DelMu( Ncomp, NIcomp, CurrGrid );
      break;

    case 4: // Model H Del_Mu formulation
      return new TimeInt_ModelH_DelMu( Ncomp, NIcomp, CurrGrid );
      break;
     
    case 10: // Model B with explicit stepper
      return new TimeInt_ModelB_Explicit( Ncomp, NIcomp, CurrGrid );
      break;

    default: // error catch
      std::cout << std::endl << "Error: Bad TimeInt key. No Integration method initialized" << std::endl;
      return NULL;
      break;
    
  }

  return NULL;

} // }}}

