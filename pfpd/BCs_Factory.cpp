/***************************************************************
*
* Factory class for the time integration schemes
* 
* DRT -- Mon 16 May 2016
*
****************************************************************/

#include "BCs_Factory.h"

// ----------------- Constructor/Destructor -----------------

BCs_Factory :: BCs_Factory()
{ // {{{
} // }}}

BCs_Factory :: ~BCs_Factory()
{ // {{{
} // }}}

// ----------------- Factory -----------------

BCs_Base* BCs_Factory::make_BCs( int key, int phivelkey, std::string filename, 
                                 int NIcomp,  Grid* CurrGrid )
{ // {{{

  switch(key) 
  {

    case 0: // Periodic

      return new BCs_Periodic( filename, NIcomp, CurrGrid );
      break;

    case 1: // Const coeffs

      if (phivelkey == 0)
      {
        return new BCs_Const_Phi( filename, NIcomp, CurrGrid );
      }
      else
      {
        return new BCs_Const_Vel( filename, NIcomp, CurrGrid );
      }
      break;

    case 2: // Time dependent conditions 

      if (phivelkey == 0)
      {
        return new BCs_TimeDep_Phi( filename, NIcomp, CurrGrid );
      }
      else
      {
        return new BCs_TimeDep_Vel( filename, NIcomp, CurrGrid );
      }
      break;

    case 3: // VIPS boundary conditions 

      if (phivelkey == 0)
      {
        std::cout << std::endl << "You chose VIPS boundary conditions." << std::endl;
        return new BCs_VIPS( filename, NIcomp, CurrGrid );
      }
      else
      {
        std::cout << std::endl << "Error: VIPS with Model H not supported." << std::endl;
        return NULL;
      }
      break;
    
    default: // error catch
      std::cout << std::endl << "Error: Bad BCs key." << std::endl;
      return NULL;
      break;
    
  }

  return NULL;

} // }}}

