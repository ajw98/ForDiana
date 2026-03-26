/***************************************************************
*
* Factory class for the operators
* 
* DRT -- Thurs 18 Jan 2017
*
****************************************************************/

#include "Operators_Factory.h"

// ----------------- Constructor/Destructor -----------------

Operators_Factory :: Operators_Factory()
{ // {{{
} // }}}

Operators_Factory :: ~Operators_Factory()
{ // {{{
} // }}}

// ----------------- Factory -----------------

Operators_Base* Operators_Factory :: make_operators( Grid * gridptr,
                                                     BCs_Base * bcsptr )
{ // {{{

  switch(gridptr->GridFlag()) 
  {

    case 0: // Psuedospectral operators
      return new Operators_PS(gridptr, bcsptr);
      break;

    case 1: // Hybrid Pseudospectral/Finite-Difference operators
      return new Operators_FD(gridptr, bcsptr);
      break;
    
    default: // error catch
      std::cout << std::endl << "Error: Bad Operators key. No Operators initialized" << std::endl;
      return NULL;
      break;
    
  }

  return NULL;

} // }}}

void Operators_Factory :: recycle_operators( Operators_Base * ThisOp )
{ // {{{

  delete ThisOp;

} // }}}

