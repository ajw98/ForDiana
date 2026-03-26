/***************************************************************
*
* Class defining a wrapper around the field class for better
* use of scratch space
* 
* DRT -- Wed, 19 Aug 2015
*
****************************************************************/

#include "FieldStack.h"

FFTlayout* FieldStack::CurrLayout;
std::vector< Field<FieldType>* > FieldStack::ptr_stack;
int FieldStack::FieldCount;

// ----------------- Non-member routines for the statics -----------------
void init_FieldStack_statics( FFTlayout* curr_layout )
{ // {{{

  FieldStack::CurrLayout = curr_layout;
  FieldStack::FieldCount = 0;

  std::cout << std::endl;
  std::cout << " * Field Stack Initialized" << std::endl;
  std::cout << "   - Stack size: " << FieldStack::ptr_stack.size() << std::endl;
  std::cout << "   - Field count: " << FieldStack::FieldCount << std::endl;

} // }}}

void delete_FieldStack_statics()
{ // {{{

  Field<FieldType>* f_ptr = NULL;
  while ( ! FieldStack::ptr_stack.empty() )
  {
    f_ptr = FieldStack::ptr_stack.back(); // get field off of the stack 
    FieldStack::ptr_stack.pop_back();
    delete f_ptr;
  }
  
  std::cout << std::endl;
  std::cout << " * Field Stack Final State" << std::endl;
  std::cout << "   - Stack size: " << FieldStack::ptr_stack.size() << std::endl;
  std::cout << "   - Max # Fields : " << FieldStack::FieldCount << std::endl;

} // }}}

// ----------------- Constructor/Destructor -----------------

FieldStack :: FieldStack()
{ // {{{
} // }}}

FieldStack :: ~FieldStack()
{ // {{{
} // }}}

// ----------------- Member functions -----------------
// get some scratch space off of the stack
void FieldStack :: FetchField( Field<FieldType>* & data_ptr )
{ // {{{

  if ( ptr_stack.empty() ) // if stack is empty, put a new field on the stack
  { 
    data_ptr = new Field<FieldType>(*CurrLayout, true, "");
    FieldCount++;
    ptr_stack.push_back(data_ptr);
  }
  data_ptr = ptr_stack.back(); // get field off of the stack 
  ptr_stack.pop_back();

  // !! Important: Initialize the field (zero and put in r/k-space locally) !!

} // }}}

// put the scratch space back in the stack
void FieldStack :: RetireField( Field<FieldType>* & data_ptr )
{ // {{{

  ptr_stack.push_back(data_ptr);
 
} // }}}

