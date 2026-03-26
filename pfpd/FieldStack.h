/***************************************************************
*
* Class defining a wrapper around the field class for better
* use of scratch space
* 
* DRT -- Wed, 19 Aug 2015
*
****************************************************************/

#ifndef _FIELDSTACK_H
#define _FIELDSTACK_H

#include "global.h"
#include "Field.h"

void init_FieldStack_statics( FFTlayout * curr_layout );
void delete_FieldStack_statics();

class FieldStack
{

  public:

    // --- constructor(s) and destructor ---

    FieldStack();
    ~FieldStack();

    friend void init_FieldStack_statics( FFTlayout * curr_layout );
    friend void delete_FieldStack_statics();

    // --- static member variables ---
    // * must assign before constructor is called!!

    static FFTlayout* CurrLayout; // FFT layout of the fields
    // scratch space (stack of field pointers)
    static std::vector< Field<FieldType>* > ptr_stack; 
    static int FieldCount;

    protected:

      void FetchField( Field<FieldType>* & data_ptr );
      void RetireField( Field<FieldType>* & data_ptr );

};

#endif //_FIELDSTACK_H
