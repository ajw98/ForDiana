/***************************************************************
*
* Class for Interleaved Banded Field Matrices
* 
* DRT -- Fri, 20 Mar 2017
*
****************************************************************/

#include "InterFieldVec.h"

// ----------------- Constructor/Destructor -----------------

InterFieldVec::InterFieldVec( int ncol ) :
_Ncol(ncol),
_DataPtr( ncol, NULL )
{ // {{{
  
  CreateField();
  
} // }}}

InterFieldVec::InterFieldVec( const InterFieldVec &rhs ):
_Ncol( rhs._Ncol ),
_DataPtr( _Ncol, NULL )
{ // {{{

  CreateField();
  
  for (UInt i=0; i<_Ncol; i++)
  {
    *(_DataPtr[i])=*(rhs._DataPtr[i]); 
  }

} // }}}

InterFieldVec::~InterFieldVec()
{ // {{{

  for (UInt i=0; i<_Ncol; i++)
  {
    FieldStack::RetireField( _DataPtr[i] );
  }

} // }}}

void InterFieldVec::CreateField()
{ // {{{

  for (UInt i=0; i<_Ncol; i++)
  {
    FieldStack::FetchField( _DataPtr[i] );
  }

} // }}}

// ----------------- Operators  -----------------

InterFieldVec & InterFieldVec::operator=( const InterFieldVec & rhs )
{ // {{{

  if (this != &rhs) 
  { 
    this->_Ncol = rhs.Ncol();
    for (UInt i = 0; i<_Ncol; i++)
    {
      *(_DataPtr[i])=*(rhs._DataPtr[i]); 
    }
  }
  return *this;

} // }}}

InterFieldVec & InterFieldVec::operator=( const Vec1dFieldType &rhs )
{ //{{{
  if (this->_Ncol != rhs.size())
  {
    std::cout << "**(Error) in InterFieldVec::operator=***" << std::endl;
    std::cout << "Unequal number of elements." << std::endl;
    exit(1);
  }
  for (UInt i=0; i<_Ncol; ++i)
  {
    *(_DataPtr[i]) = rhs[i];
  }
  return *this;
} //}}}

InterFieldVec & InterFieldVec::operator=( const FieldType &rhs )
{ //{{{
  for (UInt i=0; i<this->_Ncol; ++i)
  {
    *(_DataPtr[i]) = rhs; // uses field class assignment to constant operator
  }
  return *this;
} //}}}


InterFieldVec & InterFieldVec::operator+=( const InterFieldVec &rhs )
{ // {{{
  if (this->_Ncol != rhs._Ncol)
  {
    std::cout << "**(Error) in InterFieldVec::operator+=***" << std::endl;
    std::cout << "Unequal number of elements." << std::endl;
    exit(1);
  }
  for (UInt i=0; i<rhs._Ncol; ++i) 
    *(this->_DataPtr[i]) += *(rhs._DataPtr[i]);
  return *this;
} // }}}

InterFieldVec & InterFieldVec::operator+=( const Vec1dFieldType &rhs )
{ // {{{
  if (this->_Ncol != rhs.size())
  {
    std::cout << "**(Error) in InterFieldVec::operator+=***" << std::endl;
    std::cout << "Unequal number of elements." << std::endl;
    exit(1);
  }
  for (UInt i=0; i<rhs.size(); ++i) 
    *(this->_DataPtr[i]) += rhs[i];
  return *this;
} // }}}

InterFieldVec & InterFieldVec::operator+=( const FieldType &rhs )
{ // {{{
  for (UInt i=0; i<this->_Ncol; ++i)
    *(this->_DataPtr[i]) += rhs;
  return *this;
} // }}}


InterFieldVec & InterFieldVec::operator-=( const InterFieldVec &rhs )
{ // {{{
  if (this->_Ncol != rhs._Ncol)
  {
    std::cout << "**(Error) in InterFieldVec::operator-=***" << std::endl;
    std::cout << "Unequal number of elements." << std::endl;
    exit(1);
  }
  for (UInt i=0; i<rhs._Ncol; ++i) 
    *(this->_DataPtr[i]) -= *(rhs._DataPtr[i]);
  return *this;
} // }}}

InterFieldVec & InterFieldVec::operator-=( const Vec1dFieldType &rhs )
{ // {{{
  if (this->_Ncol != rhs.size())
  {
    std::cout << "**(Error) in InterFieldVec::operator-=***" << std::endl;
    std::cout << "Unequal number of elements." << std::endl;
    exit(1);
  }
  for (UInt i=0; i<rhs.size(); ++i) 
    *(this->_DataPtr[i]) += -rhs[i];
  return *this;
} // }}}

InterFieldVec & InterFieldVec::operator-=( const FieldType &rhs )
{ // {{{
  for (UInt i=0; i<this->_Ncol; ++i)
    *(this->_DataPtr[i]) += -rhs;
  return *this;
} // }}}


InterFieldVec & InterFieldVec::operator*=( const InterFieldVec &rhs )
{ // {{{
  if (this->_Ncol != rhs._Ncol)
  {
    std::cout << "**(Error) in InterFieldVec::operator*=***" << std::endl;
    std::cout << "Unequal number of elements." << std::endl;
    exit(1);
  }
  for (UInt i=0; i<rhs._Ncol; ++i) 
    *(this->_DataPtr[i]) *= *(rhs._DataPtr[i]);
  return *this;
} // }}}

InterFieldVec & InterFieldVec::operator*=( const Vec1dFieldType &rhs )
{ // {{{
  if (this->_Ncol != rhs.size())
  {
    std::cout << "**(Error) in InterFieldVec::operator*=***" << std::endl;
    std::cout << "Unequal number of elements." << std::endl;
    exit(1);
  }
  for (UInt i=0; i<rhs.size(); ++i) 
    *(this->_DataPtr[i]) *= rhs[i];
  return *this;
} // }}}

InterFieldVec & InterFieldVec::operator*=( const FieldType &rhs )
{ // {{{
  for (UInt i=0; i<this->_Ncol; ++i) 
    *(this->_DataPtr[i]) *= rhs;
  return *this;
} // }}}

// ----------------- Member Functions -----------------

void InterFieldVec::zero()
{ // {{{

  for (int i = 0; i < _Ncol; i++)
  {
    _DataPtr[i]->zero();
  }

} // }}}

void InterFieldVec::setflag_inrealspace( bool realspaceflag )
{ // {{{

  for (int i = 0; i < _Ncol; i++)
  {
    _DataPtr[i]->setflag_inrealspace(realspaceflag);
  }

} // }}}

bool InterFieldVec::getflag_inrealspace() const
{ // {{{

  bool rspace_flag = true;
  bool tmp_flag = true;

  for (int i = 0; i < _Ncol; i++)
  {
    tmp_flag = _DataPtr[i]->getflag_inrealspace();
    if (i > 0 and rspace_flag != tmp_flag)
    {
      std::cout << "Warning: InterFieldVec is mixed real space/fourier space. Use component functions.";
    }
    rspace_flag = tmp_flag;
  }

  return rspace_flag;

} // }}}

void InterFieldVec::fft()
{ // {{{
  for(UInt i=0; i<this->_Ncol; ++i)
  {
    (*this)[i].fft_rtok();
  }
} // }}}

void InterFieldVec::ifft()
{ // {{{
  for(UInt i=0; i<this->_Ncol; ++i)
  {
    (*this)[i].fft_ktor();
  }
} // }}}

// --- matrix-vector multiplication: c = A.b ---

void InterFieldVec::dot( InterFieldMat &A, InterFieldVec &b)
{ // {{{

  Field<FieldType>* tmp;
  FieldStack::FetchField( tmp );

  if ( this == &b ) // check for in-place multiplication
  {
    std::cout << "(Error) Cannot do in-place dot product" << std::endl;
    return;
  }

  if ( A.Ncol() != b.Ncol() )
  {
    std::cout << "*** (Error) Matrix multiplication of an " << A.Nrow() << "x" << A.Ncol();
    std::cout << " matrix with a vector with " << b.Ncol() << " elements." << std::endl;
    exit(1);
  }

  for (UInt i=0; i<A.Nrow(); ++i)
  {
    *(this->_DataPtr[i]) = FieldType(0.);
    for (UInt j=0; j<A.Ncol(); ++j)
    {
        *tmp = A[i][j];
        *tmp *= b[j];
        *(this->_DataPtr[i]) += *tmp;
    }
  }

  FieldStack::RetireField( tmp );

} // }}}

void InterFieldVec::dot( Vec2dFieldType &A, InterFieldVec &b)
{ // {{{

  Field<FieldType>* tmp;
  FieldStack::FetchField( tmp );

  if ( this == &b ) // check for in-place multiplication
  {
    std::cout << "(Error) Cannot do in-place dot product" << std::endl;
    return;
  }

  if ( A[0].size() != b.Ncol() )
  {
    std::cout << "*** (Error) Matrix multiplication of an " << A.size() << "x" << A[0].size();
    std::cout << " matrix with a vector with " << b.Ncol() << " elements." << std::endl;
    exit(1);
  }

  for (UInt i=0; i<A.size(); ++i)
  {
    *(this->_DataPtr[i]) = FieldType(0.);
    for (UInt j=0; j<A[0].size(); ++j)
    {
        *tmp = A[i][j];
        *tmp *= b[j];
        *(this->_DataPtr[i]) += *tmp;
    }
  }

  FieldStack::RetireField( tmp );

} // }}}

void InterFieldVec::dot( InterFieldMat &A, Vec1dFieldType &b)\
{ // {{{

  Field<FieldType>* tmp;
  FieldStack::FetchField( tmp );

  if ( A.Ncol() != b.size() )
  {
    std::cout << "*** (Error) Matrix multiplication of an " << A.Nrow() << "x" << A.Ncol();
    std::cout << " matrix with a vector with " << b.size() << " elements." << std::endl;
    exit(1);
  }

  for (UInt i=0; i<A.Nrow(); ++i)
  {
    *(this->_DataPtr[i]) = FieldType(0.);
    for (UInt j=0; j<A.Ncol(); ++j)
    {
        *tmp = A[i][j];
        *tmp *= b[j];
        *(this->_DataPtr[i]) += *tmp;
    }
  }

  FieldStack::RetireField( tmp );

} // }}}


// --- vector-matrix multiplication: c = a.B ---

void InterFieldVec::dot( InterFieldVec &a, InterFieldMat &B)
{ // {{{

  Field<FieldType>* tmp;
  FieldStack::FetchField( tmp );

  if ( this == &a ) // check for in-place multiplication
  {
    std::cout << "(Error) Cannot do in-place dot product" << std::endl;
    return;
  }
  
  for (UInt j=0; j<B.Ncol(); ++j)
  {
    *(this->_DataPtr[j]) = FieldType(0.);
    for (UInt i=0; i<B.Nrow(); ++i)
    {
        *tmp = B[i][j];
        *tmp *= a[i];
        *(this->_DataPtr[j])+= *tmp;
    }
  }

  FieldStack::RetireField( tmp );

} // }}}

void InterFieldVec::dot( Vec1dFieldType &a, InterFieldMat &B)
{ // {{{

  Field<FieldType>* tmp;
  FieldStack::FetchField( tmp );

  for (UInt j=0; j<B.Ncol(); ++j)
  {
    *(this->_DataPtr[j]) = FieldType(0.);
    for (UInt i=0; i<B.Nrow(); ++i)
    {
        *tmp = B[i][j];
        *tmp *= a[i];
        *(this->_DataPtr[j]) += *tmp;
    }
  }

  FieldStack::RetireField( tmp );

} // }}}

void InterFieldVec::dot( InterFieldVec &a, Vec2dFieldType &B)
{ // {{{

  Field<FieldType>* tmp;
  FieldStack::FetchField( tmp );

  if ( this == &a ) // check for in-place multiplication
  {
    std::cout << "(Error) Cannot do in-place dot product" << std::endl;
    return;
  }

  for (UInt j=0; j<B[0].size(); ++j)
  {
    *(this->_DataPtr[j]) = FieldType(0.);
    for (UInt i=0; i<B.size(); ++i)
    {
        *tmp = B[i][j];
        *tmp *= a[i];
        *(this->_DataPtr[j]) += *tmp;
    }
  }

  FieldStack::RetireField( tmp );

} // }}}


void InterFieldVec::linsolve( InterFieldMat &A, InterFieldVec &b )
{ // {{{

  if ( this == &b ) // check for in-place multiplication
  {
    std::cout << "(Error) Cannot do in-place linsolve" << std::endl;
    return;
  }

  if (this->_Ncol == 1)
  {
    *(this->_DataPtr[0]) = b[0];
    *(this->_DataPtr[0]) /= A[0][0];
  }
  else if (this->_Ncol == 2)
  {
    // (for 2x2): x = A^{-1} b
    // A^{-1} = 1/det(A) [[A11, -A01],[-A10, A00]
    // det(A) = A00*A11 - A01*A10
  
    Field<FieldType>* tmp;
    Field<FieldType>* detA;

    FieldStack::FetchField( tmp );
    FieldStack::FetchField( detA );

    InterFieldMat Ainv( A );

    *detA = A[1][1];
    *detA *= A[0][0];
    *tmp = A[0][1];
    *tmp *= A[1][0];
    *detA -= *tmp;

    Ainv[0][0] = A[1][1]; Ainv[0][0] /= *detA;
    Ainv[0][1] = A[0][1]; Ainv[0][1] *= -1.; Ainv[0][1] /= *detA;
    Ainv[1][0] = A[1][0]; Ainv[1][0] *= -1.; Ainv[1][0] /= *detA;
    Ainv[1][1] = A[0][0]; Ainv[1][1] /= *detA;

    this->dot(Ainv, b);

    FieldStack::RetireField( tmp );
    FieldStack::RetireField( detA );

  }
  else
  {

    Field<FieldType>* tmp;
    Field<FieldType>* ratio;

    FieldStack::FetchField( tmp );
    FieldStack::FetchField( ratio );

    tmp->zero();
    ratio->zero();

    // Gauss elimination
    for (int j = 0; j<_Ncol-1; j++) // columns
    {
      for (int i = j+1; i<_Ncol; i++) // rows
      {
        *ratio = A[i][j];
        *ratio /= A[j][j]; 
        for (int k = j+1; k<_Ncol; k++) // columns
        {
          *tmp = *ratio;
          *tmp *= A[j][k];
          A[i][k] -= *tmp;
        }
        *tmp = *ratio;
        *tmp *= b[j];
        b[i] -= *tmp;
      }
    }

    FieldStack::RetireField(ratio);

    // back substitution
    Field<FieldType>* sum;
    FieldStack::FetchField( sum );

    for (int i=_Ncol-1; i>=0; i--)
    {
      sum->zero();
      for (int j = i+1; j<_Ncol; j++)
      {
        *tmp = *(this->_DataPtr[j]);
        *tmp *= A[i][j];
        *sum += *tmp;
      }
      *(this->_DataPtr[i]) = b[i];
      *(this->_DataPtr[i]) -= *sum;
      *(this->_DataPtr[i]) /= A[i][i];
    }

    FieldStack::RetireField(tmp);
    FieldStack::RetireField(sum);
    
  }

} // }}}

