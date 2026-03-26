/***************************************************************
*
* Class for Interleaved Banded Field Matrices
* 
* DRT -- Fri, 20 Mar 2017
*
****************************************************************/

#include "InterFieldMat.h"

// ----------------- Constructor/Destructor -----------------

InterFieldMat::InterFieldMat( int nrow, int ncol ) :
_Nrow(nrow),
_Ncol(ncol),
_Data( nrow, InterFieldVec(ncol) )
{ // {{{
} // }}}

InterFieldMat::InterFieldMat( const InterFieldMat &rhs ) :
_Nrow( rhs._Nrow ),
_Ncol( rhs._Ncol ),
_Data( rhs._Nrow, InterFieldVec(rhs._Ncol) )
{ // {{{

  for (UInt i=0; i<_Nrow; i++)
  {
    _Data[i]=rhs._Data[i]; 
  }

} // }}}

InterFieldMat::~InterFieldMat()
{ // {{{
} // }}}

// ----------------- Operators  -----------------

InterFieldMat & InterFieldMat::operator=( const InterFieldMat &rhs )
{ // {{{
  if (this != &rhs) // check for self-assignment
  {
    if (_Nrow != rhs._Nrow ) // check for # of elements
    {
      std::cout << "*** Error in operator=( InterFieldMat &rhs ) ***\n" << std::endl;
      std::cout << "Unequal number of elements in matrix assignment." << std::endl;
      return *this;
    }
    for (UInt i=0; i<_Nrow; ++i)
    {
      this->_Data[i] = rhs._Data[i]; // uses InterFieldVec class assignment to constant vector
    }
  }
  return *this;
} // }}}

InterFieldMat & InterFieldMat::operator=( const Vec2dFieldType &rhs )
{ //{{{
  if (this->_Nrow != rhs.size())
  {
    std::cout << "*** Error in operator=( Vec2dFieldType &rhs ) ***\n" << std::endl;
    std::cout << "Unequal number of elements in matrix assignment." << std::endl;
    return *this;
  }
  for (UInt i=0; i<_Nrow; ++i)
  {
    this->_Data[i] = rhs[i]; // uses InterFieldVec class assignment to constant vector
  }
  return *this;
} //}}}

InterFieldMat & InterFieldMat::operator=( const FieldType &rhs )
{ //{{{
  for (UInt i=0; i<this->_Nrow; ++i)
    this->_Data[i] = rhs; // uses InterFieldVec class assignment to constant scalar
  return *this;
} //}}}


InterFieldMat & InterFieldMat::operator+=( const InterFieldMat &rhs )
{ // {{{
  if (this->_Nrow != rhs._Nrow) 
  {
    std::cout << "*** Error in operator+=( InterFieldMat &rhs ) ***\n" << std::endl;
    std::cout << "Unequal number of elements in matrix += operation." << std::endl;
    return *this;
  } 
  for (UInt i=0; i<this->_Nrow; ++i) 
    this->_Data[i] += rhs._Data[i];
  return *this;
} // }}}

InterFieldMat & InterFieldMat::operator+=( const Vec2dFieldType &rhs )
{ // {{{
  if (this->_Nrow != rhs.size())
  {
    std::cout << "**(Error) in InterFieldMat::operator+=***" << std::endl;
    std::cout << "Unequal number of elements." << std::endl;
    exit(1);
  }
  for (UInt i=0; i<this->_Nrow; ++i) 
    this->_Data[i] += rhs[i];
  return *this;
} // }}}

InterFieldMat & InterFieldMat::operator+=( const FieldType &rhs )
{ // {{{
  for (UInt i=0; i<this->_Nrow; ++i)
    this->_Data[i] += rhs;
  return *this;
} // }}}


InterFieldMat & InterFieldMat::operator-=( const InterFieldMat &rhs )
{ // {{{
  if (this->_Nrow != rhs._Nrow)
  {
    std::cout << "**(Error) in InterFieldMat::operator-=***" << std::endl;
    std::cout << "Unequal number of elements." << std::endl;
    exit(1);
  }
  for (UInt i=0; i<rhs._Nrow; ++i) 
    this->_Data[i] -= rhs._Data[i];
  return *this;
} // }}}

InterFieldMat & InterFieldMat::operator-=( const Vec2dFieldType &rhs )
{ // {{{
  if (this->_Nrow != rhs.size())
  {
    std::cout << "**(Error) in InterFieldMat::operator-=***" << std::endl;
    std::cout << "Unequal number of elements." << std::endl;
    exit(1);
  }
  for (UInt i=0; i<rhs.size(); ++i) 
    this->_Data[i] -= rhs[i];
  return *this;
} // }}}

InterFieldMat & InterFieldMat::operator-=( const FieldType &rhs )
{ // {{{
  for (UInt i=0; i<this->_Nrow; ++i)
    this->_Data[i] -= rhs;
  return *this;
} // }}}


InterFieldMat & InterFieldMat::operator*=( const InterFieldMat &rhs )
{ // {{{
  if (this->_Nrow != rhs._Nrow)
  {
    std::cout << "**(Error) in InterFieldMat::operator*=***" << std::endl;
    std::cout << "Unequal number of elements." << std::endl;
    exit(1);
  }
  for (UInt i=0; i<rhs._Nrow; ++i) 
    this->_Data[i] *= rhs._Data[i];
  return *this;
} // }}}

InterFieldMat & InterFieldMat::operator*=( const Vec2dFieldType &rhs )
{ // {{{
  if (this->_Nrow != rhs.size())
  {
    std::cout << "**(Error) in InterFieldMat::operator*=***" << std::endl;
    std::cout << "Unequal number of elements." << std::endl;
    exit(1);
  }
  for (UInt i=0; i<this->_Nrow; ++i) 
    this->_Data[i] *= rhs[i];
  return *this;
} // }}}

InterFieldMat & InterFieldMat::operator*=( const FieldType &rhs )
{ // {{{
  for (UInt i=0; i<this->_Nrow; ++i) 
    this->_Data[i] *= rhs;
  return *this;
} // }}}

// ----------------- Member Functions -----------------

void InterFieldMat::zero()
{ // {{{

  for (int i = 0; i < _Nrow; i++)
  {
    _Data[i].zero();
  }

} // }}}

void InterFieldMat::setflag_inrealspace( bool realspaceflag )
{ // {{{

  for (int i = 0; i < _Nrow; i++)
  {
    _Data[i].setflag_inrealspace(realspaceflag);
  }

} // }}}

bool InterFieldMat::getflag_inrealspace() const
{ // {{{

  bool rspace_flag = true;
  bool tmp_flag = true;

  for (int i = 0; i < _Nrow; i++)
  {
    tmp_flag = _Data[i].getflag_inrealspace();
    if (i > 0 and rspace_flag != tmp_flag)
    {
      std::cout << "Warning: InterFieldVec is mixed real space/fourier space. Use component functions.";
    }
    rspace_flag = tmp_flag;
  }

  return rspace_flag;

} // }}}

void InterFieldMat::fft()
{ // {{{
  for(UInt i=0; i<this->_Nrow; ++i)
  {
    (*this)[i].fft();
  }
} // }}}

void InterFieldMat::ifft()
{ // {{{
  for(UInt i=0; i<this->_Nrow; ++i)
  {
    (*this)[i].ifft();
  }
} // }}}


void InterFieldMat::eye( FieldType scalar )
{ //{{{
  for( UInt i=0; i<_Nrow; i++)
  {
    for( UInt j=0; j<_Ncol; j++)
    {
      if ( i == j ) 
      {
        _Data[i][j] = scalar;
      }
      else
      {
        _Data[i][j].zero();
      }
    }
  }
} // }}}


void InterFieldMat::dot( const InterFieldMat &A, const InterFieldMat &B)
{ //{{{

  Field<FieldType>* tmp;
  FieldStack::FetchField( tmp );

  tmp->setflag_inrealspace( B.getflag_inrealspace() );

  if (this == &A || this == &B) // check for in-place multiplication
  {
    std::cout << "(Error) Cannot do in-place matrix multiplication" << std::endl;
    return;
  }

  if (A._Ncol != B._Nrow)
  {
    std::cout << "(Error) Cannot multiply matricies: ";
    std::cout << "A = " << A._Nrow << "x" << A._Ncol << ", ";
    std::cout << "B = " << B._Nrow << "x" << B._Ncol << std::endl;
    return;
  }
  
  this->setflag_inrealspace(B.getflag_inrealspace());
  for (UInt i=0; i<A._Nrow; ++i)
  {
    for (UInt j=0; j<B._Ncol; ++j)
    {
      (*this)[i][j] = FieldType(0.);
      for (UInt k=0; k<A._Ncol; ++k)
      {
        *tmp = A[i][k];
        *tmp *= B[k][j];
        (*this)[i][j] += *tmp;
      }
    }
  }

  FieldStack::RetireField( tmp );

} //}}}

void InterFieldMat::dot( InterFieldMat &A, InterFieldMat &B)
{ //{{{

  Field<FieldType>* tmp;
  FieldStack::FetchField( tmp );

  tmp->setflag_inrealspace( B.getflag_inrealspace() );

  if (this == &A || this == &B) // check for in-place multiplication
  {
    std::cout << "(Error) Cannot do in-place matrix multiplication" << std::endl;
    return;
  }

  if (A._Ncol != B._Nrow)
  {
    std::cout << "(Error) Cannot multiply matricies: ";
    std::cout << "A = " << A._Nrow << "x" << A._Ncol << ", ";
    std::cout << "B = " << B._Nrow << "x" << B._Ncol << std::endl;
    return;
  }

  this->setflag_inrealspace(B.getflag_inrealspace());
  for (UInt i=0; i<A._Nrow; ++i)
  {
    for (UInt j=0; j<B._Ncol; ++j)
    {
      (*this)[i][j] = FieldType(0.);
      for (UInt k=0; k<A._Ncol; ++k)
      {
        *tmp = A[i][k];
        *tmp *= B[k][j];
        (*this)[i][j] += *tmp;
      }
    }
  }

  FieldStack::RetireField( tmp );

} //}}}

void InterFieldMat::dot( const Vec2dFieldType &A, const InterFieldMat &B)
{ // {{{

  Field<FieldType>* tmp;
  FieldStack::FetchField( tmp );

  tmp->setflag_inrealspace( B.getflag_inrealspace() );

  if (this == &B) // check for in-place multiplication
  {
    std::cout << "(Error) Cannot do in-place matrix multiplication" << std::endl;
    return;
  }

  if (A[0].size() != B._Nrow)
  {
    std::cout << "(Error) Cannot multiply matricies: ";
    std::cout << "A = " << A.size() << "x" << A[0].size() << ", ";
    std::cout << "B = " << B._Nrow << "x" << B._Ncol << std::endl;
    return;
  }
  
  this->setflag_inrealspace(B.getflag_inrealspace());
  for (UInt i=0; i<A.size(); ++i)
  {
    for (UInt j=0; j<B._Ncol; ++j)
    {
      (*this)[i][j] = FieldType(0.);
      for (UInt k=0; k<B._Nrow; ++k)
      {
        *tmp = A[i][k];
        *tmp *= B[k][j];
        (*this)[i][j] += *tmp;
      }
    }
  }

  FieldStack::RetireField( tmp );

} // }}}

void InterFieldMat::dot( Vec2dFieldType &A, InterFieldMat &B)
{ // {{{

  Field<FieldType>* tmp;
  FieldStack::FetchField( tmp );

  tmp->setflag_inrealspace( B.getflag_inrealspace() );

  if (this == &B) // check for in-place multiplication
  {
    std::cout << "(Error) Cannot do in-place matrix multiplication" << std::endl;
    return;
  }

  if (A[0].size() != B._Nrow)
  {
    std::cout << "(Error) Cannot multiply matricies: ";
    std::cout << "A = " << A.size() << "x" << A[0].size() << ", ";
    std::cout << "B = " << B._Nrow << "x" << B._Ncol << std::endl;
    return;
  }
  
  this->setflag_inrealspace(B.getflag_inrealspace());
  for (UInt i=0; i<A.size(); ++i)
  {
    for (UInt j=0; j<B._Ncol; ++j)
    {
      (*this)[i][j] = FieldType(0.);
      for (UInt k=0; k<B._Nrow; ++k)
      {
        *tmp = A[i][k];
        *tmp *= B[k][j];
        (*this)[i][j] += *tmp;
      }
    }
  }

  FieldStack::RetireField( tmp );

} // }}}

void InterFieldMat::dot( const InterFieldMat &A, const Vec2dFieldType &B)
{ // {{{

  Field<FieldType>* tmp;
  FieldStack::FetchField( tmp );

  tmp->setflag_inrealspace( A.getflag_inrealspace() );

  if (this == &A) // check for in-place multiplication
  {
    std::cout << "(Error) Cannot do in-place matrix multiplication" << std::endl;
    return;
  }

  if (A._Ncol != B.size())
  {
    std::cout << "(Error) Cannot multiply matricies: ";
    std::cout << "A = " << A._Nrow << "x" << A._Ncol << ", ";
    std::cout << "B = " << B.size() << "x" << B[0].size() << std::endl;
    return;
  }
  
  this->setflag_inrealspace(A.getflag_inrealspace());
  for (UInt i=0; i<A._Nrow; ++i)
  {
    for (UInt j=0; j<B[0].size(); ++j)
    {
      (*this)[i][j] = FieldType(0.);
      for (UInt k=0; k<A._Ncol; ++k)
      {
        *tmp = A[i][k];
        *tmp *= B[k][j];
        (*this)[i][j] += *tmp;
      }
    }
  }

  FieldStack::RetireField( tmp );

} // }}}

void InterFieldMat::dot( InterFieldMat &A, Vec2dFieldType &B)
{ // {{{

  Field<FieldType>* tmp;
  FieldStack::FetchField( tmp );

  tmp->setflag_inrealspace( A.getflag_inrealspace() );

  if (this == &A) // check for in-place multiplication
  {
    std::cout << "(Error) Cannot do in-place matrix multiplication" << std::endl;
    return;
  }

  if (A._Ncol != B.size())
  {
    std::cout << "(Error) Cannot multiply matricies: ";
    std::cout << "A = " << A._Nrow << "x" << A._Ncol << ", ";
    std::cout << "B = " << B.size() << "x" << B[0].size() << std::endl;
    return;
  }
  
  this->setflag_inrealspace(A.getflag_inrealspace());
  for (UInt i=0; i<A._Nrow; ++i)
  {
    for (UInt j=0; j<B[0].size(); ++j)
    {
      (*this)[i][j] = FieldType(0.);
      for (UInt k=0; k<A._Ncol; ++k)
      {
        *tmp = A[i][k];
        *tmp *= B[k][j];
        (*this)[i][j] += *tmp;
      }
    }
  }

  FieldStack::RetireField( tmp );

} // }}}


void InterFieldMat :: invert( InterFieldMat &A )
{ // {{{
  // Invert the argument A, and make &(this) = Ainv

  if (this == &A and _Nrow != 2)
  {
    std::cout << "*** (Error) Matrix inverstion cannot be in place ***\n";
    exit(1);
  }
  
  if (this->_Nrow != this->_Ncol) // not square
  {
    std::cout << "*** (Error) Matrix inverstion only works for a square matrix ***\n";
    exit(1);
  }

  if (this->_Nrow == 1) // 2x2
  {
    (*this)[0][0] = FieldType(1.);
    (*this)[0][0] /= A[0][0];
  }
  else if (this->_Nrow == 2) // 2x2
  {
    Field<FieldType>* tmp;
    Field<FieldType>* detA;

    FieldStack::FetchField( tmp );
    FieldStack::FetchField( detA );

    // A^{-1} = 1/det(A) [[A11, -A01],[-A10, A00]
    // det(A) = A00*A11 - A01*A10
    
    *detA = A[1][1];
    *detA *= A[0][0];
    *tmp = A[0][1];
    *tmp *= A[1][0];
    *detA -= *tmp;

    // check if determinant is zero
    if ( detA->l1norm() < 1e-12 )
    {
      std::cout << "Error inverting InterFieldMatrix. The determinant of the field is zero." << std::endl; 
    }

    *tmp = (*this)[0][0];
    (*this)[0][0] = A[1][1]; (*this)[0][0] /= *detA;
    (*this)[0][1] = A[0][1]; (*this)[0][1] *= FieldType(-1); (*this)[0][1] /= *detA;
    (*this)[1][0] = A[1][0]; (*this)[1][0] *= FieldType(-1); (*this)[1][0] /= *detA;
    (*this)[1][1] = A[0][0]; (*this)[1][1] /= *detA;

    FieldStack::RetireField( tmp );
    FieldStack::RetireField( detA );

  }
  else // not 2x2, square
  {

    Field<FieldType>* tmp;
    Field<FieldType>* ratio;

    FieldStack::FetchField( tmp );
    FieldStack::FetchField( ratio );

    tmp->zero();
    ratio->zero();

    // Initialize Ainv as Eye
    for (int i = 0; i<_Nrow; i++)
    {
      for (int j = 0; j<_Ncol; j++)
      {
        if (i == j)
          (*this)[i][j] = FieldType(1.);
        else
          (*this)[i][j] = FieldType(0.);
      }
    }

    // Gauss elimination
    for (int j = 0; j<_Ncol-1; j++) // columns
    {
      for (int i = j+1; i<_Nrow; i++) // rows
      {
        *ratio = A[i][j];
        *ratio /= A[j][j]; 

        for (int k = j+1; k<_Ncol; k++) // columns
        {
          *tmp = *ratio;
          *tmp *= A[j][k];
          A[i][k] -= *tmp;
        }
      
        for (int k = 0; k<_Ncol; k++) // columns
        {
          *tmp = *ratio;
          *tmp *= (*this)[j][k];
          (*this)[i][k] -= *tmp;
        }
      }
    }

    FieldStack::RetireField(ratio);

    // back substitution
    Field<FieldType>* sum;
    FieldStack::FetchField( sum );

    for (int i=_Nrow-1; i>=0; i--)
    {
      for (int k = 0; k<_Ncol; k++)
      {
        sum->zero();
        for (int j = i+1; j<_Ncol; j++)
        {
          *tmp = (*this)[j][k];
          *tmp *= A[i][j];
          *sum += *tmp;
        }
        (*this)[i][k] -= *sum;
        (*this)[i][k] /= A[i][i];
      }
    }
    
  }

} // }}}

FieldType InterFieldMat::max()
{ // {{{

  // returns element with maximum abs value
  RealType fieldmax = 0.;
  RealType localmax = 0.;
  UInt max_i = 0;
  UInt max_j = 0;

  for( UInt i=0; i<_Nrow; i++)
  {
    for( UInt j=0; j<_Ncol; j++)
    {
      localmax = std::abs( _Data[i][j].max() );
      if (localmax > fieldmax)
      {
        fieldmax = localmax;
        max_i = i;
        max_j = j;
      }
    }
  }
  return _Data[max_i][max_j].max();

} // }}}


