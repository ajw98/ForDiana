/***************************************************************
*
* Matrix Field Class
* 
* DRT
* Wed, 30 Apr 2015
*
****************************************************************/

#include "SmartFieldMat.h"

// ------------------------------
// --- Constructor/Destructor ---
// ------------------------------

SmartFieldMat :: SmartFieldMat( UInt nrow, UInt ncol ):
//{{{
_Nrow(nrow),
_Ncol(ncol),
_Data(nrow*ncol)
{
} //}}}

SmartFieldMat :: SmartFieldMat( const SmartFieldMat &rhs ):
//{{{
_Nrow(rhs._Nrow),
_Ncol(rhs._Ncol),
_Data(rhs._Nrow*rhs._Ncol)
{
  this->setflag_inrealspace(rhs.getflag_inrealspace());
  for ( UInt i = 0; i<_Nrow; i++)
  {
    for ( UInt j = 0; j<_Ncol; j++)
    {
      (*this)(i, j) = rhs(i, j); // only copy data, not addresses
    }
  }
} //}}}

SmartFieldMat :: ~SmartFieldMat()
{ // {{{
} // }}}

// -----------------
// --- Operators ---
// -----------------

// --- Assignment ---
// assignment to a matrix field
SmartFieldMat & SmartFieldMat :: operator=(const SmartFieldMat &rhs)
{ // {{{
  if (this != &rhs) // check for self-assignment
  {
    if (_Nrow != rhs._Nrow and _Ncol != rhs._Ncol) // check for # of elements
    {
      std::cout << "Error with vector assignment: unequal number of elements." << std::endl;
      return *this;
    }

    this->setflag_inrealspace(rhs.getflag_inrealspace());
    for (UInt i=0; i<_Nrow; ++i)
    {
      for (UInt j=0; j<_Ncol; ++j)
      {
        (*this)(i,j) = rhs(i,j);
      }
    }
  }
  return *this;
} // }}}

// assignment to a constant matrix (complex)
SmartFieldMat & SmartFieldMat :: operator=(const Vec2dFieldType &rhs)
{ //{{{
  if (this->_Nrow != rhs.size() and this->_Ncol != rhs[0].size())
  {
    std::cout << "Error with vector assignment: unequal number of elements." << std::endl;
    return *this;
  }
  for (UInt i=0; i<_Nrow; ++i)
  {
    for (UInt j=0; j<_Ncol; ++j)
    {
      (*this)(i,j) = rhs[i][j]; // uses SmartFieldVec class assignment to constant vector
    }
  }
  return *this;
} //}}}

// assignment to a constant matrix (real)
SmartFieldMat & SmartFieldMat :: operator=(const Vec2dReal &rhs)
{ //{{{
  if (this->_Nrow != rhs.size() and this->_Ncol != rhs[0].size())
  {
    std::cout << "Error with vector assignment: unequal number of elements." << std::endl;
    return *this;
  }
  for (UInt i=0; i<_Nrow; ++i)
  {
    for (UInt j=0; j<_Ncol; ++j)
    {
      (*this)(i,j) = rhs[i][j]; // uses SmartFieldVec class assignment to constant vector
    }
  }
  return *this;
} //}}}

// assigment to scalar field
SmartFieldMat & SmartFieldMat :: operator=(const Field<FieldType> &rhs)
{ //{{{
  this->setflag_inrealspace(rhs.getflag_inrealspace());
  for (UInt i=0; i<this->_Nrow; ++i)
  {
    for (UInt j=0; j<this->_Ncol; ++j)
    {
      (*this)(i, j) = rhs; // uses SmartFieldVec class assignment to constant scalar
    }
  }
  return *this;
} //}}}

// assigment to scalar smart field
SmartFieldMat & SmartFieldMat :: operator=(const SmartField &rhs)
{ //{{{
  this->setflag_inrealspace(rhs.getflag_inrealspace());
  for (UInt i=0; i<this->_Nrow; ++i)
  {
    for (UInt j=0; j<this->_Ncol; ++j)
    {
      (*this)(i,j) = rhs; // uses SmartFieldVec class assignment to constant scalar
    }
  }
  return *this;
} //}}}

// assignment to a constant scalar
SmartFieldMat & SmartFieldMat :: operator=(const FieldType &rhs)
{ //{{{
  for (UInt i=0; i<this->_Nrow; ++i)
  {
    for (UInt j=0; j<this->_Ncol; ++j)
    {
      (*this)(i,j) = rhs; // uses SmartFieldVec class assignment to constant scalar
    }
  }
  return *this;
} //}}}

// --- Addition ---
// in-place addition with field matrix
SmartFieldMat & SmartFieldMat :: operator+=( const SmartFieldMat &rhs )
{ // {{{
  if (this->_Nrow != rhs._Nrow and this->_Ncol != rhs._Ncol) 
  {
    std::cout << "Error with matrix assigment: unequal number of elements." << std::endl;
    return *this;
  } 
  for (UInt i=0; i<this->_Nrow; ++i) 
    for (UInt j=0; j<this->_Ncol; ++j) 
      (*this)(i,j) += rhs(i,j);
  return *this;
} // }}}

// in-place addition with constant matrix
SmartFieldMat & SmartFieldMat :: operator+=( const Vec2dFieldType &rhs)
{ // {{{
  if (this->_Nrow != rhs.size() and this->_Ncol != rhs[0].size()) 
  {
    std::cout << "**(Error) in SmartFieldMat::operator+=***" << std::endl;
    std::cout << "Unequal number of elements." << std::endl;
    exit(1);
  }
  for (UInt i=0; i<this->_Nrow; ++i) 
    for (UInt j=0; j<this->_Ncol; ++j) 
      (*this)(i,j) += rhs[i][j];
  return *this;
} // }}}

// in-place addition with scalar field
SmartFieldMat & SmartFieldMat :: operator+=( const Field<FieldType> &rhs)
{ // {{{
  for (UInt i=0; i<this->_Nrow; ++i)
    for (UInt j=0; j<this->_Ncol; ++j) 
      (*this)(i,j) += rhs;
  return *this;
} // }}}

// in-place addition with scalar smart field
SmartFieldMat & SmartFieldMat :: operator+=( const SmartField &rhs)
{ // {{{
  for (UInt i=0; i<this->_Nrow; ++i)
    for (UInt j=0; j<this->_Ncol; ++j) 
      (*this)(i,j) += rhs;
  return *this;
} // }}}

// in-place addition with constant scalar
SmartFieldMat & SmartFieldMat :: operator+=( const FieldType &rhs)
{ // {{{
  for (UInt i=0; i<this->_Nrow; ++i)
    for (UInt j=0; j<this->_Ncol; ++j) 
      (*this)(i,j) += rhs;
  return *this;
} // }}}

// --- Subtraction ---
// in-place subtraction with field matrix
SmartFieldMat & SmartFieldMat :: operator-=( const SmartFieldMat &rhs)
{ // {{{
  if (this->_Nrow != rhs._Nrow and this->_Ncol != rhs._Ncol)
  {
    std::cout << "**(Error) in SmartFieldMat::operator-=***" << std::endl;
    std::cout << "Unequal number of elements." << std::endl;
    exit(1);
  }
  for (UInt i=0; i<this->_Nrow; ++i) 
    for (UInt j=0; j<this->_Ncol; ++j) 
      (*this)(i,j) -= rhs(i,j);
  return *this;
} // }}}

// in-place subtraction with constant matrix
SmartFieldMat & SmartFieldMat :: operator -=( const Vec2dFieldType &rhs)
{ // {{{
  if (this->_Nrow != rhs.size() and this->_Ncol != rhs[0].size())
  {
    std::cout << "**(Error) in SmartFieldMat::operator-=***" << std::endl;
    std::cout << "Unequal number of elements." << std::endl;
    exit(1);
  }
  for (UInt i=0; i<this->_Nrow; ++i) 
    for (UInt j=0; j<this->_Ncol; ++j) 
      (*this)(i,j) -= rhs[i][j];
  return *this;
} // }}}

// in-place subtraction with scalar field
SmartFieldMat & SmartFieldMat :: operator-=( const Field<FieldType> &rhs)
{ // {{{
  for (UInt i=0; i<this->_Nrow; ++i)
    for (UInt j=0; j<this->_Ncol; ++j) 
      (*this)(i,j) -= rhs;
  return *this;
} // }}}

// in-place subtraction with scalar smart field
SmartFieldMat & SmartFieldMat :: operator-=( const SmartField &rhs)
{ // {{{
  for (UInt i=0; i<this->_Nrow; ++i)
    for (UInt j=0; j<this->_Ncol; ++j) 
      (*this)(i,j) -= rhs;
  return *this;
} // }}}

// in-place subtraction with constant scalar
SmartFieldMat & SmartFieldMat :: operator-=( const FieldType &rhs)
{ // {{{
  for (UInt i=0; i<this->_Nrow; ++i)
    for (UInt j=0; j<this->_Ncol; ++j) 
      (*this)(i,j) -= rhs;
  return *this;
} // }}}

// --- Multiplication ---
// in-place elementwise multiplication by a matrix field
SmartFieldMat & SmartFieldMat :: operator*=( const SmartFieldMat &rhs)
{ // {{{
  if (this->_Nrow != rhs._Nrow and this->_Ncol != rhs._Ncol)
  {
    std::cout << "**(Error) in SmartFieldMat::operator*=***" << std::endl;
    std::cout << "Unequal number of elements." << std::endl;
    exit(1);
  }
  for (UInt i=0; i<this->_Nrow; ++i) 
    for (UInt j=0; j<this->_Ncol; ++j) 
      (*this)(i,j) *= rhs(i,j);
  return *this;
} // }}}

// in-place elementwise multiplication by a constant matrix
SmartFieldMat & SmartFieldMat :: operator*=( const Vec2dFieldType &rhs)
{ // {{{
  if (this->_Nrow != rhs.size() and this->_Ncol != rhs[0].size())
  {
    std::cout << "**(Error) in SmartFieldMat::operator*=***" << std::endl;
    std::cout << "Unequal number of elements." << std::endl;
    exit(1);
  }
  for (UInt i=0; i<this->_Nrow; ++i) 
    for (UInt j=0; j<this->_Ncol; ++j) 
      (*this)(i,j) *= rhs[i][j];
  return *this;
} // }}}

// in-place elementwise multiplication by a scalar field
SmartFieldMat & SmartFieldMat :: operator*=( const Field<FieldType> &rhs)
{ // {{{
  for (UInt i=0; i<this->_Nrow; ++i) 
    for (UInt j=0; j<this->_Ncol; ++j) 
      (*this)(i,j) *= rhs;
  return *this;
} // }}}

// in-place elementwise multiplication by a smart scalar field
SmartFieldMat & SmartFieldMat :: operator*=( const SmartField &rhs)
{ // {{{
  for (UInt i=0; i<this->_Nrow; ++i) 
    for (UInt j=0; j<this->_Ncol; ++j) 
      (*this)(i,j) *= rhs;
  return *this;
} // }}}

// in-place elementwise multiplication by a complex constant scalar
SmartFieldMat & SmartFieldMat :: operator*=(FieldType rhs)
{ // {{{
  for (UInt i=0; i<this->_Nrow; ++i) 
    for (UInt j=0; j<this->_Ncol; ++j) 
      (*this)(i,j) *= rhs;
  return *this;
} // }}}

// in-place elementwise multiplication by a real constant scalar
SmartFieldMat & SmartFieldMat :: operator*=(RealType rhs)
{ // {{{
  for (UInt i=0; i<this->_Nrow; ++i) 
    for (UInt j=0; j<this->_Ncol; ++j) 
      (*this)(i,j) *= rhs;
  return *this;
} // }}}

// ------------------------
// --- Member Functions ---
// ------------------------

// like subscript operator, [], to get a whole array (for FD)
void SmartFieldMat::subscript( const SmartFieldMat & rhs, const std::vector<UInt> idx )
{ // {{{

  if (&rhs == this)
  {
    std::cout << "***(Error in SmartFieldMat.subscript)***";
    std::cout << "Cannot do in place operation.";
    exit(1);
  }
  
  for (UInt i=0; i<_Nrow; i++)
    for (UInt j=0; j<_Ncol; j++)
      (*this)(i,j).subscript( rhs(i,j), idx );
  
} // }}}

// get a slice of an array (for FD)
void SmartFieldMat::slice( const SmartFieldMat & rhs, const std::vector<UInt> idx )
{ // {{{

  if (&rhs == this)
  {
    std::cout << "***(Error in SmartFieldMat.slice)***";
    std::cout << "Cannot do in place operation.";
    exit(1);
  }
  
  for (UInt i=0; i<_Nrow; i++)
    for (UInt j=0; j<_Ncol; j++)
      (*this)(i,j).slice( rhs(i,j), idx );

} // }}}

// --- matrix-matrix multiplication: C = A.B ---

// const field matrix/const field matrix multiply
void SmartFieldMat::dot( const SmartFieldMat &A, const SmartFieldMat &B)
{ //{{{

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
      (*this)(i,j) = A(i,0);
      (*this)(i,j) *= B(0,j);
      for (UInt k=1; k<A._Ncol; ++k)
      {
        (*this)(i,j).accumulateproduct_inplace(A(i,k), B(k,j));
      }
    }
  }
} //}}}

// field matrix/field matrix multiply
void SmartFieldMat::dot( SmartFieldMat &A, SmartFieldMat &B)
{ //{{{

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
      (*this)(i,j) = A(i,0);
      (*this)(i,j) *= B(0,j);
      for (UInt k=1; k<A._Ncol; ++k)
      {
        (*this)(i,j).accumulateproduct_inplace(A(i,k), B(k,j));
      }
    }
  }
} //}}}

// const matrix/field matrix multiply
void SmartFieldMat::dot( const Vec2dFieldType &A, const SmartFieldMat &B)
{ // {{{

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
      (*this)(i,j) = B(0,j);
      (*this)(i,j) *= A[i][0];
      for (UInt k=1; k<B._Nrow; ++k)
      {
        (*this)(i,j).xpby_inplace(B(k,j), A[i][k]);
      }
    }
  }
} // }}}

// const matrix/field matrix multiply
void SmartFieldMat::dot( Vec2dFieldType &A, SmartFieldMat &B)
{ // {{{

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
      (*this)(i,j) = B(0,j);
      (*this)(i,j) *= A[i][0];
      for (UInt k=1; k<B._Nrow; ++k)
      {
        (*this)(i,j).xpby_inplace(B(k,j), A[i][k]);
      }
    }
  }
} // }}}

// field matrix/const matrix multiply
void SmartFieldMat::dot( const SmartFieldMat &A, const Vec2dFieldType &B)
{ // {{{

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
      (*this)(i,j) = A(i,0);
      (*this)(i,j) *= B[0][j];
      for (UInt k=1; k<A._Ncol; ++k)
      {
        (*this)(i,j).xpby_inplace(A(i,k), B[k][j]);
      }
    }
  }
} // }}}

// field matrix/const matrix multiply
void SmartFieldMat::dot( SmartFieldMat &A, Vec2dFieldType &B)
{ // {{{
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
      (*this)(i,j) = A(i,0);
      (*this)(i,j) *= B[0][j];
      for (UInt k=1; k<A._Ncol; ++k)
      {
        (*this)(i,j).xpby_inplace(A(i,k), B[k][j]);
      }
    }
  }
} // }}}

// outer product of two vectors: C_{ij} = a_{i}b_{j}
void SmartFieldMat::outer( const SmartFieldVec &a, const SmartFieldVec &b )
{ //{{{
  this->setflag_inrealspace(a.getflag_inrealspace());
  for (UInt i=0; i<a._Nelem; ++i)
  {
    for (UInt j=0; j<b._Nelem; ++j)
    {
      (*this)(i,j) = a[i];
      (*this)(i,j) *= b[j];
    }
  }
} //}}}

// --- Fourier Transforms ---
void SmartFieldMat :: fft(bool applyscale /*=true*/)
{ // {{{
  for(UInt i=0; i<this->_Nrow; ++i)
  {
    for(UInt j=0; j<this->_Ncol; ++j)
    {
      (*this)(i,j).fft_rtok(applyscale);
    }
  }
} // }}}

void SmartFieldMat :: ifft()
{ // {{{
  for(UInt i=0; i<this->_Nrow; ++i)
  {
    for(UInt j=0; j<this->_Ncol; ++j)
    {
      (*this)(i,j).fft_ktor();
    }
  }
} // }}}

void SmartFieldMat :: setflag_inrealspace( bool rspace_flag )
{ // {{{
  for(UInt i=0; i<this->_Nrow; ++i)
  {
    for(UInt j=0; j<this->_Ncol; ++j)
    {
      (*this)(i,j).setflag_inrealspace( rspace_flag );
    }
  }
} // }}}

bool SmartFieldMat :: getflag_inrealspace() const
{ // {{{

  bool rspace_flag = true;
  bool tmp_flag = true;

  int k = 0;
  for(UInt i=0; i<this->_Nrow; ++i)
  {
    for(UInt j=0; j<this->_Ncol; ++j)
    {
      tmp_flag = (*this)(i,j).getflag_inrealspace();
      if (k > 0 and rspace_flag != tmp_flag)
      {
        std::cout << "Warning: SmartFieldMat is mixed real space/fourier space. Use component functions.";
      }
      rspace_flag = tmp_flag;
      k+=1;
    }
  }

  return rspace_flag;

} // }}}

Vec2dFieldType SmartFieldMat :: getelement(UInt indx) const
{ // {{{

  Vec2dFieldType element(_Nrow, Vec1dFieldType(_Ncol, FieldType(0.)));

  for ( UInt i=0; i < this->_Nrow; ++i )
    for ( UInt j=0; j < this->_Ncol; ++j )
      element[i][j] = (*this)(i,j).getelement(indx);

  return element;

} // }}}

// --- field manipulation/arithmetic setters ---
void SmartFieldMat::zero()
{ //{{{
  for( UInt i=0; i<_Nrow; i++)
    for( UInt j=0; j<_Ncol; j++)
      (*this)(i,j).zero();
} // }}}

void SmartFieldMat::eye( FieldType scalar )
{ //{{{
  for( UInt i=0; i<_Nrow; i++)
  {
    for( UInt j=0; j<_Ncol; j++)
    {
      if ( i == j ) 
      {
        (*this)(i,j) = scalar;
      }
      else
      {
        (*this)(i,j).zero();
      }
    }
  }
} // }}}

// --- Frobenius norm ---
void SmartFieldMat::norm ( SmartField & normA , bool AssumeSymmetric/*=false*/) const 
{ //{{{
  norm2(normA, AssumeSymmetric);
  normA.sqrt(normA);
} //}}}

// --- Frobenius norm ---
void SmartFieldMat::norm2 ( SmartField & norm2A , bool AssumeSymmetric/*=false*/) const 
{ //{{{

  // Use head to avoid extra init of norm2A
  norm2A  = (*this)(0,0);
  norm2A *= (*this)(0,0);

  if(AssumeSymmetric)
  {
    assert(_Ncol == _Nrow);
    // Wing
    for(UInt j=1; j<_Ncol; j++)
      norm2A.accumulateproduct_inplace((*this)(0,j), (*this)(0,j), 2.0);
    // Body with both i,j >= 1
    for (UInt j=1; j<_Ncol; ++j)
    {
      // Diagonal
      norm2A.accumulateproduct_inplace((*this)(j,j), (*this)(j,j));
      // Off diagonal (excl. wing)
      for (UInt i=j+1; i<_Nrow; ++i)
        norm2A.accumulateproduct_inplace((*this)(i,j), (*this)(i,j), 2.0);
    }
  } else {
    // Wing
    for(UInt j=1; j<_Ncol; j++)
      norm2A.accumulateproduct_inplace((*this)(0,j), (*this)(0,j));
    // Wing
    for (UInt i=1; i<_Nrow; ++i)
      norm2A.accumulateproduct_inplace((*this)(i,0), (*this)(i,0));
    // Body with both i,j >= 1
    for (UInt j=1; j<_Ncol; ++j)
      for (UInt i=1; i<_Nrow; ++i)
        norm2A.accumulateproduct_inplace((*this)(i,j), (*this)(i,j));
  }

} //}}}

// --- Matrix inversion ---
void SmartFieldMat :: invert( SmartFieldMat &A )
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

  if (this->_Nrow == 2) // 2x2
  {
    SmartField detA, invdetA;

    // A^{-1} = 1/det(A) [[A11, -A01],[-A10, A00]
    // det(A) = A00*A11 - A01*A10

    detA = A(0,0);
    detA *= A(1,1);
    detA.accumulateproduct_inplace(A(0,1), A(1,0), -1.0); // final argument switches sign. detA = A00*A11 - A01*A10
    // check if determinant is zero
    if ( detA.l1norm() < 1e-12 )
    {
      std::cout << "Error inverting SmartFieldMatrix. The determinant of the field is zero." << std::endl; 
    }
    invdetA = FieldType(1.0);
    invdetA.setflag_inrealspace(detA.getflag_inrealspace());
    invdetA /= detA;

    (*this)(0,0) = A(1,1);
    (*this)(0,0) *= invdetA;
    (*this)(1,1) = A(0,0);
    (*this)(1,1) *= invdetA;
    // Off diagonals need to be negated
    invdetA *= RealType(-1.);
    (*this)(0,1) = A(0,1);
    (*this)(0,1) *= invdetA;
    (*this)(1,0) = A(1,0);
    (*this)(1,0) *= invdetA;
  }
  else // not 2x2, square
  {

    SmartField ratio( FieldType(0.) );
    SmartField tmp( FieldType(0.) );

    // Initialize Ainv as Eye
    for (int i = 0; i<_Nrow; i++)
    {
      for (int j = 0; j<_Ncol; j++)
      {
        if (i == j)
          (*this)(i,j) = FieldType(1.);
        else
          (*this)(i,j) = FieldType(0.);
      }
    }

    // Gauss elimination
    for (int j = 0; j<_Ncol-1; j++) // columns
    {
      for (int i = j+1; i<_Nrow; i++) // rows
      {
        ratio = A(i,j);
        ratio /= A(j,j);

        for (int k = j+1; k<_Ncol; k++) // columns
        {
          A(i,k).accumulateproduct_inplace(ratio, A(j,k), -1.);
        }
      
        for (int k = 0; k<_Ncol; k++) // columns
        {
          (*this)(i,k).accumulateproduct_inplace(ratio, (*this)(j,k), -1.);
        }
      }
    }

    // back substitution
    SmartField sum; 
    for (int i=_Nrow-1; i>=0; i--)
    {
      for (int k = 0; k<_Ncol; k++)
      {
        // Unroll j=i+1 term from loop below to avoid extra init op
        if (i+1 < _Nrow) // Required since we unrolled the first loop (which doesn't then get bounds checked)
        {
          sum = (*this)(i+1,k);
          sum *= A(i,i+1);
          for (int j = i+2; j<_Ncol; j++) // Lower limit becomes i+2
          {
            sum.accumulateproduct_inplace((*this)(j,k), A(i,j));
          }
        }
        else
        {
          sum.zero();
        }
        (*this)(i,k) -= sum;
        (*this)(i,k) /= A(i,i);
      }
    }

  }

} // }}}

// --- Matrix transpose (in place) ---
void SmartFieldMat :: transpose()
{ // {{{

  SmartField tmp;

  if (this->_Nrow != this->_Ncol)
  {
    std::cout << "***(Error) Transpose function only works for a square matrix ***";
  }

  for (UInt n = 0; n < this->_Nrow; ++n)
  {
    for (UInt m = 0; m < n; ++m)
    {
      // swap row and col
      tmp = (*this)(n,m);
      (*this)(n,m) = (*this)(m,n);
      (*this)(m,n) = tmp;
    }
  }
      

} // }}}

void SmartFieldMat::cholesky()
{ // {{{
  // Note that this algorithm is O(N^3), but N is expected to be small
  if(_Nrow != _Ncol)
    codeerror_abort("Cholesky decomposition of non-square matrix not supported",__FILE__,__LINE__);

#ifdef DEBUG
  SmartFieldMat initial(*this); // Keep a copy to compare result
#endif

  UInt N(_Nrow);
  SmartFieldVec diag(N);
  SmartField sum;
  for(UInt i=0; i<N; i++)
  {
    for(UInt j=i; j<N; j++)
    {
      sum=(*this)(i,j);
      for(int k=i-1; k>=0; k--)
        sum.accumulateproduct_inplace((*this)(i,k),(*this)(j,k),-1.0);
      if(i==j)
        diag[i].sqrt(sum);
      else
      {
        (*this)(j,i) = sum;
        (*this)(j,i) /= diag[i];
      }
    }
  }
  // Now replace the diagonal of (*this) and zero the upper triangle
  for(UInt i=0; i<N; i++)
  {
    (*this)(i,i) = diag[i];
    for(UInt j=i+1; j<N; j++)
      (*this)(i,j).zero();
  }

#ifdef DEBUG
  std::cout << "DEBUG: TESTING CHOLESKY DECOMPOSITION" << std::endl << std::endl;
  // Check that M-LL^T is zero
  // Note that on entry we only required the upper triangle to be filled, so only test that
  SmartField elem;
  elem.setflag_inrealspace(initial.getflag_inrealspace());
  for(UInt i=0; i<N; i++)
  {
    for(UInt j=i; j<N; j++)
    {
      elem.zero();
      for(UInt k=0; k<N; k++)
        elem.accumulateproduct_inplace((*this)(i,k),(*this)(j,k));
      elem -= initial(i,j);
      if(!CompareEqual(elem.l2norm(),0.0))
        codeerror_abort("Error in Cholesky factorization",__FILE__,__LINE__);
    }
  }
#endif

} // }}}

void test_SmartFieldMat()
{ // {{{

  std::cout << "\n\n***Testing SmartFieldMat matrix solve***\n";

  SmartFieldMat A( 3, 3 );
  SmartFieldVec b( 3 );
  SmartFieldVec x( 3 );

  A(0,0) = FieldType(1.);
  A(0,1) = FieldType(2.);
  A(0,2) = FieldType(3.);
  A(1,0) = FieldType(4.);
  A(1,1) = FieldType(5.);
  A(1,2) = FieldType(4.);
  A(2,0) = FieldType(3.);
  A(2,1) = FieldType(2.);
  A(2,2) = FieldType(1.);

  b[0] = FieldType(14.);
  b[1] = FieldType(26.);
  b[2] = FieldType(10.);

  x.linsolve(A, b);

  std::cout << "x[0] = " << x[0].getelement(0) << std::endl; // x[0] = 1
  std::cout << "x[1] = " << x[1].getelement(0) << std::endl; // x[1] = 2
  std::cout << "x[2] = " << x[2].getelement(0) << std::endl; // x[2] = 3

  A(0,0) = FieldType(1.);
  A(0,1) = FieldType(2.);
  A(0,2) = FieldType(3.);
  A(1,0) = FieldType(4.);
  A(1,1) = FieldType(5.);
  A(1,2) = FieldType(4.);
  A(2,0) = FieldType(3.);
  A(2,1) = FieldType(2.);
  A(2,2) = FieldType(1.);

  b[0] = FieldType(14.);
  b[1] = FieldType(26.);
  b[2] = FieldType(10.);

  SmartFieldMat Ainv( 3, 3 );
  Ainv.invert( A );

  std::cout << "\nA inverse\n";
  std::cout << Ainv(0,0).getelement(0) << "\t";
  std::cout << Ainv(0,1).getelement(0) << "\t";
  std::cout << Ainv(0,2).getelement(0) << std::endl; 
  std::cout << Ainv(1,0).getelement(0) << "\t";
  std::cout << Ainv(1,1).getelement(0) << "\t";
  std::cout << Ainv(1,2).getelement(0) << std::endl; 
  std::cout << Ainv(2,0).getelement(0) << "\t";
  std::cout << Ainv(2,1).getelement(0) << "\t";
  std::cout << Ainv(2,2).getelement(0) << std::endl;

  x.dot(Ainv, b);

  std::cout << "x[0] = " << x[0].getelement(0) << std::endl;
  std::cout << "x[1] = " << x[1].getelement(0) << std::endl;
  std::cout << "x[2] = " << x[2].getelement(0) << std::endl;

  std::cout << "\n***Finished SmartFieldMat Test***\n";

} // }}}

