/***************************************************************
*
* Class member functions for vector operations on fields
* 
* DRT
* Wed, 30 Apr 2015
*
****************************************************************/

#include "SmartFieldVec.h"

// ------------------------------
// --- Constructor/Destructor ---
// ------------------------------

SmartFieldVec :: SmartFieldVec( UInt nelem ) : 
//{{{
_Nelem(nelem),
_Data(nelem)
{
} //}}}

SmartFieldVec :: SmartFieldVec( const SmartFieldVec &rhs ) : 
//{{{
_Nelem(rhs._Nelem),
_Data(rhs._Nelem)
{ 
  for ( UInt i = 0; i < _Nelem; i++)
  {
    _Data[i] = rhs._Data[i]; // only copy data, not address!
  }
} //}}}

SmartFieldVec :: ~SmartFieldVec()
{ // {{{
} // }}}

// -----------------
// --- Operators ---
// -----------------

// --- Assignment ---
// assignment to a vector field
SmartFieldVec & SmartFieldVec :: operator=(const SmartFieldVec &rhs)
{ // {{{
  if (this != &rhs) // check for self-assignment
  {
    if (_Nelem != rhs._Nelem ) // check for # of elements
    {
      std::cout << "**(Error) in SmartFieldVec::operator=***" << std::endl;
      std::cout << "Unequal number of elements." << std::endl;
      exit(1);
    }
    setflag_inrealspace(rhs.getflag_inrealspace());
    for ( UInt i = 0; i < _Nelem; i++)
    {
      this->_Data[i] = rhs._Data[i];
    }
  }
  return *this;
} // }}}

// assignment to a constant vector (complex)
SmartFieldVec & SmartFieldVec :: operator=(const Vec1dFieldType &rhs)
{ //{{{
  if (this->_Nelem != rhs.size())
  {
    std::cout << "**(Error) in SmartFieldVec::operator=***" << std::endl;
    std::cout << "Unequal number of elements." << std::endl;
    exit(1);
  }
  for (UInt i=0; i<_Nelem; ++i)
  {
    this->_Data[i] = rhs[i]; // uses field class assignment to constant operator
  }
  return *this;
} //}}}

// assignment to a constant vector (real)
SmartFieldVec & SmartFieldVec :: operator=(const Vec1dReal &rhs)
{ //{{{
  if (this->_Nelem != rhs.size())
  {
    std::cout << "**(Error) in SmartFieldVec::operator=***" << std::endl;
    std::cout << "Unequal number of elements." << std::endl;
    exit(1);
  }
  for (UInt i=0; i<_Nelem; ++i)
  {
    this->_Data[i] = rhs[i]; // uses field class assignment to constant operator
  }
  return *this;
} //}}}

/// assignment to a scalar field
SmartFieldVec & SmartFieldVec :: operator=(const Field<FieldType> &rhs)
{ //{{{
  setflag_inrealspace(rhs.getflag_inrealspace());
  for (UInt i=0; i<this->_Nelem; ++i)
  {
    this->_Data[i] = rhs; // uses field class assignment to constant operator
  }
  return *this;
} //}}}

/// assignment to a scalar field
SmartFieldVec & SmartFieldVec :: operator=(const SmartField &rhs)
{ //{{{
  for (UInt i=0; i<this->_Nelem; ++i)
  {
    this->_Data[i] = rhs; // uses field class assignment to constant operator
  }
  return *this;
} //}}}

// assignment to a constant scalar
SmartFieldVec & SmartFieldVec :: operator=(const FieldType &rhs)
{ //{{{
  for (UInt i=0; i<this->_Nelem; ++i)
  {
    this->_Data[i] = rhs; // uses field class assignment to constant operator
  }
  return *this;
} //}}}


// --- Addition ---
// in-place addition with vector field
SmartFieldVec & SmartFieldVec :: operator+=( const SmartFieldVec &rhs )
{ // {{{
  if (this->_Nelem != rhs._Nelem)
  {
    std::cout << "**(Error) in SmartFieldVec::operator+=***" << std::endl;
    std::cout << "Unequal number of elements." << std::endl;
    exit(1);
  }
  for (UInt i=0; i<rhs._Nelem; ++i) 
    this->_Data[i] += rhs._Data[i];
  return *this;
} // }}}

// in-place addition with constant vector
SmartFieldVec & SmartFieldVec :: operator+=( const Vec1dFieldType &rhs)
{ // {{{
  if (this->_Nelem != rhs.size())
  {
    std::cout << "**(Error) in SmartFieldVec::operator+=***" << std::endl;
    std::cout << "Unequal number of elements." << std::endl;
    exit(1);
  }
  for (UInt i=0; i<rhs.size(); ++i) 
    this->_Data[i] += rhs[i];
  return *this;
} // }}}

// in-place addition with scalar field
SmartFieldVec & SmartFieldVec :: operator+=( const Field<FieldType> &rhs)
{ // {{{
  for (UInt i=0; i<this->_Nelem; ++i)
    this->_Data[i] += rhs;
  return *this;
} // }}}

// in-place addition with scalar field
SmartFieldVec & SmartFieldVec :: operator+=( const SmartField &rhs)
{ // {{{
  for (UInt i=0; i<this->_Nelem; ++i)
    this->_Data[i] += rhs;
  return *this;
} // }}}

// in-place addition with constant scalar
SmartFieldVec & SmartFieldVec :: operator+=( const FieldType &rhs)
{ // {{{
  for (UInt i=0; i<this->_Nelem; ++i)
    this->_Data[i] += rhs;
  return *this;
} // }}}


// --- Subtraction ---
// in-place subtraction with vector field
SmartFieldVec & SmartFieldVec :: operator-=( const SmartFieldVec &rhs)
{ // {{{
  if (this->_Nelem != rhs._Nelem)
  {
    std::cout << "**(Error) in SmartFieldVec::operator-=***" << std::endl;
    std::cout << "Unequal number of elements." << std::endl;
    exit(1);
  }
  for (UInt i=0; i<rhs._Nelem; ++i) 
    this->_Data[i] -= rhs._Data[i];
  return *this;
} // }}}

// in-place subtraction with constant vector
SmartFieldVec & SmartFieldVec :: operator-=( const Vec1dFieldType &rhs)
{ // {{{
  if (this->_Nelem != rhs.size())
  {
    std::cout << "**(Error) in SmartFieldVec::operator-=***" << std::endl;
    std::cout << "Unequal number of elements." << std::endl;
    exit(1);
  }
  for (UInt i=0; i<rhs.size(); ++i) 
    this->_Data[i] += -rhs[i];
  return *this;
} // }}}

// in-place subtraction with scalar field
SmartFieldVec & SmartFieldVec :: operator-=( const SmartField &rhs)
{ // {{{
  for (UInt i=0; i<this->_Nelem; ++i)
    this->_Data[i] -= rhs;
  return *this;
} // }}}

// in-place subtraction with scalar field
SmartFieldVec & SmartFieldVec :: operator-=( const Field<FieldType> &rhs)
{ // {{{
  for (UInt i=0; i<this->_Nelem; ++i)
    this->_Data[i] -= rhs;
  return *this;
} // }}}

// in-place subtraction with constant scalar
SmartFieldVec & SmartFieldVec :: operator-=( const FieldType &rhs)
{ // {{{
  for (UInt i=0; i<this->_Nelem; ++i)
    this->_Data[i] += -rhs;
  return *this;
} // }}}


// --- Multiplication ---
// in-place elementwise multiplication by a vector field
SmartFieldVec & SmartFieldVec :: operator*=( const SmartFieldVec &rhs)
{ // {{{
  if (this->_Nelem != rhs._Nelem)
  {
    std::cout << "**(Error) in SmartFieldVec::operator*=***" << std::endl;
    std::cout << "Unequal number of elements." << std::endl;
    exit(1);
  }
  for (UInt i=0; i<rhs._Nelem; ++i) 
    this->_Data[i] *= rhs._Data[i];
  return *this;
} // }}}

// in-place elementwise multiplication by a constant vector
SmartFieldVec & SmartFieldVec :: operator*=( const Vec1dFieldType &rhs)
{ // {{{
  if (this->_Nelem != rhs.size())
  {
    std::cout << "**(Error) in SmartFieldVec::operator*=***" << std::endl;
    std::cout << "Unequal number of elements." << std::endl;
    exit(1);
  }
  for (UInt i=0; i<rhs.size(); ++i) 
    this->_Data[i] *= rhs[i];
  return *this;
} // }}}

// in-place elementwise multiplication by a scalar field
SmartFieldVec & SmartFieldVec :: operator*=( const Field<FieldType> &rhs)
{ // {{{
  for (UInt i=0; i<this->_Nelem; ++i) 
    this->_Data[i] *= rhs;
  return *this;
} // }}}

// in-place elementwise multiplication by a scalar field
SmartFieldVec & SmartFieldVec :: operator*=( const SmartField &rhs)
{ // {{{
  for (UInt i=0; i<this->_Nelem; ++i) 
    this->_Data[i] *= rhs;
  return *this;
} // }}}

// in-place elementwise multiplication by a complex constant scalar
SmartFieldVec & SmartFieldVec :: operator*=(FieldType rhs)
{ // {{{
  for (UInt i=0; i<this->_Nelem; ++i) 
    this->_Data[i] *= rhs;
  return *this;
} // }}}

// in-place elementwise multiplication by a real constant scalar
SmartFieldVec & SmartFieldVec :: operator*=(RealType rhs)
{ // {{{
  for (UInt i=0; i<this->_Nelem; ++i) 
    this->_Data[i] *= rhs;
  return *this;
} // }}}

// ------------------------
// --- Member Functions ---
// ------------------------

// like subscript operator, [], to get a whole array (for FD)
void SmartFieldVec::subscript( const SmartFieldVec & rhs, const std::vector<UInt> idx )
{ // {{{

  if (&rhs == this)
  {
    std::cout << "***(Error in SmartFieldVec.subscript)***";
    std::cout << "Cannot do in place operation.";
    exit(1);
  }
  
  for (UInt i=0; i<_Nelem; i++)
    _Data[i].subscript( rhs[i], idx );
  
} // }}}

// get a slice of an array (for FD)
void SmartFieldVec::slice( const SmartFieldVec & rhs, const std::vector<UInt> idx )
{ // {{{

  if (&rhs == this)
  {
    std::cout << "***(Error in SmartFieldVec.slice)***";
    std::cout << "Cannot do in place operation.";
    exit(1);
  }
  
  for (UInt i=0; i<_Nelem; i++)
    _Data[i].slice( rhs[i], idx );
  
} // }}}

// -------------------
// --- Dot Product ---
// -------------------

// --- vector-vector multiply ---
// field vector/field vector multiply: c = a.b (non-member friend)
void dot( const SmartFieldVec &a, const SmartFieldVec &b, SmartField &c)
{ // {{{

  // Partial unroll
  c = a._Data[0];
  c *= b._Data[0];
  for (UInt i=1; i<a._Nelem; ++i)
  {
    c.accumulateproduct_inplace(a[i], b[i]);
  }

} // }}}

// const vector/field matrix multiply: c = a.b (non-member friend)
void dot( const Vec1dFieldType &a, const SmartFieldVec &b, SmartField &c)
{ // {{{

  // Partial unroll
  c = b._Data[0];
  c *= a[0];
  for (UInt i=1; i<a.size(); ++i) 
  {
    c.xpby_inplace(b[i], a[i]);
  }

} // }}}

// field vector/const matrix multiply: c = a.b (non-member friend)
void dot( const SmartFieldVec &a, const Vec1dFieldType &b, SmartField &c)
{ // {{{

  c = a._Data[0];
  c *= b[0];
  for (UInt i=1; i<a._Nelem; ++i) 
  {
    c.xpby_inplace(a[i], b[i]);
  }

} // }}}

// --- matrix-vector multiplication: c = A.b ---
// field matrix/field vector multiply
void SmartFieldVec :: dot( SmartFieldMat &A, SmartFieldVec &b)
{ // {{{

  bool realspaceflag = this->getflag_inrealspace();

  if ( A.getflag_inrealspace() != realspaceflag or
       b.getflag_inrealspace() != realspaceflag   )
  {
    std::cout << "*** (Error) in SmartFieldVec::dot ***"<< std::endl;
    std::cout << "Matrix/vector multiplication not in compatible Fourier representation\n";
    exit(1);
  }

  if ( this == &b ) // check for in-place multiplication
  {
    std::cout << "(Error) Cannot do in-place dot product" << std::endl;
    return;
  }

  if ( A._Ncol != b.get_Nelem() )
  {
    std::cout << "*** (Error) Matrix multiplication of an " << A._Nrow << "x" << A._Ncol << " matrix with a vector with " << b.get_Nelem() << " elements." << std::endl;
    exit(1);
  }

  for (UInt i=0; i<A._Nrow; ++i)
  {
    this->_Data[i] = A(i,0);
    this->_Data[i] *= b[0];
    for (UInt j=1; j<A._Ncol; ++j)
    {
      this->_Data[i].accumulateproduct_inplace(A(i,j), b[j]);
    }
  }

} // }}}

// field matrix/field vector multiply
void SmartFieldVec :: dot( const SmartFieldMat &A, const SmartFieldVec &b)
{ // {{{

  bool realspaceflag = this->getflag_inrealspace();

  if ( A.getflag_inrealspace() != realspaceflag or
       b.getflag_inrealspace() != realspaceflag   )
  {
    std::cout << "*** (Error) in SmartFieldVec::dot ***"<< std::endl;
    std::cout << "Matrix/vector multiplication not in compatible Fourier representation\n";
    exit(1);
  }

  if ( this == &b ) // check for in-place multiplication
  {
    std::cout << "(Error) Cannot do in-place dot product" << std::endl;
    return;
  }

  if ( A._Ncol != b.get_Nelem() )
  {
    std::cout << "*** (Error) Matrix multiplication of an " << A._Nrow << "x" << A._Ncol << " matrix with a vector with " << b.get_Nelem() << " elements." << std::endl;
    exit(1);
  }

  for (UInt i=0; i<A._Nrow; ++i)
  {
    this->_Data[i] = A(i,0);
    this->_Data[i] *= b[0];
    for (UInt j=1; j<A._Ncol; ++j)
    {
        this->_Data[i].accumulateproduct_inplace(A(i,j), b[j]);
    }
  }
} // }}}

// const matrix/field vector multiply
void SmartFieldVec :: dot( Vec2dFieldType &A, SmartFieldVec &b)
{ // {{{

  if ( this == &b ) // check for in-place multiplication
  {
    std::cout << "(Error) Cannot do in-place dot product" << std::endl;
    return;
  }

  if ( A[0].size() != b.get_Nelem() )
  {
    std::cout << "*** (Error) Matrix multiplication of an " << A.size() << "x" << A[0].size() << " matrix with a vector with " << b.get_Nelem() << " elements." << std::endl;
    exit(1);
  }

  for (UInt i=0; i<A.size(); ++i)
  {
    this->_Data[i] = b[0];
    this->_Data[i] *= A[i][0];
    for (UInt j=1; j<A[0].size(); ++j)
    {
        this->_Data[i].xpby_inplace(b[j], A[i][j]);
    }
  }

} // }}}

// field matrix/const vector multiply
void SmartFieldVec :: dot( SmartFieldMat &A, Vec1dFieldType &b)
{ // {{{

  if ( A._Ncol != b.size() )
  {
    std::cout << "*** (Error) Matrix multiplication of an " << A._Nrow << "x" << A._Ncol << " matrix with a vector with " << b.size() << " elements." << std::endl;
    exit(1);
  }

  for (UInt i=0; i<A._Nrow; ++i)
  {
    this->_Data[i] = A(i,0);
    this->_Data[i] *= b[0];
    for (UInt j=1; j<A._Ncol; ++j)
    {
        this->_Data[i].xpby_inplace(A(i,j), b[j]);
    }
  }

} // }}}

// --- vector-matrix multiplication: c = a.B ---

// field vector/field matrix multiply
void SmartFieldVec :: dot( SmartFieldVec &a, SmartFieldMat &B)
{ // {{{

  if ( this == &a ) // check for in-place multiplication
  {
    std::cout << "(Error) Cannot do in-place dot product" << std::endl;
    return;
  }

  for (UInt j=0; j<B._Ncol; ++j)
  {
    this->_Data[j] = B(0,j);
    this->_Data[j] *= a[0];
    for (UInt i=1; i<B._Nrow; ++i)
    {
        this->_Data[j].accumulateproduct_inplace(B(i,j),a[i]);
    }
  }

} // }}}

// field vector/field matrix multiply
void SmartFieldVec :: dot( const SmartFieldVec &a, const SmartFieldMat &B)
{ // {{{

  if ( this == &a ) // check for in-place multiplication
  {
    std::cout << "(Error) Cannot do in-place dot product" << std::endl;
    return;
  }

  for (UInt j=0; j<B._Ncol; ++j)
  {
    this->_Data[j] = B(0,j);
    this->_Data[j] *= a[0];
    for (UInt i=1; i<B._Nrow; ++i)
    {
        this->_Data[j].accumulateproduct_inplace(B(i,j),a[i]);
    }
  }

} // }}}

// const vector/field matrix multiply
void SmartFieldVec :: dot( Vec1dFieldType &a, SmartFieldMat &B)
{ // {{{

  for (UInt j=0; j<B._Ncol; ++j)
  {
    this->_Data[j] = B(0,j);
    this->_Data[j] *= a[0];
    for (UInt i=1; i<B._Nrow; ++i)
    {
      this->_Data[j].xpby_inplace(B(i,j), a[i]);
    }
  }

} // }}}

// field vector/const matrix multiply
void SmartFieldVec :: dot( SmartFieldVec &a, Vec2dFieldType &B)
{ // {{{

  if ( this == &a ) // check for in-place multiplication
  {
    std::cout << "(Error) Cannot do in-place dot product" << std::endl;
    return;
  }

  for (UInt j=0; j<B[0].size(); ++j)
  {
    this->_Data[j] = a[0];
    this->_Data[j] *= B[0][j];
    for (UInt i=1; i<B.size(); ++i)
    {
      this->_Data[j].xpby_inplace(a[i], B[i][j]);
    }
  }

} // }}}

// --------------------------------------
// ---        Other Routines          ---
// --------------------------------------

// solve Ax=b for x
void SmartFieldVec :: linsolve( SmartFieldMat &A, SmartFieldVec &b )
{ // {{{

  if ( this == &b ) // check for in-place multiplication
  {
    std::cout << "(Error) Cannot do in-place linsolve" << std::endl;
    return;
  }

  if (this->_Nelem == 2)
  {
    // (for 2x2): x = A^{-1} b
    // A^{-1} = 1/det(A) [[A11, -A01],[-A10, A00]
    // det(A) = A00*A11 - A01*A10

    SmartField detA, invdetA;

    detA = A(0,0);
    detA *= A(1,1);
    detA.accumulateproduct_inplace(A(0,1), A(1,0), -1.0); // final argument switches sign. detA = A00*A11 - A01*A10
    invdetA = FieldType(1.0);
    invdetA.setflag_inrealspace(detA.getflag_inrealspace());
    invdetA /= detA;

    (*this)[0] = b[0];
    (*this)[0] *= A(1,1);
    (*this)[0].accumulateproduct_inplace(A(0,1), b[1], -1.0);
    (*this)[0] *= invdetA;
    
    (*this)[1] = b[1];
    (*this)[1] *= A(0,0);
    (*this)[1].accumulateproduct_inplace(A(1,0), b[0], -1.0);
    (*this)[1] *= invdetA;

  }
  else
  {
    SmartField ratio( FieldType(0.) );

    // Gauss elimination
    for (int j = 0; j<_Nelem-1; j++) // columns
    {
      for (int i = j+1; i<_Nelem; i++) // rows
      {
        ratio = A(i,j);
        ratio /= A(j,j);
        for (int k = j+1; k<_Nelem; k++) // columns
        {
          A(i,k).accumulateproduct_inplace(ratio, A(j,k),-1.);
        }
        b[i].accumulateproduct_inplace(b[j], ratio, -1.0);
      }
    }

    // back substitution
    SmartField sum;
    for (int i=_Nelem-1; i>=0; i--)
    {
      if(i+1 < _Nelem) // Required since we unrolled the first loop (which doesn't then get bounds checked)
      {
        // Unroll j=i+1 term from loop below to avoid extra init op
        sum = this->_Data[i+1];
        sum *= A(i,i+1);
        for (int j = i+2; j<_Nelem; j++) // Lower limit becomes i+2
        {
          sum.accumulateproduct_inplace(this->_Data[j], A(i,j));
        }
      } 
      else 
      {
        sum.zero();
	sum.setflag_inrealspace(b[i].getflag_inrealspace());
      }
      this->_Data[i] = b[i];
      this->_Data[i] -= sum;
      this->_Data[i] /= A(i,i);
    }
  }

} // }}}

void SmartFieldVec :: fft(bool applyscale /*=true*/)
{ // {{{
  for(UInt i=0; i<this->_Nelem; ++i)
  {
    (*this)[i].fft_rtok(applyscale);
  }
} // }}}

void SmartFieldVec :: ifft()
{ // {{{
  for(UInt i=0; i<this->_Nelem; ++i)
  {
    (*this)[i].fft_ktor();
  }
} // }}}

void SmartFieldVec :: setflag_inrealspace( bool rspace_flag )
{ // {{{
  for(UInt i=0; i<this->_Nelem; ++i)
  {
    (*this)[i].setflag_inrealspace(rspace_flag);
  }
} // }}}

bool SmartFieldVec :: getflag_inrealspace() const
{ // {{{

  bool rspace_flag = true;
  bool tmp_flag = true;

  for(UInt i=0; i<this->_Nelem; ++i)
  {
    tmp_flag = this->_Data[i].getflag_inrealspace();
    if (i > 0 and rspace_flag != tmp_flag)
    {
      std::cout << "Warning: SmartFieldVec is mixed real space/fourier space. Use component functions.";
    }
    rspace_flag = tmp_flag;
  }

  return rspace_flag;

} // }}}

// field manipulation/arithmetic setters
void SmartFieldVec::zero()
{ //{{{
  for( UInt m=0; m<_Data.size(); m++)
  {
    _Data[m].zero();
  }
} // }}}

void SmartFieldVec::zeroreal()
{ //{{{
  for( UInt m=0; m<_Data.size(); m++)
  {
    _Data[m].zeroreal();
  }
} // }}}

void SmartFieldVec::zeroimag()
{ //{{{
  for( UInt m=0; m<_Data.size(); m++)
  {
    _Data[m].zeroimag();
  }
} // }}}

void SmartFieldVec::zero_k0()
{ //{{{
  for( UInt m=0; m<_Data.size(); m++)
  {
    _Data[m].zero_k0();
  }
} // }}}

void SmartFieldVec :: axpby_inplace(const SmartFieldVec &y, RealType a, RealType b)
{ // {{{
  for( UInt i=0; i<_Nelem; i++)
  {
    _Data[i].axpby_inplace(y[i], a, b);
  }
} // }}}

void SmartFieldVec :: axpby_inplace(const SmartFieldVec &y, Vec1dReal a, Vec1dReal b)
{ // {{{
  for( UInt i=0; i<_Nelem; i++)
  {
    _Data[i].axpby_inplace(y[i], a, b);
  }
} // }}}

void SmartFieldVec :: xpby_inplace(const SmartFieldVec &y, RealType b)
{ // {{{
  for( UInt i=0; i<_Nelem; i++)
  {
    _Data[i].xpby_inplace(y[i], b);
  }
} // }}}

void SmartFieldVec :: xpby_inplace(const SmartFieldVec &y, Vec1dReal b)
{ // {{{
  for( UInt i=0; i<_Nelem; i++)
  {
    _Data[i].xpby_inplace(y[i], b);
  }
} // }}}

void SmartFieldVec :: axpy_inplace(const SmartFieldVec &y, RealType a)
{ // {{{
  for( UInt i=0; i<_Nelem; i++)
  {
    _Data[i].axpy_inplace(y[i], a);
  }
} // }}}

void SmartFieldVec :: axpy_inplace(const SmartFieldVec &y, Vec1dReal a)
{ // {{{
  for( UInt i=0; i<_Nelem; i++)
  {
    _Data[i].axpy_inplace(y[i], a);
  }
} // }}}

RealType SmartFieldVec :: maxl2norm()
{ // {{{

  RealType this_max;

  this_max = 0.;
  for(UInt i=0; i<this->_Nelem; ++i)
  {
    this_max = std::max( (*this)[i].l2norm(), this_max );
  }

  return this_max;
} // }}}

RealType SmartFieldVec :: maxl1norm()
{ // {{{

  RealType this_max;

  this_max = 0.;
  for(UInt i=0; i<this->_Nelem; ++i)
  {
    this_max = std::max( (*this)[i].l1norm(), this_max );
  }

  return this_max;
} // }}}

