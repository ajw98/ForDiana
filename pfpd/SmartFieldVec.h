/***************************************************************
*
* Defines a class for vector operations with fields
* 
* DRT -- Wed, 30 Apr 2015
*
****************************************************************/

#ifndef _SMARTFIELDVEC_H
#define _SMARTFIELDVEC_H

#include "global.h"
#include "SmartField.h"
#include "SmartFieldMat.h"

class SmartFieldMat;

//template< typename T>
// will need to generalize this in a moment to be a template class, for now let it be
class SmartFieldVec
{

  public:

    // --- Constructor/Destructor ---
    SmartFieldVec( UInt nelem );
    SmartFieldVec( const SmartFieldVec &rhs ); // copy from another SmartFieldVec
    ~SmartFieldVec();

    // --- Operators ---
    // indexing operator
    SmartField & operator[](const UInt i) { return _Data[i]; };
    const SmartField & operator[](const UInt i) const { return _Data[i]; };
    // assignment (a = rhs)
    SmartFieldVec & operator=(const SmartFieldVec &rhs);          // other field vector
    SmartFieldVec & operator=(const Vec1dFieldType &rhs);    // const vector
    SmartFieldVec & operator=(const Vec1dReal &rhs);    // const vector
    SmartFieldVec & operator=(const Field<FieldType> &rhs);  // const scalar
    SmartFieldVec & operator=(const SmartField &rhs);  // const scalar
    SmartFieldVec & operator=(const FieldType &rhs);         // const scalar
    // in-place addition (a = a+rhs)
    SmartFieldVec & operator += (const SmartFieldVec &rhs);       // add field vector
    SmartFieldVec & operator += (const Vec1dFieldType &rhs); // add const vector
    SmartFieldVec & operator += (const Field<FieldType> &rhs);// add scalar field
    SmartFieldVec & operator += (const SmartField &rhs);// add scalar field
    SmartFieldVec & operator += (const FieldType &rhs);      // add const scalar
    // in-place subtraction (a = a-rhs)
    SmartFieldVec & operator -= (const SmartFieldVec &rhs);       // subtract field vector
    SmartFieldVec & operator -= (const Vec1dFieldType &rhs); // subtract const vector
    SmartFieldVec & operator -= (const Field<FieldType> &rhs);      // subtract const scalar
    SmartFieldVec & operator -= (const SmartField &rhs);      // subtract const scalar
    SmartFieldVec & operator -= (const FieldType &rhs);      // subtract const scalar
    // in-place elementwise/scalar multiplication (a = a*rhs)
    SmartFieldVec & operator *= (const SmartFieldVec &rhs);
    SmartFieldVec & operator *= (const Vec1dFieldType &rhs);
    SmartFieldVec & operator *= (const Field<FieldType> &rhs);
    SmartFieldVec & operator *= (const SmartField &rhs);
    SmartFieldVec & operator *= (FieldType rhs);
    SmartFieldVec & operator *= (RealType rhs);

    // --- Member functions ---
    // like a subscript operator, [], but evaluates a Smartfield vector at a vector of indices
    // ( indices of the **field** (i.e. for finite differences) not indices of the SmartFieldVec )
    void subscript(const SmartFieldVec & rhs, const std::vector<UInt> idx);
    void slice(const SmartFieldVec & rhs, const std::vector<UInt> idx);
    
    // --- Dot products ---
    // vector-vector multiply: 
    //  this = a.b (a = row vector, b = column vector)
    //  a or b can be "this" for in-place op
    // (non-member friend function)
    friend void dot( const SmartFieldVec &a, const SmartFieldVec &b, SmartField &c );
    friend void dot( const Vec1dFieldType &a, const SmartFieldVec &b, SmartField &c );
    friend void dot( const SmartFieldVec &a, const Vec1dFieldType &b, SmartField &c );

    // matrix-vector multiply:
    //  this = A.b (b = column vector)
    //  b can be "this" for in-place op
    void dot( SmartFieldMat &A, SmartFieldVec &b );
    void dot( const SmartFieldMat &A, const SmartFieldVec &b );
    void dot( Vec2dFieldType &A, SmartFieldVec &b );
    void dot( SmartFieldMat &A, Vec1dFieldType &b );

    // vector-matrix multiply:
    //  this = a.B (a = row vector)
    //  a can be "this" for in-place op
    void dot( SmartFieldVec &a, SmartFieldMat &B );
    void dot( const SmartFieldVec &a, const SmartFieldMat &B );
    void dot( Vec1dFieldType &a, SmartFieldMat &B );
    void dot( SmartFieldVec &a, Vec2dFieldType &B );

    // --- Other Routines ---
    // solve Ax=b for x
    void linsolve( SmartFieldMat &A, SmartFieldVec &b );

    // fourier transform and inverse
    void fft(bool applyscale=true);
    void ifft();

    // set/get real space flag
    void setflag_inrealspace( bool rspace_flag );
    bool getflag_inrealspace() const;

    // field manipulation/arithmetic setters
    void zero();
    void zeroreal();
    void zeroimag();
    void zero_k0(); // zero out the k=0 mode in the PS fields

    // data manipulation
    void axpby_inplace(const SmartFieldVec &y, RealType a, RealType b);
    void axpby_inplace(const SmartFieldVec &y, Vec1dReal a, Vec1dReal b);
    void xpby_inplace(const SmartFieldVec &y, RealType b);
    void xpby_inplace(const SmartFieldVec &y, Vec1dReal b);
    void axpy_inplace(const SmartFieldVec &y, RealType a);
    void axpy_inplace(const SmartFieldVec &y, Vec1dReal a);

    // maximum norm of the fields in the vector
    RealType maxl1norm();
    RealType maxl2norm();

    // --- Getters and Setters ---
    inline UInt get_Nelem() const{ return this->_Nelem; };

    friend class SmartFieldMat;

  private:

    UInt _Nelem;
    // Don't use a pointer array. It becomes a mess when trying to delete it, 
    // while still preserving the SmartField memory.
    std::vector< SmartField > _Data;
    
};

#endif //_SMARTFIELDVEC_H

