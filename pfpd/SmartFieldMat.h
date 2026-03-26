/***************************************************************
*
* Defines a class for Matrix operations on fields
* 
* DRT -- Wed, 30 Apr 2015
*
****************************************************************/

#ifndef _SMARTFIELDMAT_H
#define _SMARTFIELDMAT_H

#include "global.h"
#include "FFTlayout.h"
#include "SmartField.h"
#include "SmartFieldVec.h"

class SmartFieldVec;
void test_SmartFieldMat();

class SmartFieldMat
{

  public:

    // --- Constructor/Destructor ---
    SmartFieldMat( UInt nrow, UInt ncol );
    SmartFieldMat( const SmartFieldMat &rhs ); // copy from another SmartFieldMat
    ~SmartFieldMat();

    // --- Operators ---
    // indexing operator
    SmartField & operator()(const UInt i, const UInt j) { return _Data[i*_Ncol+j]; };
    const SmartField & operator()(const UInt i, const UInt j) const { return _Data[i*_Ncol+j]; };
    // assignment (a = rhs)
    SmartFieldMat & operator=(const SmartFieldMat &rhs);          // from other field matrix
    SmartFieldMat & operator=(const Vec2dFieldType &rhs);  // from complex 2d vector
    SmartFieldMat & operator=(const Vec2dReal &rhs);       // from real 2d vector
    SmartFieldMat & operator=(const Field<FieldType> &rhs);
    SmartFieldMat & operator=(const SmartField &rhs);
    SmartFieldMat & operator=(const FieldType &rhs);         // const scalar
    // in-place addition (a = a+rhs)
    SmartFieldMat & operator += (const SmartFieldMat &rhs);       // add field vector
    SmartFieldMat & operator += (const Vec2dFieldType &rhs); // add const vector
    SmartFieldMat & operator += (const Field<FieldType> &rhs);
    SmartFieldMat & operator += (const SmartField &rhs);
    SmartFieldMat & operator += (const FieldType &rhs);      // add const scalar
    // in-place subtraction (a = a-rhs)
    SmartFieldMat & operator -= (const SmartFieldMat &rhs);       // subtract field vector
    SmartFieldMat & operator -= (const Vec2dFieldType &rhs); // subtract const vector
    SmartFieldMat & operator -= (const Field<FieldType> &rhs);
    SmartFieldMat & operator -= (const SmartField &rhs);
    SmartFieldMat & operator -= (const FieldType &rhs);      // subtract const scalar
    // in-place elementwise/scalar multiplication (a = a*rhs)
    SmartFieldMat & operator *= (const SmartFieldMat &rhs);
    SmartFieldMat & operator *= (const Vec2dFieldType &rhs);
    SmartFieldMat & operator *= (const Field<FieldType> &rhs);
    SmartFieldMat & operator *= (const SmartField &rhs);
    SmartFieldMat & operator *= (FieldType rhs);
    SmartFieldMat & operator *= (RealType rhs);

    // --- Member Functions ---

    // like a subscript operator, [], but evaluates a Smartfield vector at a vector of indices
    // ( indices of the **field** (i.e. for finite differences) not indices of the SmartFieldVec )
    void subscript(const SmartFieldMat & rhs, const std::vector<UInt> idx);
    void slice(const SmartFieldMat & rhs, const std::vector<UInt> idx);
    
    // matrix matrix multiply C = A.B
    void dot( const SmartFieldMat &A, const SmartFieldMat &B);
    void dot( SmartFieldMat &A, SmartFieldMat &B);

    void dot( const Vec2dFieldType &A, const SmartFieldMat &B);
    void dot( Vec2dFieldType &A, SmartFieldMat &B);

    void dot( const SmartFieldMat &A, const Vec2dFieldType &B);
    void dot( SmartFieldMat &A, Vec2dFieldType &B);

    // outer product of two vectors: Cij = ai bj
    void outer( const SmartFieldVec &a, const SmartFieldVec &b );

    // Matrix inversion (makes this = A^{-1})
    void invert( SmartFieldMat &A );
    // in place matrix transpose (turns A into A^{T})
    void transpose();

    // Cholesky decomposition:
    //  On entry: (*this) contains the upper triangle of the matrix to be decomposed
    //  On exit   (*this) the lower triangle and diagonal contains the L decomposition; upper triangle is zero.
    void cholesky();

    // fourier transform and inverse
    void fft(bool applyscale=true);
    void ifft();

    // getters/setters
    void setflag_inrealspace( bool rspace_flag );
    bool getflag_inrealspace() const;
    Vec2dFieldType getelement( UInt indx ) const;

    // field manipulation/arithmetic setters
    void zero();
    void eye( FieldType scalar );

    // norm of the matrix of fields, returns a scalar field normA
    void norm ( SmartField & normA , bool AssumeSymmetric=false) const;
    // norm^2 of the matrix of fields, returns a scalar field norm2A
    void norm2 ( SmartField & norm2A , bool AssumeSymmetric=false) const;

    // --- Getters and Setters ---
    inline  UInt get_Nrow() const{ return this->_Nrow; };
    inline  UInt get_Ncol() const{ return this->_Ncol; };

    friend class SmartFieldVec;

  private:

    UInt _Nrow;
    UInt _Ncol;
    std::vector< SmartField > _Data; // the matrix is a flattened vector of fields

};

#endif //_SMARTFIELDMAT_H

