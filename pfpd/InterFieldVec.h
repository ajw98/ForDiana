/***************************************************************
*
* Class for the psuedospectral operators
* 
* DRT -- Fri, 13 Mar 2015
*
****************************************************************/

#ifndef _INTERFIELDVEC_H
#define _INTERFIELDVEC_H

#include "global.h"
#include "FieldStack.h"
#include "InterFieldMat.h"

class InterFieldMat;

class InterFieldVec : public FieldStack
{

  public:

    // --- Constructor/Destructor ---

    InterFieldVec( int ncol );
    InterFieldVec( const InterFieldVec &rhs ); // copy constructor
    ~InterFieldVec();

    void CreateField();

    // --- operators (field algebra) ---

    // indexing operator
    Field<FieldType> & operator[](const UInt i) { return *(_DataPtr[i]); };
    const Field<FieldType> & operator[](const UInt i) const { return *(_DataPtr[i]); };

    // assignment (a = rhs)
    InterFieldVec & operator=(const InterFieldVec &rhs);  // other InterFieldVec
    InterFieldVec & operator=(const Vec1dFieldType &rhs); // const vector
    InterFieldVec & operator=(const FieldType &rhs);      // const scalar

    // in-place addition (a = a+rhs)
    InterFieldVec & operator += (const InterFieldVec &rhs);       // add field vector
    InterFieldVec & operator += (const Vec1dFieldType &rhs); // add const vector
    InterFieldVec & operator += (const FieldType &rhs);      // add const scalar

    // in-place subtraction (a = a-rhs)
    InterFieldVec & operator -= (const InterFieldVec &rhs);       // subtract field vector
    InterFieldVec & operator -= (const Vec1dFieldType &rhs); // subtract const vector
    InterFieldVec & operator -= (const FieldType &rhs);      // subtract const scalar

    // in-place elementwise/scalar multiplication (a = a*rhs)
    InterFieldVec & operator *= (const InterFieldVec &rhs);
    InterFieldVec & operator *= (const Vec1dFieldType &rhs);
    InterFieldVec & operator *= (const FieldType &rhs);

    // --- Member Functions ---

    // getters/setters
    inline UInt Ncol() const { return _Ncol; };

    void zero();
    void setflag_inrealspace( bool realspaceflag );
    bool getflag_inrealspace() const;

    // fourier transform and inverse
    void fft();
    void ifft();

    // matrix-vector multiply:
    //  this = A.b (b = column vector, not in-place)
    void dot( InterFieldMat &A, InterFieldVec &b );
    void dot( Vec2dFieldType &A, InterFieldVec &b );
    void dot( InterFieldMat &A, Vec1dFieldType &b );

    // vector-matrix multiply:
    //  this = a.B (a = row vector, not in-place)
    void dot( InterFieldVec &a, InterFieldMat &B );
    void dot( Vec1dFieldType &a, InterFieldMat &B );
    void dot( InterFieldVec &a, Vec2dFieldType &B );

    // solve Ax=b for x
    void linsolve( InterFieldMat &A, InterFieldVec &b );

  private:

    // --- Member Variables --

    int _Ncol; // number of components
    std::vector< Field<FieldType>* > _DataPtr; // data

};


#endif // _INTERFIELDVEC_H
 
