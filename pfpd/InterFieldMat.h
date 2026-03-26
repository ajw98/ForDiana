/***************************************************************
*
* Class for the psuedospectral operators
* 
* DRT -- Fri, 13 Mar 2015
*
****************************************************************/

#ifndef _INTERFIELDMAT_H
#define _INTERFIELDMAT_H

#include "global.h"
#include "FieldStack.h"
#include "InterFieldVec.h"

class InterFieldVec;

class InterFieldMat : public FieldStack
{

  public:

    // --- Constructor/Destructor ---

    InterFieldMat( int nrow, int ncol );
    InterFieldMat( const InterFieldMat &rhs ); // copy constructor
    ~InterFieldMat();

    void CreateField();

    // --- operators (field algebra) ---

    // indexing operator
    InterFieldVec & operator[](const UInt i) { return _Data[i]; };
    const InterFieldVec & operator[](const UInt i) const { return _Data[i]; };

    // assignment (a = rhs)
    InterFieldMat & operator=(const InterFieldMat &rhs);  // other InterFieldMat
    InterFieldMat & operator=(const Vec2dFieldType &rhs); // const vector
    InterFieldMat & operator=(const FieldType &rhs);      // const scalar

    // in-place addition (a = a+rhs)
    InterFieldMat & operator += (const InterFieldMat &rhs);       // add field vector
    InterFieldMat & operator += (const Vec2dFieldType &rhs); // add const vector
    InterFieldMat & operator += (const FieldType &rhs);      // add const scalar

    // in-place subtraction (a = a-rhs)
    InterFieldMat & operator -= (const InterFieldMat &rhs);       // subtract field vector
    InterFieldMat & operator -= (const Vec2dFieldType &rhs); // subtract const vector
    InterFieldMat & operator -= (const FieldType &rhs);      // subtract const scalar

    // in-place elementwise/scalar multiplication (a = a*rhs)
    InterFieldMat & operator *= (const InterFieldMat &rhs);
    InterFieldMat & operator *= (const Vec2dFieldType &rhs);
    InterFieldMat & operator *= (const FieldType &rhs);

    // --- Member Functions ---

    // getters/setters
    inline UInt Nrow() const { return _Nrow; };
    inline UInt Ncol() const { return _Ncol; };

    void zero();
    void setflag_inrealspace( bool realspaceflag );
    bool getflag_inrealspace() const;
    FieldType max();

    // fourier transform and inverse
    void fft();
    void ifft();

    // matrix operations

    void eye( FieldType scalar );

    // matrix matrix multiply C = A.B
    void dot( const InterFieldMat &A, const InterFieldMat &B);
    void dot( InterFieldMat &A, InterFieldMat &B);

    void dot( const Vec2dFieldType &A, const InterFieldMat &B);
    void dot( Vec2dFieldType &A, InterFieldMat &B);

    void dot( const InterFieldMat &A, const Vec2dFieldType &B);
    void dot( InterFieldMat &A, Vec2dFieldType &B);

    // Matrix inversion (makes this = A^{-1})
    void invert( InterFieldMat &A ); 

  private:

    // --- Member Variables --

    int _Nrow; // number of rows
    int _Ncol; // number of cols
    std::vector< InterFieldVec > _Data;

};


#endif // _INTERFIELDMAT_H
 
