/***************************************************************
*
* Class for the psuedospectral operators
* 
* DRT -- Fri, 13 Mar 2015
*
****************************************************************/

#ifndef _SMARTFIELDOPMAT_H
#define _SMARTFIELDOPMAT_H

#include "global.h"
#include "SmartFieldMat.h"
#include "SmartFieldOp.h"

// Note: Need to initialize SmartFields *before*
//       making any SmartFieldOp objects!

typedef std::vector< SmartFieldOp > Vec1dSFO;
typedef std::vector< Vec1dSFO > Vec2dSFO;

class SmartFieldOpMat
{

  public:

    // --- Constructor/Destructor ---

    SmartFieldOpMat( UInt nrow, UInt ncol );
    SmartFieldOpMat( const SmartFieldOpMat &rhs ); // copy constructor
    ~SmartFieldOpMat();

    // --- Operators ---
    class ArrayProxy
    {
      public:

        ArrayProxy( Vec1dSFO* array_ptr ) : _ArrayPtr( array_ptr ) {};

        SmartFieldOp & operator[] (const UInt i) { return (*_ArrayPtr)[i]; };

      private:
        Vec1dSFO * _ArrayPtr;
    };

    ArrayProxy operator[](const UInt i) { return ArrayProxy( &(_Data[i]) ); };

    //SmartFieldOpVec & operator[](const UInt i) { return _Data[i]; };
    //const SmartFieldOpVec & operator[](const UInt i) const { return _Data[i]; };

    SmartFieldOpMat & operator=(const SmartFieldOpMat &rhs); // assignment operator
    SmartFieldOpMat & operator+=(const SmartFieldOpMat &rhs); // compound addition operator
    SmartFieldOpMat & operator*=(const FieldType rhs); // Scalar Multiply

    // --- Member Functions ---

    void AddBand( int bandnum, const SmartFieldMat & Band );
    void AddBand( int bandnum, const Vec2dFieldType & Band );
    void AddBand( int bandnum, const SmartFieldMat & Band, FieldType scalar );
    void GetBand( int bandnum, SmartFieldMat & Band ) const;
    void zero();

    void setflag_inrealspace( bool realspaceflag );
    bool getflag_inrealspace() const;

    // Perform the operation on a field: Operator . In = Out
    void OpDot( const SmartFieldVec & in, SmartFieldVec & out );

    // --- Getters/Setters ---

    inline UInt Nrow() const { return _Nrow; };
    inline UInt Ncol() const { return _Ncol; };
    inline UInt BW() const { return _Data[0][0].BW(); };
    inline UInt FDSize() const { return _Data[0][0].FDSize(); };
    inline UInt Nband() const { return _Data[0][0].Nband(); };
    inline std::vector< int > GetBandStencil() const { return _Data[0][0].GetBandStencil(); }

    // --- Reporters ---

    void ReportState();

  private:

    // --- Member Variables --

    UInt _Nrow;
    UInt _Ncol;
    Vec2dSFO _Data; // data

};


#endif // _SMARTFIELDOPMAT_H
 
