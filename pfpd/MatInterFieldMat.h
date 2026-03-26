/***************************************************************
*
* Class for the psuedospectral operators
* 
* DRT -- Fri, 13 Mar 2015
*
****************************************************************/

#ifndef _MATINTERFIELDMAT_H
#define _MATINTERFIELDMAT_H

#include "global.h"
#include "SmartFieldOpMat.h"
#include "InterFieldVec.h"
#include "InterFieldMat.h"
#include "VecInterFieldVec.h"

typedef std::vector< InterFieldMat > Vec1dIFM;
typedef std::vector< Vec1dIFM > Vec2dIFM;

class MatInterFieldMat
{

  public:

    // --- Constructor/Destructor ---

    MatInterFieldMat( int nrow, int ncol, int nfmrow, int nfmcol );
    MatInterFieldMat( const MatInterFieldMat &rhs ); // copy constructor
    ~MatInterFieldMat();

    // --- operators (field algebra) ---

    class ArrayProxy
    {
      public:

        ArrayProxy( Vec1dIFM* array_ptr ) : _ArrayPtr( array_ptr ) {};

        InterFieldMat & operator[] (const UInt i) { return (*_ArrayPtr)[i]; };

      private:
        Vec1dIFM * _ArrayPtr;
    };

    ArrayProxy operator[](const UInt i) { return ArrayProxy( &(_Data[i]) ); };
    
    MatInterFieldMat & operator=(const MatInterFieldMat &rhs); // assignment operator

    // --- Member Functions ---
    void zero();
    void setflag_inrealspace( bool realspaceflag );
    bool getflag_inrealspace() const;

    // --- Matrix solver functions ---
    void smart_to_inter( SmartFieldOpMat & T );
    void inter_to_smart( SmartFieldOpMat & T );
    void band_LU();
    void band_solve( VecInterFieldVec & b, VecInterFieldVec & x );
    void band_solve( MatInterFieldMat & B, MatInterFieldMat & X );

    void matmul( MatInterFieldMat & A, MatInterFieldMat & B );
    void transpose();
    FieldType max();

    // --- Getters/Setters ---

    inline UInt Nrow() const { return _Nrow; };
    inline UInt Ncol() const { return _Ncol; };
    inline UInt NFMrow() const { return _NFMrow; };
    inline UInt NFMcol() const { return _NFMcol; };
    inline UInt BW() const { return _BW; };
    inline UInt NFD() const { return _NFD; };
    inline UInt NIcomp() const { return _NIcomp; };
    inline bool LU() const { return _LU; };

    // --- Reporters ---
    void ReportState();

  private:

    // --- Member Variables --

    int _Nrow; // in banded matrix is 2*BW+1
    int _Ncol; // in banded matrix is number of FD points
    int _NFMrow; // in banded matrix is number of components (NIcomp)
    int _NFMcol; // in banded matrix is number of components (NIcomp)
    int _BW; // in banded matrix is (_Nrow-1)/2
    int _NFD; // in banded matrix is _Ncol
    int _NIcomp; // in banded matrix is _NFMrow, _NFMcol
    bool _LU; // if _Data is LU decomposed then this is true
    Vec2dIFM _Data; // data

};


#endif // _MATINTERFIELDMAT_H
 
