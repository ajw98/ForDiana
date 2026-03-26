/***************************************************************
*
* Class for the psuedospectral operators
* 
* DRT -- Fri, 13 Mar 2015
*
****************************************************************/

#ifndef _VECINTERFIELDVEC_H
#define _VECINTERFIELDVEC_H

#include "global.h"
#include "SmartFieldVec.h"
#include "InterFieldVec.h"
#include "InterFieldMat.h"

// Note: Need to initialize SmartFields *before*
//       making any BandField objects!

typedef std::vector< InterFieldVec > Vec1dIFV;
typedef std::vector< Vec1dIFV> Vec2dIFV;

class VecInterFieldVec
{

  public:

    // --- Constructor/Destructor ---

    VecInterFieldVec( int nfd, int ncomp );
    VecInterFieldVec( const VecInterFieldVec &rhs ); // copy constructor
    ~VecInterFieldVec();

    // --- operators (field algebra) ---

    VecInterFieldVec & operator=(const VecInterFieldVec &rhs); // assignment operator

    InterFieldVec & operator[](const UInt i) { return _Data[i]; };
    const InterFieldVec & operator[](const UInt i) const { return _Data[i]; };

    // --- Member Functions ---
    void zero();
    void setflag_inrealspace( bool realspaceflag );
    bool getflag_inrealspace() const;

    // --- Matrix solver functions ---
    void smart_to_inter( SmartFieldVec & v );
    void inter_to_smart( SmartFieldVec & v );

    // --- Getters/Setters ---

    inline UInt NFD() const { return _NFD; };
    inline UInt NIcomp() const { return _NIcomp; };

    // --- Reporters ---
    void ReportState( std::string filename );

  private:

    // --- Member Variables --

    int _NFD; // number of FD points
    int _NIcomp; // number of components
    Vec1dIFV _Data; // data

};


#endif // _VECINTERFIELDVEC_H
 
