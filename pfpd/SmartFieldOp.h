/***************************************************************
*
* Class for the psuedospectral operators
* 
* DRT -- Fri, 13 Mar 2015
*
****************************************************************/

#ifndef _SMARTFIELDOP_H
#define _SMARTFIELDOP_H

#include "global.h"
#include "SmartField.h"

// Note: Need to initialize SmartFields *before*
//       making any SmartFieldOp objects!

class SmartFieldOp
{

  public:

    // --- Constructor/Destructor ---

    SmartFieldOp();
    SmartFieldOp( const SmartFieldOp &rhs ); // copy constructor
    ~SmartFieldOp();

    // --- Operators ---
    SmartFieldOp & operator=(const SmartFieldOp &rhs); // assignment operator
    SmartFieldOp & operator+=(const SmartFieldOp &rhs); // compound addition operator
    SmartFieldOp & operator*=(const FieldType rhs); // Scalar Multiply

    // --- Member Functions ---

    void SetBand( int bandnum, SmartField & Band );
    void AddBand( int bandnum, const SmartField & Band );
    void AddBand( int bandnum, FieldType Band );
    void GetBand( int bandnum, SmartField & Band ) const;

    // Perform the operation on a field: Operator . In = Out
    void OpDot( const SmartField & in, SmartField & out );

    void zero();
    void setflag_inrealspace( bool realspaceflag );
    bool getflag_inrealspace() const;

    // --- Getters/Setters ---

    inline UInt BW() const { return _BW; };
    inline UInt Nband() const { return _BandStencil.size(); };
    inline std::vector< int > GetBandStencil() const { return _BandStencil; }
    int GetStencilPos( int bandnum ) const;
    inline UInt FDSize() const { return _Data[0].GetSize(); };

    // --- Reporters ---

    void ReportState();
    void writefield(std::string filename, bool omitimaginary);
    inline void writerealfield(std::string filename) { writefield(filename, true); };
    inline void writecplxfield(std::string filename) { writefield(filename, false); };

  private:

    // --- Member Variables --

    int _BW; // band width

    std::vector< int > _BandStencil; // band stencil
    std::vector< SmartField > _Data; // data

};


#endif // _SMARTFIELDOP_H
 
