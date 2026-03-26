/***************************************************************
*
* Base class for the time integration schemes
* 
* DRT -- Mon 16 May 2016
*
****************************************************************/

#ifndef _GRID_H
#define _GRID_H

#include "global.h"
#include "Field.h"
#include "FFTlayout.h"
#include "pllhandler.h"

class Grid
{

  public:

    // --- constructor/destructor ---

    Grid( pllhandler & pll, bool is_periodic );
    ~Grid();

    void write_PS_grid();
    void write_FD_grid();

    // getters/setters
    inline Cell* GetCell() { return _Cell; };
    inline FFTlayout* GetLayout() { return _CurrLayout; };

    inline UInt GridFlag() { return _GridFlag; }
    inline UInt FDSize() { return _FDSize; }
    inline UInt Dim() { return _Dim; }
    inline UInt PSDim() { return _PSDim; }
    inline UInt Nx() { return _Nx; }
    inline UInt Ny() { return _Ny; }
    inline UInt Nz() { return _Nz; }
    inline RealType Lx() { return _Lx; }
    inline RealType Ly() { return _Ly; }
    inline RealType Lz() { return _Lz; }
    inline RealType Dx() { return _Dx; }
    inline UInt NPW() { return _NPW; }
    inline RealType DxDim( UInt i ) { return _DxDim[i]; }

  protected:

    // pll, cell, layout, model objects
    // (need to keep these variables or they go out of scope!)
    pllhandler *_Pll; 
    Cell *_Cell;
    FFTlayout *_CurrLayout;
 
    // Grid variables
    bool _IsPeriodic;
    UInt _GridFlag, _Dim, _PSDim, _FDSize;
    UInt _Nx, _Ny, _Nz, _NPW;
    RealType _Lx, _Ly, _Lz, _Dx;
    Vec1dReal _DxDim; 

};

#endif //_GRID_H
 
