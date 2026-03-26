 /***************************************************************
*
* Base class for the time integration schemes
* 
* DRT -- Mon 16 May 2016
*
****************************************************************/

#ifndef _TIMEINT_BASE_H
#define _TIMEINT_BASE_H

#include "global.h"
#include "SmartField.h"
#include "SmartFieldVec.h"
#include "SmartFieldMat.h"
#include "Grid.h"
#include "Model_Visc_Base.h"
#include "Model_Energy_Base.h"
#include "Model_Reactions_Base.h"
#include "Model_Mobility_Base.h"
#include "Model_Factory.h"
#include "BCs_Base.h"
#include "BCs_Factory.h"
#include "Operators_Base.h"
#include "Operators_Factory.h"

class TimeInt_Base
{

  public:

    // --- constructor/destructor ---
    TimeInt_Base( int Ncomp, int NIcomp, Grid* CurrGrid );
    virtual ~TimeInt_Base();

    // --- Setup/Cleanup ---

    // setup_/cleanup_opearators are defaults for _PhiOp only.
    // Override these in ModelH
    virtual void setup_operators();
    virtual void cleanup_operators();

    void setup_models( int EnergyKey, int MobilityKey, int ViscKey, int ReactionsKey );
    void cleanup_models();

    // --- Input/Output ---
    virtual void set_params() = 0;
    virtual void read_initial_state() = 0;
    virtual void write_final_state() = 0;

    // --- Main integration loop ---
    virtual void outer_loop() = 0;

  protected:

    // Number of components
    int _Ncomp, _NIcomp;

    // BCs variables
    int _NBCs; // Number of boundary conditions
    int _BCFlag; // Which BCs are we using?

    // Grid variables
    Grid* _CurrGrid;
    UInt _Dim;

    // ICs
    int _BinaryInFlag, _BinaryOutFlag; // is IO binary or formatted?
    int _PSFlag; // Is the simulation pure PS?
    int _WriteVelFlag; // Are there velocities?
    int _ICFlag; // Which IC are we using?

    // Factories for models and operators
    BCs_Factory         *_BCFactory; // for making operators
    Operators_Factory   *_OpFactory; // for making operators
    Model_Factory       *_MFactory; // for making models

    // Model components
    Model_Energy_Base   *_EnergyModel;  // a free energy model pointer
    Model_Mobility_Base *_MobilityModel;  // a free energy model pointer
    Model_Visc_Base     *_ViscModel; // a viscosity model pointer
    Model_Reactions_Base *_ReactionsModel; // a reaction model pointer

    // Phi Operators (which every TimeInt must have)
    Operators_Base *_PhiOp; // Operators for Phi fields
    BCs_Base *_PhiBCs; // BCs object

    // Other
    bool check_phi( const SmartFieldVec &phi);
    bool check_phi( const SmartFieldVec &phi, const RealType &glassMax );
    void capField ( SmartField &F, const RealType k );
    void floorField ( SmartField &F, const RealType k );
    void floorField ( SmartFieldVec &F, const RealType k );

    // Fluctuation functions
    void calc_noise_phi ( SmartFieldMat const &M, RealType const rFac, RealType const dt, SmartFieldVec &noise );
    void calc_noise_vel ( SmartField const &eta, RealType const rFac, SmartFieldVec &noise);

};

#endif //_TIMEINT_BASE_H
 
