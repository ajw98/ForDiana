/***************************************************************
*
* Class for the time integration scheme
* 
* DRT -- Fri, 13 Mar 2015
*
****************************************************************/

#ifndef _TIMEINT_MODELH_H
#define _TIMEINT_MODELH_H

#include "global.h"
#include "FFTlayout.h"
#include "pllhandler.h"
#include "Field.h"
#include "SmartFieldVec.h"
#include "SmartFieldMat.h"
#include "SmartField.h"
#include "TimeInt_Base.h"
#include "BCs_Base.h"
#include "BCs_Factory.h"
#include "histogram.h"

class TimeInt_ModelH : public TimeInt_Base
{

  public:

    // --- constructor/destructor ---
    TimeInt_ModelH( int Ncomp, int NIcomp, Grid* CurrGrid );
    ~TimeInt_ModelH();

    // --- Setup/Cleanup ---
    void setup_operators();
    void cleanup_operators();

    // --- Input/Output ---
    void set_params();
    void read_initial_state();
    void write_final_state();

    // Main integration loop
    void outer_loop();
    void outer_loop_const();
    void outer_loop_adap();

  private:

    // --- Member Functions ---

    // I/O
    void write_curr_state( int n_disp );

    // Integration routines
    void step_phi( SmartFieldVec &phi, SmartFieldVec const &vel, SmartFieldVec const &mu,
                   SmartFieldVec const &mu_lin, SmartFieldOpMat const &mu_lin_np1,
                   RealType dt, RealType t );
    void step_vel( SmartFieldVec const &phi, SmartFieldVec &vel, SmartFieldVec const &mu,
                   UInt &iter, RealType t );
    RealType update_timestep( SmartFieldVec &vel, 
                              RealType phi_err, 
                              RealType dt );
    void f_sigmoid(const SmartFieldVec &phi, SmartFieldVec &BodyForce, int comp, RealType phiT, RealType w, 
                    Vec1dReal scale,SmartField &finvert);

    // --- Member Variables --

    // TimeInt parameters
    int _TimeIntFlag; // flag for TimeInt driver
    RealType _TMax; // maximum time
    RealType _Dt0, _DtDisp; // initial timestep, display timestep
    RealType _DtMin, _DtMax; // timestep max/min
    RealType _PhiErrTol, _VelErrTol; // error tolerances for phi and velocity
    int _VarViscFlag; // 0 if constant viscosity, 1 if variable viscosity
    int _OutputIntervalFlag; // 0 if linear, 1 if logarithmic 
    UInt _NtMax; // maximum number of possible timesteps
    UInt _NStepDisp; // number of steps between job.log output
    UInt _PhiIterMax; // max # of times to iterate when phi doesn't converge
    UInt _VarDtFlag; // 1 if variable timestep, 0 if not
    UInt _CapPhiFlag; // 1 if we artificially cap phi
    RealType _GlassMax; // maximum phi[0] concentration for glassy asymptote
    
    int _BodyForceFlag;
    bool _EnableFluctuations;
    RealType _ReductionFactor;
    bool _AdvFlag;

    RealType _W;
    RealType _PhiT;
    Vec1dReal _Scale;
    int bf_comp;
    SmartField finvert;

    // Fields
    SmartFieldVec _Phi;
    SmartFieldVec _Vel;
    SmartFieldVec _BodyForce;

    // tmp vecs with persistent storage for performance - do not rely on data persistence
    // between levels of scope!
    SmartFieldVec _tmpvec, _tmpvec2;
    SmartFieldMat _Mobility;

    // Operators (in addition to _Phi operators from base class)
    Operators_Base* _VelOp; // velocity operators
    BCs_Base* _VelBCs;

};

#endif // _TIMEINT_MODELH_H
 
