/***************************************************************
*
* Class for the time integration scheme:
*   - Model B (diffusive only)
*   - pseudospectral method (periodic BCs in all directions)
* 
* DRT -- Fri, 13 Mar 2015
*
****************************************************************/

#ifndef _TIMEINT_MODELB_H
#define _TIMEINT_MODELB_H

#include "global.h"
#include "FFTlayout.h"
#include "Field.h"
#include "SmartFieldVec.h"
#include "SmartFieldMat.h"
#include "SmartField.h"
#include "Model_Energy_Base.h"
#include "Model_Mobility_Base.h"
#include "Model_Reactions_Base.h"
#include "TimeInt_Base.h"
#include "histogram.h"

class TimeInt_ModelB : public TimeInt_Base
{

  public:

    // --- constructor/destructor ---
    TimeInt_ModelB( int Ncomp, int NIcomp, Grid* CurrGrid );
    ~TimeInt_ModelB();

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
    void step_phi( SmartFieldVec &phi, SmartFieldVec const &mu,
                   SmartFieldVec const &mu_lin, SmartFieldOpMat const &mu_lin_np1,
                   RealType dt, RealType t );
    void calc_lhs ( SmartFieldVec const &phi, Vec2dFieldType const &M_max,
                    SmartFieldOpMat const &mu_lin_np1, RealType const &dt, 
		    SmartFieldOpMat &lhs );
    RealType update_timestep( RealType phi_err, RealType dt );

    // --- Member Variables --

    // parameters
    int _TimeIntFlag; // flag for TimeInt driver
    RealType _TMax; // maximum time
    RealType _Dt0, _DtDisp; // initial timestep, display timestep
    RealType _DtMin, _DtMax; // timestep max/min
    int _OutputIntervalFlag; // 0 if linear, 1 if logarithmic 
    UInt _NtMax; // maximum number of possible timesteps
    UInt _NStepDisp; // number of steps between job.log output
    UInt _PhiIterMax; // max. # of times to iterate when phi doesn't converge
    UInt _VarDtFlag; // 1 if variable timestep, 0 if not
    UInt _CapPhi0Flag; // 1 if we artificially cap phi
    RealType _Phi0Max; // glassy maximum for phi[0]---make more general later!
    // parameters for adaptive tolerance-adaptive time stepping
    int _PhiErrTolFlag;
    RealType _PhiErrTol;
    RealType _PhiErrTolMax;
    RealType _PhiErrTolMin;
    RealType _PhiErrTolW;
    RealType _PhiErrTolMid;
    bool _EnableFluctuations;
    RealType _ReductionFactor;

    // Fields
    SmartFieldVec _Phi;

};

#endif // _TIMEINT_MODELB_H
 
