/***************************************************************
*
* Class for the time integration scheme
* 
* DRT -- Fri, 13 Mar 2015
*
****************************************************************/
//
//
// d phi_i / d t + v.del phi_i = sum_j del . ( M_ij del mu_j )
// 0 = - del p + del . (eta (del v + del v^T)) - Nr sum_i del mu_i
// 
//

#include "TimeInt_ModelH.h"

// ----------------- Constructor/Destructor -----------------

TimeInt_ModelH :: TimeInt_ModelH( int Ncomp, int NIcomp,
                                  Grid* CurrGrid ):
TimeInt_Base( Ncomp, NIcomp, CurrGrid ),
_Phi( NIcomp ),
_Vel( CurrGrid->Dim() ),
_BodyForce( CurrGrid->Dim() ),
_VelOp( NULL ),
_VelBCs( NULL ),
_tmpvec( NIcomp ),
_tmpvec2( NIcomp ),
_Mobility (NIcomp, NIcomp)
// {{{
{

  std::cout << std::endl;
  std::cout << " * Time Integration derived class initialization\n";
  std::cout << "   - Model H (convection/diffusion)" << std::endl;

//  // -- Overall Model Params ---
//  // *** TODO Get Nr from Energy Model Params
//  _Pe = 1.; // Pe = 1 for periodic BCs
//  _Ca = 1./Nr; // Ca = 1/Nr for periodic BCs
//
//  std::cout << "   - General Model Parameters:\n" << std::scientific;
//  std::cout << "      Pe = " << _Pe << "\n";
//  std::cout << "      Ca = " << _Ca << "\n";

} // }}}

TimeInt_ModelH :: ~TimeInt_ModelH()
{ // {{{
} // }}}

// ----------------- Setup/Cleanup -----------------

void TimeInt_ModelH :: setup_operators()
{ // {{{

  // set up the BCs and operators for phi
  {
    std::string filename( "params_BCs_phi.in" );
    
    if (jsonFile(filename.c_str()))
    {
      // Read in BC key from file
      std::ifstream f2(filename.c_str());
      if ( f2.is_open() )
      {
        rapidjson::IStreamWrapper isw(f2);
        rapidjson::Document d;
        d.ParseStream(isw);
        _BCFlag = d["BCFlag"].GetInt();
      }
      else
      {
        std::cout << "*** Error in TimeInt_ModelH :: setup_operators() ***\n";
        std::cout << "*** Cannot open " << filename << " ***\n";
        exit(1);
      }
    }
    else
    {//Legacy file input
      std::string tmp_str;
      std::ifstream f2(filename.c_str());
      if ( f2.is_open() )
      {
        std::getline(f2, tmp_str, '#');
        _BCFlag = atoi(tmp_str.c_str());
        std::getline(f2, tmp_str, '\n');
      }
      else
      {
        std::cout << "*** Error in TimeInt_ModelH :: setup_operators() ***\n";
        std::cout << "*** Cannot open " << filename << " ***\n";
        exit(1);
      }
    }

    // call factory to create new BC object based on the kind of BCs
    _PhiBCs=_BCFactory->make_BCs( _BCFlag, 0, filename, _NIcomp, _CurrGrid );

    // call factory to make new operators object for phi
    _PhiOp = _OpFactory->make_operators( _CurrGrid, _PhiBCs );
  }

  // set up the BCs and operators for vel
  {
    std::string filename( "params_BCs_vel.in" );
    
    if (jsonFile(filename.c_str()))
    {
      // Read in BC key from file
      std::ifstream f2(filename.c_str());
      if ( f2.is_open() )
      {
        rapidjson::IStreamWrapper isw(f2);
        rapidjson::Document d;
        d.ParseStream(isw);
        _BCFlag = d["BCFlag"].GetInt();
      }
      else
      {
        std::cout << "*** Error in TimeInt_ModelH :: setup_operators() ***\n";
        std::cout << "*** Cannot open " << filename << " ***\n";
        exit(1);
      }
    }
    else
    {//Legacy file format
      std::string tmp_str;
      std::ifstream f2(filename.c_str());
      if ( f2.is_open() )
      {
        std::getline(f2, tmp_str, '#');
        _BCFlag = atoi(tmp_str.c_str());
        std::getline(f2, tmp_str, '\n');
      }
      else
      {
        std::cout << "*** Error in TimeInt_ModelH :: setup_operators() ***\n";
        std::cout << "*** Cannot open " << filename << " ***\n";
        exit(1);
      }
    }
    // call factory to create new BC object based on the kind of BCs
    _VelBCs=_BCFactory->make_BCs( _BCFlag, 1, filename, _NIcomp, _CurrGrid );

    // call factory to make new operators object for phi
    _VelOp = _OpFactory->make_operators( _CurrGrid, _VelBCs );
  }

} // }}}

void TimeInt_ModelH :: cleanup_operators()
{ // {{{

  delete _PhiBCs;
  _OpFactory->recycle_operators( _PhiOp );

  delete _VelBCs;
  _OpFactory->recycle_operators( _VelOp );

} // }}}

// ----------------- Input/Output -----------------

void TimeInt_ModelH :: set_params()
{//{{{
  
  //Define default parameter values first
  _EnableFluctuations = false; 
  _ReductionFactor = 1.;
  _AdvFlag = true;

  if (jsonFile("params_TimeInt.in"))
  {
    std::ifstream f1("params_TimeInt.in");
    if ( f1.is_open() )
    {
      rapidjson::IStreamWrapper isw(f1);
      rapidjson::Document d;
      d.ParseStream(isw);

      _TimeIntFlag = d["TimeIntFlag"].GetInt();
      
      _TMax = d["t_max"].GetDouble();

      _Dt0 = d["dt0"].GetDouble();

      _DtMax = d["max_dt"].GetDouble();

      _DtMin = d["min_dt"].GetDouble();

      _DtDisp = d["dt_disp"].GetDouble();

      _PhiErrTol = d["phi_err_tol_max"].GetDouble();

      _VelErrTol = d["v_err_tol"].GetDouble();

      _VarViscFlag = d["var_visc_flag"].GetInt();

      _OutputIntervalFlag = d["output_interval_flag"].GetInt();

      _NtMax = d["Nt_max"].GetInt();

      _NStepDisp = d["N_step_disp"].GetInt();

      _PhiIterMax = d["phi_iter_max"].GetInt();

      _VarDtFlag = d["var_dt_flag"].GetInt();

      _CapPhiFlag = d["CapPhi0Flag"].GetInt();

      _GlassMax = d["phi0_max"].GetDouble();

      _BodyForceFlag = d["BodyForce_flag"].GetInt();
      
      int fluct = d["enablefluctuations"].GetInt();
      _EnableFluctuations = (fluct == 1) ? true:false;	
      
      _ReductionFactor = d["reductionFactor"].GetDouble();
      
      int adv_in = d["advection_flag"].GetInt();
      _AdvFlag = (adv_in == 1) ? true:false;
      
      _PhiT = d["bf_phiT"].GetDouble();

      _W = d["bf_w"].GetDouble();
      bf_comp = d["bf_comp"].GetInt();
      
      //_Scale = d["bf_scale"].GetDouble();  
    
      _Scale.resize(_Dim, 0.);
      rapidjson::Value& Scale_in = d["bf_scale"];
      for( int i=0; i<_Dim; i++ )
      {
        _Scale[i] = Scale_in[i].GetDouble();
      }
    }
    else
    {
      std::cout << "cannot open params_TimeInt.in" << std::endl;
      exit(1);
    }
    f1.close();
  }
  else
  {//Legacy file format
    std::string tmp_str;

    std::ifstream f1("params_TimeInt.in");
    if ( f1.is_open() )
    {
      std::getline(f1, tmp_str, '#');
      _TimeIntFlag = atof(tmp_str.c_str());
      std::getline(f1, tmp_str, '\n');

      std::getline(f1, tmp_str, '#');
      _TMax = atof(tmp_str.c_str());
      std::getline(f1, tmp_str, '\n');

      std::getline(f1, tmp_str, '#');
      _Dt0 = atof(tmp_str.c_str());
      std::getline(f1, tmp_str, '\n');

      std::getline(f1, tmp_str, '#');
      _DtMax = atof(tmp_str.c_str());
      std::getline(f1, tmp_str, '\n');

      std::getline(f1, tmp_str, '#');
      _DtMin = atof(tmp_str.c_str());
      std::getline(f1, tmp_str, '\n');

      std::getline(f1, tmp_str, '#');
      _DtDisp = atof(tmp_str.c_str());
      std::getline(f1, tmp_str, '\n');
      
      std::getline(f1, tmp_str, '#');
      _PhiErrTol = atof(tmp_str.c_str());
      std::getline(f1, tmp_str, '\n');

      std::getline(f1, tmp_str, '#');
      _VelErrTol = atof(tmp_str.c_str());
      std::getline(f1, tmp_str, '\n');

      std::getline(f1, tmp_str, '#');
      _VarViscFlag = atoi(tmp_str.c_str());
      std::getline(f1, tmp_str, '\n');

      std::getline(f1, tmp_str, '#');
      _OutputIntervalFlag = atoi(tmp_str.c_str());
      std::getline(f1, tmp_str, '\n');

      std::getline(f1, tmp_str, '#');
      _NtMax = atoi(tmp_str.c_str());
      std::getline(f1, tmp_str, '\n');

      std::getline(f1, tmp_str, '#');
      _NStepDisp = atoi(tmp_str.c_str());
      std::getline(f1, tmp_str, '\n');

      std::getline(f1, tmp_str, '#');
      _PhiIterMax = atoi(tmp_str.c_str());
      std::getline(f1, tmp_str, '\n');

      std::getline(f1, tmp_str, '#');
      _VarDtFlag = atoi(tmp_str.c_str());
      std::getline(f1, tmp_str, '\n');

      std::getline(f1, tmp_str, '#');
      _CapPhiFlag = atoi(tmp_str.c_str());
      std::getline(f1, tmp_str, '\n');

      std::getline(f1, tmp_str, '#');
      _GlassMax = atof(tmp_str.c_str());
      std::getline(f1, tmp_str, '\n');

      std::getline(f1, tmp_str, '#');
      _BodyForceFlag = atoi(tmp_str.c_str());
      std::getline(f1, tmp_str, '\n');
      
      // Read flag for fluctuating simulations.
      // Tolerate old input files not containing the flag.
      std::getline(f1, tmp_str, '#');
      if(!f1.eof())
      {
        int fluct = atoi(tmp_str.c_str());
        std::getline(f1, tmp_str, '\n');
        if(fluct==1)
          _EnableFluctuations = true;
        else
          _EnableFluctuations = false;	
      }
      
      std::getline(f1, tmp_str, '#');
      if(!f1.eof())
      {
        _ReductionFactor = atof(tmp_str.c_str());
        std::getline(f1, tmp_str, '\n');
      }
    }
    else
    {
      std::cout << "cannot open params_TimeInt.in" << std::endl;
      exit(1);
    }
    f1.close();
  }

  // Get BodyForce from file 
  if ( _BodyForceFlag ==1)
  {
    char filename[40];
    char dimension[3];
    dimension[0] = 'x';
    dimension[1] = 'y';
    dimension[2] = 'z';

    _BodyForce.setflag_inrealspace(true);
    for (int i=0; i<_Dim; i++)
    {
      sprintf( filename, "bf%c.in", dimension[i] );
      _BodyForce[i].readfield(filename);
    }
  }
  else if (_BodyForceFlag == 2)
  {
    finvert.readfield("invert.in");
  }
  std::cout << "   - Parameters:\n" << std::scientific;
  std::cout << "      TimeIntFlag        = " << _TimeIntFlag << std::endl;
  std::cout << "      TMax               = " << _TMax << std::endl;
  std::cout << "      Dt0                = " << _Dt0 << std::endl;
  std::cout << "      DtMax              = " << _DtMax << std::endl;
  std::cout << "      DtMin              = " << _DtMin << std::endl;
  std::cout << "      DtDisp             = " << _DtDisp << std::endl;
  std::cout << "      PhiErrTol          = " << _PhiErrTol << std::endl;
  std::cout << "      VelErrTol          = " << _VelErrTol << std::endl;
  std::cout << "      VarViscFlag        = " << _VarViscFlag << std::endl;
  std::cout << "      OutputIntervalFlag = " << _OutputIntervalFlag << std::endl;
  std::cout << "      NtMax              = " << _NtMax << std::endl;
  std::cout << "      NStepDisp          = " << _NStepDisp << std::endl;
  std::cout << "      PhiIterMax         = " << _PhiIterMax << std::endl;
  std::cout << "      VarDtFlag          = " << _VarDtFlag << std::endl;
  std::cout << "      BodyForceFlag      = " << _BodyForceFlag << std::endl;
  std::cout << "      EnableFluctuations = " << _EnableFluctuations << std::endl;
  std::cout << "      ReductionFactor    = " << _ReductionFactor << std::endl;
  std::cout << "      AdvectionFlag      = " << _AdvFlag << std::endl;

  // Warning print out if fluctuations enabled with standard advection term formulation
  // Tests have shown this combination breaks conservation of mass
  if (_EnableFluctuations and _AdvFlag)
  {
    std::cout << std::endl;
    std::cout << "  *** Warning: fluctuations enabled with std. advection formulation!" << std::endl;
    std::cout << "  *** Check for composition drift in compositions.dat " << std::endl;
  }

  return;

} // }}}

void TimeInt_ModelH::read_initial_state()
{ // {{{

  // read in initial conditions
  char filename[40];
  for (int i=0; i<_NIcomp; i++)
  {
    sprintf( filename, "phi%d.in", i+1 );
    _Phi[i].readfield(filename);
  }

  char dimension[3];
  dimension[0] = 'x';
  dimension[1] = 'y';
  dimension[2] = 'z';

  for (int i=0; i<_Dim; i++)
  {
    sprintf( filename, "v%c.in", dimension[i] );
    _Vel[i].readfield(filename);
  }

} // }}}

void TimeInt_ModelH::write_final_state()
{ // {{{

  // write out final state
  char filename[40];
  for (int i=0; i<_NIcomp; i++)
  {
    sprintf( filename, "phi%d.out", i+1 );
    _Phi[i].writerealfield(filename);
  }

  char dimension[3];
  dimension[0] = 'x';
  dimension[1] = 'y';
  dimension[2] = 'z';

  for (int i=0; i<_Dim; i++)
  {
    sprintf( filename, "v%c.out", dimension[i] );
     _Vel[i].writerealfield(filename);
  }

} // }}}

void TimeInt_ModelH :: write_curr_state( int n_disp )
{ // {{{

  // write out current state
  char filename[40];
  for (int i=0; i<_NIcomp; i++)
  {
    if (n_disp>=0)
      sprintf( filename, "phi%d_%05d.dat", i+1, n_disp );
    else
      sprintf( filename, "phi%d_ERR.dat", i+1 );
    _Phi[i].writerealfield(filename);
  }

  char dimension[3];
  dimension[0] = 'x';
  dimension[1] = 'y';
  dimension[2] = 'z';

  for (int i=0; i<_Dim; i++)
  {
    if (n_disp>=0)
      sprintf( filename, "v%c_%05d.dat", dimension[i], n_disp );
    else
      sprintf( filename, "v%c_ERR.dat", dimension[i] );

     _Vel[i].writerealfield(filename);
  }

} // }}}

// ----------------- Integration Routines -----------------
void TimeInt_ModelH :: outer_loop()
{ // {{{
  // wrapper around the 2 cases:
  // (i) adaptive time-stepping by DRT
  // (ii) constant time-stepping 

  if (_VarDtFlag)
  {
    outer_loop_adap();
  } 
  else
  {
    outer_loop_const();
  } 

} // }}}

void TimeInt_ModelH :: outer_loop_const()
{ // {{{
  // constant time-stepping
  SmartFieldVec vel_prev( _Vel );
  SmartFieldVec vel_slope( _Vel );

  RealType t = 0.; // total time
  RealType dt = _Dt0; // time step
  RealType dt_new = _Dt0; // time step
  RealType t_check = 0; // time when write to disk
  bool is_physical; // make sure phi is between zero and one
  
  UInt Nt = 0; // number of timesteps taken
  UInt n_disp = 0; // number of fields written to file
  UInt phi_iter = 0; // number of iterations for phi
  UInt vel_iter = 0; // number of iterations from velocity step

  // record which times are written
  std::ofstream timefile;
  timefile.open("time.dat");
  timefile << "# Columns: Ndisp time Nt\n" << std::scientific;

  // diagnostics of step-doubling
  std::ofstream stepfile; 
  stepfile.open("stepsize.dat");
  stepfile << "# Columns: N t dt phi_iter vel_iter\n";
  stepfile << std::scientific;

  std::ofstream compfile;
  compfile.open("compositions.dat");
  compfile << std::scientific;
  compfile << "# Columns: N phi[0].avg ... phi[N].avg phi[0].min ... phi[N].min phi[0].max ... phi[N].max\n";

  // print to screen
  std::cout << std::endl;
  std::cout << " * Beginning Time Integration Loop\n";
  std::cout << " * Constant time-stepping; no error estimates for phi\n";
  printf("%8s", "Nt");
  printf("%12s", "t");
  printf("%12s", "dt");
  printf("%12s", "phi_iter");
  printf("%12s", "vel_iter");
  std::cout << std::endl;

  // write initial state
  stepfile << Nt << "\t" << t << "\t" << dt << "\t";
  stepfile << phi_iter << "\t" << vel_iter << std::endl;
  printf("%8lu%12.3e%12.3e%12lu%12lu", Nt, t, dt, 
                                             phi_iter, vel_iter);
  std::cout << std::endl;

  printf("      * writing ndisp = %05lu to disk\n", n_disp);
  write_curr_state( n_disp );
  timefile << n_disp << "\t" << t << "\t" << Nt << std::endl;
  n_disp++;

  if ( _OutputIntervalFlag == 0 )
  {
    t_check = n_disp*_DtDisp;
  }
  else
  {
    t_check = pow(10., n_disp*_DtDisp); // minimum output is at t = 1
  }

  vel_slope = FieldType(0.);

  // Thermally averaged structure factor for all pairs
  SmartFieldVec SK((_NIcomp*(_NIcomp+1))/2); // Species-species structure factors in a packed storage format
  UInt numSKsamples(0);
  histogram<FieldType, double> SKhist(*_CurrGrid->GetLayout(), true);
  SK.setflag_inrealspace(false);
  SK.zero();

  while( t < _TMax && Nt < _NtMax )
  {
    // Keep a copy of mu_n and mu_n+1/2 to avoid repeat evaluations in stepper routines
    SmartFieldVec mu_current(_NIcomp);
    SmartFieldVec mu_current_linexp(_NIcomp);
    SmartFieldOpMat mu_current_linimp(_NIcomp,_NIcomp);
    mu_current.setflag_inrealspace(false);

    // Compute and store current mu before any adaptive time stepper looping
    _Phi.fft();
    _EnergyModel->calc_mu(_Phi, mu_current);
    _EnergyModel->calc_mu_lin_exp(_Phi, mu_current_linexp);
    _EnergyModel->calc_mu_lin_imp(_Phi, mu_current_linimp);
    _Phi.ifft();

    step_phi( _Phi, _Vel, mu_current, mu_current_linexp, mu_current_linimp, dt, t ); 
    is_physical = check_phi(_Phi,_GlassMax);
    if ( !is_physical )
    {
      std::cout << "*** (ERROR) TimeInt_ModelH::outer_loop_const ***\n";
      std::cout << "Phi becomes unphysical at this timestep (dt = ";
      std::cout << dt << ")." << std::endl;
      exit(-1);
    } 
    
    if (_VarViscFlag)
    { 
      // Save a copy of old _Vel before modification
      vel_prev = _Vel;

      // Use 1st order continuation for velocity guess
      _Vel.xpby_inplace(vel_slope, dt);
        
      // Chemical potential needs to be updated for the velocity calculation
      _Phi.fft();
      _EnergyModel->calc_mu(_Phi, mu_current);
      _Phi.ifft();

      // Update _Vel
      step_vel( _Phi, _Vel, mu_current, vel_iter, t );
       
      // Update vel_slope for future continuation in next time-step
      vel_slope = _Vel;
      vel_slope.axpby_inplace(vel_prev, 1./dt, -1./dt);
    }
    else  
    {
      //constant-viscosity case does not need the 1st order continuation;
      //no iterations required within step_vel
      
      // Chemical potential needs to be updated for the velocity calculation
      _Phi.fft();
      _EnergyModel->calc_mu(_Phi, mu_current);
      _Phi.ifft();

      // Update _Vel
      step_vel( _Phi, _Vel, mu_current, vel_iter, t );
    }

    t = t + dt;
    Nt++;

    bool enableSK(false);
    if(enableSK)
    {
      numSKsamples++;
      double V(_CurrGrid->GetCell()->getVol());
      SmartField tmp;
      for(int i=0; i<_NIcomp; i++)
        _Phi[i].fft_rtok();
      // S(k) = <rho_k rho_{-k}> = <|rho_k|^2>
      int k=0;
      for(int i=0; i<_NIcomp; i++)
        for(int j=i; j<_NIcomp; j++, k++)
        {
          tmp = _Phi[i];
          tmp.prodconjg(_Phi[j]);
          SK[k].xpby_inplace(tmp,V);
        }
      // Revert phi fields to (r) space
      for(int i=0; i<_NIcomp; i++)
        _Phi[i].fft_ktor();
    }

    // --- write diagnostics and data ---
    // step doubling diagnostics
    if ( Nt % _NStepDisp == 0 || t >= t_check || t == _TMax || Nt == _NtMax )
    {
      stepfile << Nt << "\t" << t << "\t" << dt << "\t";
      stepfile << phi_iter << "\t" << vel_iter << std::endl;
      printf("%8lu%12.3e%12.3e%12lu%12lu", Nt, t, dt, phi_iter, vel_iter);
      std::cout << std::endl;

      // Average and extrema compositions---useful for debugging
      // Note: organzied this way so that the average compositions are the left-most cols. 
      compfile << Nt << "\t";
      for(int i=0; i<_NIcomp; i++)
        compfile << _Phi[i].integrate().real() << "\t";
      for(int i=0; i<_NIcomp; i++)
        compfile << _Phi[i].minsigned().real() << "\t";
      for(int i=0; i<_NIcomp; i++)
        compfile << _Phi[i].maxsigned().real() << "\t";
      compfile << std::endl;
    }

    // write data to file
    if ( t >= t_check || t == _TMax || Nt == _NtMax )
    {
      printf("      * writing ndisp = %05lu to disk\n", n_disp);
      write_curr_state( n_disp );
      timefile << n_disp << "\t" << t << "\t" << Nt << std::endl;
      n_disp++;

      if ( _OutputIntervalFlag == 0 )
      {
        t_check = n_disp*_DtDisp;
      }
      else
      {
        t_check = pow(10., n_disp*_DtDisp); // minimum output is at t = 1
      }
    }

    // change dt_new so that output is on the exact time we want
    if ( (t+dt_new) > t_check ) 
    {
      dt_new = t_check-t;
      // if dt < _DtMin, it will crash
      dt_new = std::max(dt_new, _DtMin); 
    }

  }

  stepfile.close();
  timefile.close();

  std::cout << std::endl;
  std::cout << " * Finished time integration\n";
  std::cout << "   - Final Time: " << t << std::endl;
  std::cout << "   - No. of Steps: " << Nt << std::endl;

} // }}}

void TimeInt_ModelH :: outer_loop_adap()
{ // {{{

  SmartFieldVec phi_half( _Phi );
  SmartFieldVec phi_whole( _Phi );
  SmartFieldVec vel_half( _Vel );
  SmartFieldVec vel_slope( _Vel );

  RealType t = 0.; // total time
  RealType dt = _Dt0; // time step
  RealType dt_new = _Dt0; // time step
  RealType t_check = 0; // time when write to disk
  RealType phi_err = 1.; // max truncation error
  bool is_physical; // make sure phi is between zero and one
  
  UInt Nt = 0; // number of timesteps taken
  UInt n_disp = 0; // number of fields written to file
  UInt phi_iter = 0; // number of iterations for phi
  UInt vel_iter = 0; // number of iterations from velocity step

  // record which times are written
  std::ofstream timefile;
  timefile.open("time.dat");
  timefile << "# Columns: Ndisp time Nt\n" << std::scientific;

  // diagnostics of step-doubling
  std::ofstream stepfile; 
  stepfile.open("stepsize.dat");
  stepfile << "# Columns: N t dt phi_err phi_iter vel_iter\n";
  stepfile << std::scientific;

  std::ofstream compfile;
  compfile.open("compositions.dat");
  compfile << std::scientific;
  compfile << "# Columns: N phi[0].avg ... phi[N].avg phi[0].min ... phi[N].min phi[0].max ... phi[N].max\n";

  // print to screen
  std::cout << std::endl;
  std::cout << " * Beginning Time Integration Loop\n";
  printf("%8s", "Nt");
  printf("%12s", "t");
  printf("%12s", "dt");
  printf("%12s", "phi_err");
  printf("%12s", "phi_iter");
  printf("%12s", "vel_iter");
  std::cout << std::endl;

  // write initial state
  stepfile << Nt << "\t" << t << "\t" << dt << "\t";
  stepfile << phi_err << "\t" << phi_iter << "\t" << vel_iter << std::endl;
  printf("%8lu%12.3e%12.3e%12.3e%12lu%12lu", Nt, t, dt, phi_err, 
                                             phi_iter, vel_iter);
  std::cout << std::endl;

  printf("      * writing ndisp = %05lu to disk\n", n_disp);
  write_curr_state( n_disp );
  timefile << n_disp << "\t" << t << "\t" << Nt << std::endl;
  n_disp++;

  if ( _OutputIntervalFlag == 0 )
  {
    t_check = n_disp*_DtDisp;
  }
  else
  {
    t_check = pow(10., n_disp*_DtDisp); // minimum output is at t = 1
  }

  vel_slope = FieldType(0.);

  // Thermally averaged structure factor for all pairs
  SmartFieldVec SK((_NIcomp*(_NIcomp+1))/2); // Species-species structure factors in a packed storage format
  UInt numSKsamples(0);
  histogram<FieldType, double> SKhist(*_CurrGrid->GetLayout(), true);
  SK.setflag_inrealspace(false);
  SK.zero();

  while( t < _TMax && Nt < _NtMax )
  {
    // Keep a copy of mu_n and mu_n+1/2 to avoid repeat evaluations in stepper routines
    SmartFieldVec mu_current(_NIcomp);
    SmartFieldVec mu_current_linexp(_NIcomp);
    SmartFieldOpMat mu_current_linimp(_NIcomp,_NIcomp);
    SmartFieldVec mu_half(_NIcomp);
    SmartFieldVec mu_half_linexp(_NIcomp);
    SmartFieldOpMat mu_half_linimp(_NIcomp,_NIcomp);
    mu_current.setflag_inrealspace(false);
    mu_half.setflag_inrealspace(false);

    // Compute and store current mu before any adaptive time stepper looping
    _Phi.fft();
    _EnergyModel->calc_mu(_Phi, mu_current);
    _EnergyModel->calc_mu_lin_exp(_Phi, mu_current_linexp);
    _EnergyModel->calc_mu_lin_imp(_Phi, mu_current_linimp);
    _Phi.ifft();

    // --- step doubling for phi and vel ---
    phi_err = 1.;
    phi_iter = 0;
    while( phi_err >= _PhiErrTol && phi_iter <= _PhiIterMax )
    {

      // Step doubling steps 
      // (Make sure that phi stays bounded between 0 and 1)
      is_physical = false;
      while (!is_physical)
      {

        if (_VarDtFlag)
        {
          dt = dt_new;
        }

        if ( dt < _DtMin )
        {
          std::cout << "*** (ERROR) Insufficient progress ***\n";
          std::cout << "Timestep below _DtMin. Aborting Calculation.\n";
          std::cout << std::endl;
          write_curr_state(-1);
          exit(1);
        }

        // first half step for phi
        phi_half = _Phi;
        step_phi( phi_half, _Vel, mu_current, mu_current_linexp, mu_current_linimp, dt/2., t );
        is_physical = check_phi(phi_half,_GlassMax);
        if ( !is_physical ) 
        { 
          if (_VarDtFlag)
          {
            dt_new = 0.8*dt;
            continue;
          }
          else
          {
            std::cout << "*** (ERROR) TimeInt_ModelH::outer_loop ***\n";
            std::cout << "Phi becomes unphysical at this timestep (dt = ";
            std::cout << dt << ")." << std::endl;
            exit(-1);
          }
        }

        // Compute and store mu at half step
        phi_half.fft();
        _EnergyModel->calc_mu(phi_half, mu_half);
        _EnergyModel->calc_mu_lin_exp(phi_half, mu_half_linexp);
        _EnergyModel->calc_mu_lin_imp(phi_half, mu_half_linimp);
        phi_half.ifft();

        // half step for vel - estimate using continuation
        vel_half = _Vel;
        vel_half.xpby_inplace(vel_slope, 0.5*dt);
        // Solve v
        step_vel( phi_half, vel_half, mu_half, vel_iter, t );
        // Update slope = 2*(vel_half-_Vel)/dt
        vel_slope = vel_half;
        vel_slope.axpby_inplace(_Vel, 2./dt, -2./dt);

        // second half step for phi
        step_phi( phi_half, vel_half, mu_half, mu_half_linexp, mu_half_linimp, dt/2., t );
        is_physical = check_phi(phi_half,_GlassMax);
        if ( !is_physical ) 
        { 
          if (_VarDtFlag)
          {
            dt_new = 0.8*dt;
            continue;
          }
          else
          {
            std::cout << "*** (ERROR) TimeInt_ModelH::outer_loop ***\n";
            std::cout << "Phi becomes unphysical at this timestep (dt = ";
            std::cout << dt << ")." << std::endl;
            exit(-1);
          }
        }

        // whole step for phi
        phi_whole = _Phi;
        step_phi( phi_whole, _Vel, mu_current, mu_current_linexp, mu_current_linimp, dt, t ); 
        is_physical = check_phi(phi_whole,_GlassMax);
        if ( !is_physical )
        {
          if (_VarDtFlag)
          {
            dt_new = 0.8*dt;
            continue;
          }
          else
          {
            std::cout << "*** (ERROR) TimeInt_ModelH::outer_loop ***\n";
            std::cout << "Phi becomes unphysical at this timestep (dt = ";
            std::cout << dt << ")." << std::endl;
            exit(-1);
          }
        }
      }

      // find the difference between the two stepping schemes
      phi_whole -= phi_half;

      // estimate truncation error, step doubling dt
      phi_err = 0.;
      for( UInt m=0; m<_NIcomp; ++m )
      {
        phi_err = std::max( phi_err, std::abs(phi_whole[m].max().real()) );
      }

      phi_iter++;
      dt_new = update_timestep( _Vel, phi_err, dt);

    }

    // Generate noise for fluctuating simulations
    // TODO: this method only works for PS algorithms. Either generalize it, or test for FD and quit.
    if(_EnableFluctuations)
    { //{{{
      // We need to reconstruct the semi-implicit stepping operators here.
      // TODO: The mobility and semi-implicit operators are recomputed here
      //       even though they were computed in step_phi.
      //       Review order of operations and optimize this away so that this task is not repeated.
      //
      // Get Mobility
      SmartFieldMat Mobility( _NIcomp, _NIcomp );
      _MobilityModel->set_mobility( _Phi, Mobility );

      // Cholesky decompose the mobility matrix
      SmartFieldMat L(_NIcomp, _NIcomp);
      // Begin by copying the upper triangle of M to L
      for(UInt i=0; i<_NIcomp; i++)
        for(UInt j=i; j<_NIcomp; j++)
          L(i,j) = Mobility(i,j);
      // Cholesky decomposition: the lower triangle of L will be replaced
      L.cholesky();

      // NOISE - generate one noise field per non-implicit phi field
      // This must be done for every field before adding noise to any field,
      // because any off-diagonal contributions in
      // the mobility matrix will mix the noise fields together
      SmartFieldVec noise_tmp(_NIcomp);
      SmartFieldVec noise(_NIcomp);
      SmartField    noisecoef;
      double dV = _CurrGrid->GetCell()->getVol() / _CurrGrid->NPW();
      // Fill a field with sqrt(2*dt/dV)*|k|
      noisecoef.setflag_inrealspace(false);
      noisecoef = -2.*dt/dV/sqrt(_EnergyModel->mu_scale())/_ReductionFactor/_ReductionFactor;//divided 2x due to sqrt later.
      _PhiOp->Del2_f_ex_inplace(noisecoef); // * (-k^2)
      noisecoef.sqrt(noisecoef);
      // Generate NIcomp noise fields with variance 2*k^2*dt/dV
      for(UInt i=0; i<_NIcomp; i++)
      {
        noise[i].setflag_inrealspace(true);
        noise[i].fillnoise(3+_CurrGrid->Dim());
        noise[i].fft_rtok();
        noise[i] *= noisecoef;
        noise[i].fft_ktor();
      }

      // Couple the noise fields together with the Cholesky decomposed mobility matrix
      for(UInt i=0; i<_NIcomp; i++)
      {
        noise_tmp[i] = noise[0];
        noise_tmp[i] *= L(i,0);
        for(UInt j=1; j<=i; j++)
          noise_tmp[i].accumulateproduct_inplace(noise[j],L(i,j));
      }
      // We need to solve the linear system for the noise application (in k space).
      // Note also that the local approximation we used to treat the mobility
      // is not conserving.  Enforce conservation manually
      noise_tmp.fft();
      for(UInt i=0; i<_NIcomp; i++)
        noise_tmp[i].setaverage(0.);


      // -----------------------------
      // Find lhs (all implicit terms)
      // lhs = I - dt Mij_max del^2 mu_lin_op_j
      // -----------------------------
      Vec2dFieldType M_max = _MobilityModel->get_mobility_max( Mobility );
      SmartFieldOpMat lhs( _NIcomp, _NIcomp );
      lhs.zero();
      lhs.setflag_inrealspace( false );
      // I
      _PhiOp->Eye_im( lhs );
      // implicit, linear part of mu
      SmartFieldOpMat tmpop( _NIcomp, _NIcomp );
      tmpop.zero();
      tmpop.setflag_inrealspace( false );
      _PhiOp->A_Del2_F_im( M_max, mu_current_linimp, tmpop ); // tmpop <- M_max del^2 mu_lin_np1
      tmpop *= (-dt);
      lhs += tmpop;

      // -----------------------------
      // Solve for phi_n+1
      // lhs (phi_n+1) = rhs
      //    lhs = implicit ops
      //    rhs = noise
      //    (phi_n+1) must be accumulated into phi - > place into noisecoefs for convenience
      // -----------------------------
      // noise_tmp is in Fourier space. Solve with the implicit operator.
      _PhiOp->Solve_im_matrix( lhs, noise_tmp, noise, t );

      // Add the noise to phi in real space
      noise.ifft();
      phi_half += noise;
    } // }}}

    // --- time step ---
    _Phi = phi_half; // phi_half has taken 2x dt/2 steps, so this is a whole step

    // step the final vel
    _Vel = vel_slope;
    _Vel.axpy_inplace(vel_half, 0.5*dt);

    // Chemical potential needs to be refreshed for the two-half-step phi
    phi_half.fft(); // phi_half will be discarded - used to calc_mu with correct Fourier rep.
    _EnergyModel->calc_mu(phi_half, mu_half);

    // Update _Vel
    step_vel( _Phi, _Vel, mu_half, vel_iter, t );

    // Update vel_slope for future continuations
    vel_slope = _Vel;
    vel_slope.axpby_inplace(vel_half, 2./dt, -2./dt);

    t = t + dt;
    Nt++;

    bool enableSK(false);
    if(enableSK)
    {
      numSKsamples++;
      double V(_CurrGrid->GetCell()->getVol());
      SmartField tmp;
      for(int i=0; i<_NIcomp; i++)
        _Phi[i].fft_rtok();
      // S(k) = <rho_k rho_{-k}> = <|rho_k|^2>
      int k=0;
      for(int i=0; i<_NIcomp; i++)
        for(int j=i; j<_NIcomp; j++, k++)
        {
          tmp = _Phi[i];
          tmp.prodconjg(_Phi[j]);
          SK[k].xpby_inplace(tmp,V);
        }
      // Revert phi fields to (r) space
      for(int i=0; i<_NIcomp; i++)
        _Phi[i].fft_ktor();
    }

    // --- write diagnostics and data ---
    // step doubling diagnostics
    if ( Nt % _NStepDisp == 0 || t >= t_check || t == _TMax || Nt == _NtMax )
    {
      stepfile << Nt << "\t" << t << "\t" << dt << "\t";
      stepfile << phi_err << "\t" << phi_iter << "\t" << vel_iter << std::endl;
      printf("%8lu%12.3e%12.3e%12.3e%12lu%12lu", Nt, t, dt, phi_err, phi_iter, vel_iter);
      std::cout << std::endl;

      // Average and extrema compositions---useful for debugging
      // Note: organzied this way so that the average compositions are the left-most cols. 
      compfile << Nt << "\t";
      for(int i=0; i<_NIcomp; i++)
        compfile << _Phi[i].integrate().real() << "\t";
      for(int i=0; i<_NIcomp; i++)
        compfile << _Phi[i].minsigned().real() << "\t";
      for(int i=0; i<_NIcomp; i++)
        compfile << _Phi[i].maxsigned().real() << "\t";
      compfile << std::endl;
    }

    // write data to file
    if ( t >= t_check || t == _TMax || Nt == _NtMax )
    {
      printf("      * writing ndisp = %05lu to disk\n", n_disp);
      write_curr_state( n_disp );
      timefile << n_disp << "\t" << t << "\t" << Nt << std::endl;
      n_disp++;

      if ( _OutputIntervalFlag == 0 )
      {
        t_check = n_disp*_DtDisp;
      }
      else
      {
        t_check = pow(10., n_disp*_DtDisp); // minimum output is at t = 1
      }
    }

    // change dt_new so that output is on the exact time we want
    if ( (t+dt_new) > t_check ) 
    {
      dt_new = t_check-t;
      // if dt < _DtMin, it will crash
      dt_new = std::max(dt_new, _DtMin); 
    }

  }

  stepfile.close();
  timefile.close();

  std::cout << std::endl;
  std::cout << " * Finished time integration\n";
  std::cout << "   - Final Time: " << t << std::endl;
  std::cout << "   - No. of Steps: " << Nt << std::endl;

} // }}}

void TimeInt_ModelH :: step_phi( SmartFieldVec &phi, 
                                 SmartFieldVec const &vel, 
                                 SmartFieldVec const &mu, 
                                 SmartFieldVec const &mu_lin, 
                                 SmartFieldOpMat const &mu_lin_np1,
                                 RealType dt,
                                 RealType t )
{ // {{{

  // Solves the following (1 step of Model H, convection-diffusion):
  //
  // ( I - dt Mij_max del^2 mu_lin_op_j  ) phi_n+1 = 
  //     phi_n - dt vel . del phi_n + 
  //     dt del . (Mij del . mu) - dt*Mij_max del^2 mu_lin_j
  //

  // Get Mobility
  //SmartFieldMat Mobility( _NIcomp, _NIcomp );
  _MobilityModel->set_mobility( phi, _Mobility );

  // Get maximum mobility
  Vec2dFieldType M_max = _MobilityModel->get_mobility_max( _Mobility );

  // -----------------------------
  // Find rhs (all explicit terms)
  // rhs = phi_n -
  //       dt vel . del phi_n +
  //       dt del . (Mij del . mu) -
  //       dt Mij_max del^2 mu_lin_j
  // -----------------------------

  SmartFieldVec rhs( _NIcomp );

  // phi_n
  phi.fft();
  rhs = phi;

  if (_AdvFlag)
  { // advection term (calculated as vel. grad-phi)
    // standard formulation in Tree et al, Soft Matter 2017
    SmartFieldMat grad_phi_n( _NIcomp, _Dim );
    _PhiOp->Del_f_ex( phi, grad_phi_n );
    grad_phi_n.ifft();
    _tmpvec.setflag_inrealspace( true );
    // tmpvec <- vel . (del phi_n)
    _tmpvec.dot( grad_phi_n, vel );
    _tmpvec.fft();
   
    //added for debugging
    //RealType avg;
    //for(UInt i=0; i<_NIcomp; i++)
    //{
     // avg = _tmpvec[i].getaverage().real();
      //std::cout << "Avg advection for component " << i << " = " << avg << std::endl; 
    //}
    rhs.xpby_inplace(_tmpvec, -dt);
  }
  else
  { // advection term (calculated as div (phi*vel))
    // alternative formulation to solve composition drift with fluctuations ON
    phi.ifft();
    for(UInt i=0; i< _NIcomp; i++)
    {
      //Real-space calculations
      _tmpvec = vel;
      _tmpvec *= phi[i];
      
      //_tmpvec2[i] <-- div(phi[i]*v) (k-space)
      _tmpvec.fft();
      _PhiOp->Div_f_ex( _tmpvec, _tmpvec2[i] );
    }
    phi.fft();
    
    //added for debugging
    //RealType avg;
    //for(UInt i=0; i<_NIcomp; i++)
    //{
      //avg = _tmpvec2[i].getaverage().real();
      //std::cout << "Avg advection term for component " << i << " = " << avg << std::endl; 
    //}
    
    rhs.xpby_inplace(_tmpvec2, -dt);
  
  }

  // get full mu term
  //_tmpvec.setflag_inrealspace( false );
  _PhiOp->Del_A_Del_f_ex( mu, _Mobility, _tmpvec ); // tmpvec <- del.(Mij del. mu_n)
  rhs.xpby_inplace(_tmpvec, dt);

  // get explicit, linear term of mu
  _PhiOp->Del2_f_ex( mu_lin, _tmpvec ); // tmpvec <- del^2 mu_lin
  _tmpvec2.dot(M_max, _tmpvec);             // mu_lin <- M_max . tmpvec
  rhs.xpby_inplace( _tmpvec2, (-1.*dt) );
  
  // add full noise term only in constant time-stepping
  // for adaptive time-stepping, noise is stepped outside of step_phi for performance
  // note: this is based on code design from Model B, still need to optimize for Model H
  if(_EnableFluctuations && !(_VarDtFlag))
  {
    calc_noise_phi ( _Mobility, _ReductionFactor, dt, _tmpvec );// tmp <- noise
    rhs.xpby_inplace( _tmpvec, 1.0); //effect of dt already calculated in previous step
    //added for debugging
    //RealType avg;
    //for(UInt i=0; i<_NIcomp; i++)
    //{
    //  avg = _tmpvec[i].getaverage().real();
    //  std::cout << "Avg noise for component " << i << " = " << avg << std::endl; 
    //}
  }

  // -----------------------------
  // Find lhs (all implicit terms)
  // lhs = I - dt Mij_max del^2 mu_lin_op_j
  // -----------------------------
  SmartFieldOpMat lhs( _NIcomp, _NIcomp );
  lhs.setflag_inrealspace( false );
  lhs.zero();

  // I
  _PhiOp->Eye_im( lhs );

  // implicit, linear part of mu
  SmartFieldOpMat tmpop( _NIcomp, _NIcomp );
  tmpop.setflag_inrealspace( false );
  tmpop.zero();
  _PhiOp->A_Del2_F_im( M_max, mu_lin_np1, tmpop ); // tmpop <- M_max del^2 mu_lin_np1
  tmpop *= (-dt);
  lhs += tmpop;

  // get reaction terms if required
  if (_ReactionsModel->GetReactionsFlag() != 0)
  {
    _ReactionsModel->calc_reactions( phi, _tmpvec );
    rhs.xpby_inplace(_tmpvec, dt);
    _ReactionsModel->calc_reactions_lin_exp( phi, _tmpvec );
    rhs.xpby_inplace( _tmpvec, (-1.*dt) );
    _ReactionsModel->calc_reactions_lin_imp( phi, tmpop );
    tmpop *= (-dt);
    lhs += tmpop;
  }

  // -----------------------------
  // Solve for phi_n+1
  // lhs (phi_n+1) = rhs
  // -----------------------------

  _PhiOp->Solve_im_matrix( lhs, rhs, phi, t );
  phi.ifft();
  phi.zeroimag();

  if ( (_CapPhiFlag == 1) && (phi[0].maxsigned().real() > _GlassMax) )
  {
    capField(phi[0], _GlassMax);
  }

  return;

} // }}}

void TimeInt_ModelH :: step_vel( SmartFieldVec const &phi, 
                                 SmartFieldVec &vel, 
                                 SmartFieldVec const &mu,
                                 UInt &iter,
                                 RealType t )
{ // {{{  

  SmartFieldVec div_Pi( _Dim ); // divergence of osmotic stress
  SmartField eta; // viscosity
  eta.setflag_inrealspace(true);
  eta.zero();

  // --- get eta and eta_max ---
  _ViscModel->set_eta( phi, eta );
  FieldType eta_max = std::abs(eta.max());
  eta += -eta_max; // eta -> delta_eta

  // --- get divergence of osmotic stress ---
  // (Don't include it in iteration; 
  //  it is not a function of velocty)
  // del . Pi = Nr sum_i phi_i del mu_i
  {
    SmartFieldMat del_mu( _NIcomp, _Dim ); // grad of chem pot.
    del_mu.setflag_inrealspace(false);

    // take gradient of mu
    _VelOp->Del_f_ex( mu, del_mu);

    // div_Pi <- sum_i phi_i del mu_i (dot product over NIcomp)
    del_mu.ifft(); // needs to be in real space for product
    div_Pi.setflag_inrealspace(true);
    div_Pi.dot( phi, del_mu );

    // div_Pi *= Nr
    div_Pi *= FieldType( _EnergyModel->mu_scale() );

    // add body force if required
    
    if ( _BodyForceFlag )
    {
      if (_BodyForceFlag == 2 )
      {
        f_sigmoid(phi,_BodyForce,bf_comp,_PhiT,_W,_Scale,finvert); 
      }
      div_Pi += _BodyForce;
    }

    // fluctuations if required 
    if (_EnableFluctuations)
    {
      SmartField eta_tmp;//need full eta, not delta, as set earlier
      _ViscModel->set_eta( phi, eta_tmp );
      SmartFieldVec noise ( _Dim ); 
      calc_noise_vel(eta_tmp, _ReductionFactor, noise);
      noise.ifft();

      div_Pi += noise;
    }
   
    // take Fourier transform
    div_Pi.fft();

  }

  // --- get velocity ---
  // eta_m del^4 v = (del^2 - del del) . div_Pi

  if (_VarViscFlag) // variable viscosity
  { //{{{

    UInt no = 0; // outer loop counter for iteration
    UInt n = 0; // counter for iteration
    RealType vel_err = 1.;
    UInt velo_iter_max = 1000; // max # of iterations for outer loop
    UInt vel_iter_max = 1500;

    SmartFieldVec VarVisc( _Dim ); // variable viscosity term
    SmartFieldVec rhs( _Dim );
    SmartFieldVec vel_star( _Dim );
    SmartFieldVec tmp(_Dim);
    SmartFieldVec vslope(_Dim); // velocity slope to improve inner loop guess
    vslope = FieldType(0.);
    vslope.setflag_inrealspace ( true );
    SmartFieldVec vo(_Dim); // velocity fed to inner loop

    // lhs = del^2
    SmartFieldOpMat lhs(_Dim, _Dim);
    lhs.zero();
    lhs.setflag_inrealspace( false );
    Vec2dFieldType Eye( eye(_Dim, _Dim, 1.) );
    _VelOp->A_Del2_im(Eye, lhs); // lhs = del^2

    SmartFieldVec v_n( _Dim ); // original velocity
    SmartFieldVec g_n( _Dim ); // G(v_n)
    SmartFieldVec g_nm1( _Dim ); // G(v_nm1)
    SmartFieldVec g_nm2( _Dim ); // G(v_nm1)
    SmartFieldVec d_n( _Dim ); // G(v_n) - v_n
    SmartFieldVec d_nm1( _Dim ); // G(v_nm1) - v_nm1
    SmartFieldVec d_nm2( _Dim ); // G(v_nm2) - v_nm2
    SmartFieldVec d_01( _Dim ); // d_n - d_nm1
    SmartFieldVec d_02( _Dim ); // d_n - d_nm2

    SmartFieldMat grad_vel( _Dim, _Dim ); // velocity graident tensor
    SmartFieldMat grad_vel_T( _Dim, _Dim ); // transpose of grad v

    FieldType a11, a12, a21, a22, b1, b2, det_a, c1, c2;

    // do Picard iteration of v = G(v)
    while (vel_err > _VelErrTol && no < velo_iter_max)
    {
      vo   = vel;
      vel  = vo;
      vel += vslope;
      while (vel_err > _VelErrTol && n < vel_iter_max)
      {
        //{{{
        v_n = vel;
        g_n = vel;

        // Variable Viscosity Term
        g_n.fft();
        grad_vel.setflag_inrealspace( false );
        _VelOp->Del_f_ex( g_n, grad_vel ); // gradient of velocity (d/dx_i v_j)

        grad_vel_T = grad_vel;
        grad_vel_T.transpose(); // transpose of gradient (d/dx_j v_i)

        grad_vel += grad_vel_T; // del v + del v.T (2*symmetric part)
        grad_vel.ifft();
        grad_vel *= eta;  // eta * (dvj/dxi + dvi/dxj)
        grad_vel.fft();
        
        VarVisc.setflag_inrealspace( false );
        _VelOp->Div_f_ex( grad_vel, VarVisc ); // d/dx_{i} . (eta * Sij)

        // get total forcing
        // rhs: (del^2 - del del) . (div_Pi - VarVisc)
        VarVisc *= RealType(-1.);
        VarVisc += div_Pi;

        // solve eta del^2 vel^{*} = forcing
        vel_star.zero();
        vel_star.setflag_inrealspace( false );
        _VelOp->Solve_im_matrix( lhs, VarVisc, vel_star, t );

        // If periodic, need to zero out the 
        // zero-k mode, since it is undefined
        if (_VelBCs->BCFlag() == 0 )
        {
          vel_star.zero_k0();
        }
        vel_star *= (1./eta_max);

        // now get divergence free part of vel_star
        // del^2 g_n = (del^2 I - del del ) . vel_star
        // del^2 g_n = rhs
        rhs.zero();
        rhs.setflag_inrealspace( false );
        _VelOp -> transverse_f_ex( vel_star, rhs );
        _VelOp->Solve_im_matrix( lhs, rhs, g_n, t);

        // If periodic, need to zero out the 
        // zero-k mode, since it is undefined
        if (_VelBCs->BCFlag() == 0 )
        {
          g_n.zero_k0();
        }
        g_n.ifft();
        g_n.zeroimag();

        // update velocity using Anderson mixing to accelerate convergence
        d_n = g_n; d_n -= v_n; // d_n = g_n - v_n

        if (n%3 == 0 and n >= 2)
        {

          // -- Anderson mixing --
          // find minimum coeffs from:
          // [a11, a12] [c0] = [b1] 
          // [a21, a22] [c1] = [b2]
          // where, 
          // a11 = d_01 . d_01
          // a12 = d_01 . d_02
          // a21 = d_02 . d_01
          // a22 = d_02 . d_02
          // b1 = d_n . d_01
          // b2 = d_n . d_02
          // solve by Cramers rule:
          // c1 = (b1*a22 - b2*a12)/det(a)
          // c2 = (a11*b2 - b1*a21)/det(a)
          // det(a) = a11*a22 - a21*a12

          d_01 = d_n; d_01 -= d_nm1; // d_01 = d_n - d_nm1
          d_02 = d_n; d_02 -= d_nm2; // d_02 = d_n - d_nm2

          a11 = 0.;
          tmp = d_01; tmp *= d_01;
          for (UInt i=0; i<_Dim; i++) 
            a11 += tmp[i].sumelem();

          a12 = 0.;
          tmp = d_01; tmp *= d_02;
          for (UInt i=0; i<_Dim; i++) 
            a12 += tmp[i].sumelem();

          a21 = a12;

          a22 = 0.;
          tmp = d_02; tmp *= d_02;
          for (UInt i=0; i<_Dim; i++) 
            a22 += tmp[i].sumelem();

          b1 = 0.;
          tmp = d_n; tmp *= d_01;
          for (UInt i=0; i<_Dim; i++) 
            b1 += tmp[i].sumelem();

          b2 = 0.;
          tmp = d_n; tmp *= d_02;
          for (UInt i=0; i<_Dim; i++) 
            b2 += tmp[i].sumelem();

          det_a = a11*a22-a12*a21;
          if (det_a.real() == 0) 
          {
            std::cout << "Error: Singular matrix in Anderson mixing." << std::endl;
            exit(1);
          }
          c1 = (b1*a22 - b2*a12)/det_a;
          c2 = (b2*a11 - b1*a21)/det_a;
          
        }
        else
        {
          c1 = FieldType(0.);
          c2 = FieldType(0.);
        }

        // vel = (1-c1-c2)*g_n + c1*g_nm1 + c2*g_nm2
        vel = g_n; vel *= FieldType(1.)-c1-c2;
        if (n>0)
        {
          tmp = g_nm1; tmp *= c1;
          vel += tmp;
        }
        if (n>1)
        {
          tmp = g_nm2; tmp *= c2;
          vel += tmp;
        }

        // find error, do we need to iterate?
        v_n -= vel;
        vel_err = v_n.maxl2norm()/vel.maxl2norm();
      
        // get history for anderson mixing
        if (n > 0)
        {
          g_nm2 = g_nm1;
          d_nm2 = d_nm1;
        }
        g_nm1 = g_n;
        d_nm1 = d_n;

        // increment the iteration counter
        n++;
        //}}}
      } // end inner while
      vslope  = vel;
      vslope -= vo;
      if (vel_err < _VelErrTol)
      {
        no *= vel_iter_max;
        iter = no;
        iter +=n;
      }
      no++;
      n = 0;
    } //end of outer while

    //iter = n;

    if (vel_err > _VelErrTol)
    {
      std::cout << "(WARNING) velocity did not converge to tolerance: ";
      std::cout << vel_err << " steps: " << n << std::endl;
    }

  } // }}}
  else // constant viscosity
  {

    SmartFieldVec vel_star( _Dim );
    vel_star.zero();
    vel_star.setflag_inrealspace( false );
    
    // lhs: del^2
    SmartFieldOpMat lhs(_Dim, _Dim);
    lhs.zero();
    lhs.setflag_inrealspace( false );
    Vec2dFieldType Eye( eye(_Dim, _Dim, 1.) );
    _VelOp->A_Del2_im(Eye, lhs);
  
    // solve eta del^2 vel^{*} = div . Pi
    _VelOp->Solve_im_matrix( lhs, div_Pi, vel_star, t );

    // If periodic, need to zero out the 
    // zero-k mode, since it is undefined
    if (_VelBCs->BCFlag() == 0 )
    {
      vel_star.zero_k0();
    }
    vel_star *= (1./eta_max);

    // now get divergence free part of vel_star
    // del^2 vel = (del^2 I - del del ) . vel_star
    // del^2 vel = rhs
    SmartFieldVec rhs( _Dim );
    rhs.setflag_inrealspace( false );
    _VelOp -> transverse_f_ex( vel_star, rhs );
    vel.setflag_inrealspace( false );
    _VelOp->Solve_im_matrix( lhs, rhs, vel, t );

    // If periodic, need to zero out the 
    // zero-k mode, since it is undefined
    if (_VelBCs->BCFlag() == 0 )
    {
      vel.zero_k0();
    }

    vel.ifft();
    vel.zeroimag();
    iter = 0;
  }

} // }}}

RealType TimeInt_ModelH :: update_timestep( SmartFieldVec & vel, 
                                            RealType phi_err, 
                                            RealType dt )
{ // {{{

  static bool prev_step_accept = true;
  RealType theta = 0.95; // safety factor so we don't reject steps if possible
  RealType dt_new = 0.;   // new step size

  RealType dt_cfl = 1.; // cfl time step size
  RealType vmax; // maximum of absolute value of fields

  // estimate CFL stability criteria for convection
  // dt_cfl = { sum [ max(abs(v_i))/dx_i ] }^{-1}
  dt_cfl = 0.;
  for( UInt m=0; m<_Dim; ++m )
  { 
    //dx = (_CurrLayout->getCell().getBoxTensor()[m][m])/(_CurrLayout->getGridDim(m));
    vmax = std::abs( vel[m].max() );
    dt_cfl += vmax/_CurrGrid->DxDim(m);
  }
  dt_cfl = 1./dt_cfl;
  
  // only step bigger if (1) big enough and (2) didn't reject prev step; always step smaller
  if ( (prev_step_accept && ((dt_new - dt)/dt > 0.01)) || dt_new < dt )
  {
    // dt_step_dbl = dt * sqrt{ _PhiErrTol/ max[abs(phi_err)] }
    dt_new = dt*sqrt( theta*_PhiErrTol/phi_err );
    dt_new = std::min(dt_new, dt_cfl);
    dt_new = std::min(dt_new, _DtMax);
    dt_new = std::max(dt_new, _DtMin);
  }
  else
  {
    dt_new = dt;
  }
  
  if( phi_err <= _PhiErrTol )
  {
    prev_step_accept = true;
  }
  else
  {
    prev_step_accept = false;
  }

  return dt_new;
  
} // }}}

void TimeInt_ModelH :: f_sigmoid(const  SmartFieldVec &phi, SmartFieldVec &BodyForce, int comp, RealType phiT, RealType w, Vec1dReal scale,SmartField &finvert)
{ // {{{

// (in) phi: incoming concentration
// (in) phiT: threshold
// (in) w: characteristic width of transition in phi
// (out) phi: sigmoid
//  
// f(phi) = 1/(1 + exp((phi-phiT)/w))

  SmartField tmp;
  tmp = phi[comp];
  tmp += -phiT;
  tmp.exponentiate(tmp, -1./w);
  tmp += 1.; 
  for(int i=0;i<_Dim;i++)
  {
    BodyForce[i] = finvert;
    BodyForce[i] /= tmp;
    BodyForce[i] *= scale[i];
//    BodyForce[i] *= finvert;
  }
} // }}}

