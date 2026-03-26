/***************************************************************
*
* Class for the time integration scheme for Model B
*   - semi-implicit, Euler scheme
*   - variable time-step (step-doubling)
*
* DRT -- Fri, 13 Mar 2015
*
****************************************************************/
//
//
// d phi_i / d t = sum_j del . ( M_ij del mu_j )
// 
//

#include "TimeInt_ModelB.h"

// ----------------- Constructor/Destructor -----------------

TimeInt_ModelB :: TimeInt_ModelB( int Ncomp, int NIcomp, Grid* CurrGrid ):
TimeInt_Base( Ncomp, NIcomp, CurrGrid ),
_EnableFluctuations(false),
_Phi( NIcomp )
// {{{
{

  std::cout << std::endl;
  std::cout << " * Time Integration derived class initialization\n";
  std::cout << "   - Model B (diffusion only)" << std:: endl;

} // }}}

TimeInt_ModelB :: ~TimeInt_ModelB()
{ // {{{
} // }}}

// ----------------- Input/Output -----------------

void TimeInt_ModelB :: set_params()
{//{{{

  
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

      _OutputIntervalFlag = d["output_interval_flag"].GetInt();

      _NtMax = d["Nt_max"].GetInt();

      _NStepDisp = d["N_step_disp"].GetInt();

      _PhiIterMax = d["phi_iter_max"].GetInt();

      _VarDtFlag = d["var_dt_flag"].GetInt();

      _CapPhi0Flag = d["CapPhi0Flag"].GetInt();
      
      _Phi0Max = d["phi0_max"].GetDouble();

      _PhiErrTolFlag = d["phi_err_tol_flag"].GetInt();
      
      _PhiErrTolMax = d["phi_err_tol_max"].GetDouble();
      
      _PhiErrTolMin = d["phi_err_tol_min"].GetDouble();
      
      _PhiErrTolW = d["phi_err_tol_w"].GetDouble();

      _PhiErrTolMid = d["phi_err_tol_mid"].GetDouble();

      // Read flag for fluctuating simulations.
      // Tolerate old input files not containing the flag.
     
      int fluct = d["enablefluctuations"].GetInt();
      _EnableFluctuations = (fluct == 1) ? true:false;	
      
      _ReductionFactor = d["reductionFactor"].GetDouble();

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
      _CapPhi0Flag = atoi(tmp_str.c_str());
      std::getline(f1, tmp_str, '\n');
      
      std::getline(f1, tmp_str, '#');
      _Phi0Max = atof(tmp_str.c_str());
      std::getline(f1, tmp_str, '\n');

      std::getline(f1, tmp_str, '#');
      _PhiErrTolFlag = atoi(tmp_str.c_str());
      std::getline(f1, tmp_str, '\n');
      
      std::getline(f1, tmp_str, '#');
      _PhiErrTolMax = atof(tmp_str.c_str());
      std::getline(f1, tmp_str, '\n');
    
      std::getline(f1, tmp_str, '#');
      _PhiErrTolMin = atof(tmp_str.c_str());
      std::getline(f1, tmp_str, '\n');
      
      std::getline(f1, tmp_str, '#');
      _PhiErrTolW = atof(tmp_str.c_str());
      std::getline(f1, tmp_str, '\n');

      std::getline(f1, tmp_str, '#');
      _PhiErrTolMid = atof(tmp_str.c_str());
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

  std::cout << "   - Parameters:\n" << std::scientific;
  std::cout << "      TimeIntFlag        = " << _TimeIntFlag << std::endl;
  std::cout << "      TMax               = " << _TMax << std::endl;
  std::cout << "      Dt0                = " << _Dt0 << std::endl;
  std::cout << "      DtMax              = " << _DtMax << std::endl;
  std::cout << "      DtMin              = " << _DtMin << std::endl;
  std::cout << "      DtDisp             = " << _DtDisp << std::endl;
  std::cout << "      OutputIntervalFlag = " << _OutputIntervalFlag << std::endl;
  std::cout << "      NtMax              = " << _NtMax << std::endl;
  std::cout << "      NStepDisp          = " << _NStepDisp << std::endl;
  std::cout << "      PhiIterMax         = " << _PhiIterMax << std::endl;
  std::cout << "      VarDtFlag          = " << _VarDtFlag << std::endl;
  std::cout << "      Phi0Max            = " << _Phi0Max << std::endl;
  std::cout << "      PhiErrTolFlag      = " << _PhiErrTolFlag << std::endl;
  std::cout << "      PhiErrTolMax       = " << _PhiErrTolMax  << std::endl;
  std::cout << "      PhiErrTolMin       = " << _PhiErrTolMin << std::endl;
  std::cout << "      PhiErrTolW         = " << _PhiErrTolW  << std::endl;
  std::cout << "      PhiErrTolMid       = " << _PhiErrTolMid << std::endl;
  std::cout << "      EnableFluctuations = " << _EnableFluctuations << std::endl;
  std::cout << "      ReductionFactor    = " << _ReductionFactor << std::endl;

  return;

} // }}}

void TimeInt_ModelB::read_initial_state()
{ // {{{

  // read in initial conditions
  char filename[40];
  for (int i=0; i<_NIcomp; i++)
  {
    sprintf( filename, "phi%d.in", i+1 );
    _Phi[i].readfield(filename);
  }

} // }}}

void TimeInt_ModelB::write_final_state()
{ // {{{

  // write out final state
  char filename[40];
  for (int i=0; i<_NIcomp; i++)
  {
    sprintf( filename, "phi%d.out", i+1 );
    _Phi[i].writerealfield(filename);
  }

} // }}}

void TimeInt_ModelB :: write_curr_state( int n_disp )
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

} // }}}

// ----------------- Integration Routines -----------------
void TimeInt_ModelB :: outer_loop()
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

void TimeInt_ModelB :: outer_loop_const()
{ // {{{
  // constant time-stepper

  // Parameters here were  based on adaptive time-stepping
  // Could clean some parameters for constant time-stepping
  RealType t = 0.; // total time
  RealType dt = _Dt0; // time step
  RealType dt_new = _Dt0; // time step
  RealType t_check = 0; // time when write to disk
  bool is_physical; // make sure phi is between zero and one
  UInt Nt = 0; // number of timesteps taken
  UInt n_disp = 0; // number of fields written to file

  // record which times are written
  std::ofstream timefile;
  timefile.open("time.dat");
  timefile << "# Columns: Ndisp time Nt\n" << std::scientific;

  // record compositions for sanity checks and debugging
  std::ofstream compfile;
  compfile.open("compositions.dat");
  compfile << std::scientific;
  compfile << "# Columns: N phi[0].avg ... phi[N].avg phi[0].min ... phi[N].min phi[0].max ... phi[N].max\n";
  
  // print to screen (or log file)
  // header of time integration entries
  std::cout << std::endl;
  std::cout << " * Beginning Time Integration Loop\n";
  std::cout << " * Constant time stepping selected. No error estimates provided.\n";
  std::cout << " * Time stepping is adjusted only for printing out data at requested times.\n";
  printf("%8s", "Nt");
  printf("%12s", "t");
  printf("%12s", "dt");
  std::cout << std::endl;
  
  // print out initial composition
  printf("      * writing ndisp = %05lu to disk\n", n_disp);
  write_curr_state( n_disp );
  timefile << n_disp << "\t" << t << "\t" << Nt << std::endl;
  n_disp++;
  
  // set t_check
  if ( _OutputIntervalFlag == 0 )
  {
    t_check = n_disp*_DtDisp;
  }
  else
  {
    t_check = pow(10., n_disp*_DtDisp); // minimum output is at t = 1
  }

  // Thermally averaged structure factor for all pairs
  SmartFieldVec SK((_NIcomp*(_NIcomp+1))/2); // Species-species structure factors in a packed storage format
  UInt numSKsamples(0);
  histogram<FieldType, double> SKhist(*_CurrGrid->GetLayout(), true);
  SK.setflag_inrealspace(false);
  SK.zero();

  // Time integration loop
  while( t < _TMax && Nt < _NtMax )
  {
    
    SmartFieldVec mu_current(_NIcomp);
    SmartFieldVec mu_current_linexp(_NIcomp);
    SmartFieldOpMat mu_current_linimp(_NIcomp,_NIcomp);
    mu_current.setflag_inrealspace(false);

    // Compute arguments for step_phi
    _Phi.fft();
    _EnergyModel->calc_mu(_Phi, mu_current);
    _EnergyModel->calc_mu_lin_exp(_Phi, mu_current_linexp);
    _EnergyModel->calc_mu_lin_imp(_Phi, mu_current_linimp);
    _Phi.ifft();
    
    // Step phi (noise addition done within function) 
    step_phi( _Phi, mu_current, mu_current_linexp, mu_current_linimp, dt, t );
    
    // Check if phi is still within  physical bounds after solution
    is_physical = check_phi(_Phi,_Phi0Max);
    if ( !is_physical )
    {
      std::cout << "Phi becomes unphysical at this timestep ";
      std::cout << "(dt = " << dt << ")." << std::endl;
      exit(-1);
    }

    // --- time step ---
    t = t + dt;
    Nt++;

    bool enableSK(false);
    if(enableSK)
    { // {{{
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
          //tmp *= _CurrGrid->GetCell()->getVol();
          //SK[k] += tmp;
          SK[k].xpby_inplace(tmp,V);
        }
      // Revert phi fields to (r) space
      for(int i=0; i<_NIcomp; i++)
        _Phi[i].fft_ktor();
    } // }}}

    // --- write diagnostics and data ---
    if ( Nt % _NStepDisp == 0 || t >= t_check || t == _TMax || Nt == _NtMax )
    {
      //screen (log file) reports
      printf("%8lu%12.3e%12.3e", Nt, t, dt);
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

      if(enableSK)
      {
        std::stringstream fname;
        int k=0;
        for(int i=0; i<_NIcomp; i++)
          for(int j=i; j<_NIcomp; j++, k++)
          {
            SKhist.mapfromfield(SK[k].GetField(0), numSKsamples);
            fname.str("");
            fname<<"SK"<<i<<j<<".dat";
            SKhist.writehistogram(fname.str().c_str());
          }
      }
    }

    // write data to file
    if ( t - t_check >= -1e-14 || t == _TMax || Nt == _NtMax )
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

  timefile.close();

  std::cout << std::endl;
  std::cout << " * Finished time integration\n";
  std::cout << "   - Final Time: " << t << std::endl;
  std::cout << "   - No. of Steps: " << Nt << std::endl;
	
} // }}}

void TimeInt_ModelB :: outer_loop_adap()
{ // {{{
  // adaptive time-stepper
  // step-doubling algorithm by DRT	

  SmartFieldVec phi_half( _Phi );
  SmartFieldVec phi_whole( _Phi );

  RealType t = 0.; // total time
  RealType dt = _Dt0; // time step
  RealType dt_new = _Dt0; // time step
  RealType t_check = 0; // time when write to disk
  RealType phi_err = 1.; // max truncation error
  bool is_physical; // make sure phi is between zero and one
  UInt Nt = 0; // number of timesteps taken
  UInt n_disp = 0; // number of fields written to file
  UInt phi_iter = 0; // number of iterations for phi

  // record which times are written
  std::ofstream timefile;
  timefile.open("time.dat");
  timefile << "# Columns: Ndisp time Nt\n" << std::scientific;

  // diagnostics of step-doubling
  std::ofstream stepfile; 
  stepfile.open("stepsize.dat");
  stepfile << "# Columns: N t dt phi_err phi_iter phi_max phi_err_tol\n";
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
  printf("%12s", "phi_max");
  printf("%12s", "phi_err_tol");
  std::cout << std::endl;
  
  _PhiErrTol = _PhiErrTolMax;// initialize error tolerance with max tolerance

  // write initial state
  stepfile << Nt << "\t" << t << "\t" << dt << "\t";
  stepfile << phi_err << "\t" << phi_iter << "\t" << _PhiErrTol << std::endl;
  printf("%8lu%12.3e%12.3e%12.3e%12lu%12.3e", Nt, t, dt, phi_err, phi_iter, _PhiErrTol);
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

    if ( _PhiErrTolFlag > 0 )
    {
      _PhiErrTol =  _PhiErrTolMin + (_PhiErrTolMax - _PhiErrTolMin)/2.*(
        std::tanh( (_Phi0Max - _Phi[0].maxsigned().real() - _PhiErrTolMid)/_PhiErrTolW ) + 1);
    }

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
          std::cout << "(ERROR) Insufficient progress: ";
          std::cout << "timestep below _DtMin. Aborting Calculation.";
          std::cout << std::endl;
          write_curr_state(-1);
          exit(1);
        }

        // initialize for step doubling
        phi_half = _Phi;
        phi_whole = _Phi;

        // first half step for phi
        step_phi( phi_half, mu_current, mu_current_linexp, mu_current_linimp, dt/2., t );
        is_physical = check_phi(phi_half, _Phi0Max);
        if ( !is_physical )
        {
          if (_VarDtFlag)
          {
            dt_new = 0.8*dt;
            std::cout << "Phi unphysical in first half step." << std::endl;
            std::cout << "(dt = " << dt << ")." << std::endl;
            continue;
          }
          else
          {
            std::cout << "Phi becomes unphysical at this timestep ";
            std::cout << "(dt = " << dt << ")." << std::endl;
            exit(-1);
          }
        }

        // Compute and store mu at half step
        phi_half.fft();
        _EnergyModel->calc_mu(phi_half, mu_half);
        _EnergyModel->calc_mu_lin_exp(phi_half, mu_half_linexp);
        _EnergyModel->calc_mu_lin_imp(phi_half, mu_half_linimp);
        phi_half.ifft();

        // second half step for phi
        step_phi( phi_half, mu_half, mu_half_linexp, mu_half_linimp, dt/2., t );
        is_physical = check_phi(phi_half,_Phi0Max);
        if ( !is_physical )
        {
          if (_VarDtFlag)
          {
            dt_new = 0.8*dt;
            continue;
          }
          else
          {
            std::cout << "Phi becomes unphysical at this timestep ";
            std::cout << "(dt = " << dt << ")." << std::endl;
            exit(-1);
          }
        }

        // whole step for phi
        step_phi( phi_whole, mu_current, mu_current_linexp, mu_current_linimp, dt, t );
        is_physical = check_phi(phi_whole,_Phi0Max);
        if ( !is_physical )
        {
          if (_VarDtFlag)
          {
            dt_new = 0.8*dt;
            continue;
          }
          else
          {
            std::cout << "Phi becomes unphysical at this timestep ";
            std::cout << "(dt = " << dt << ")." << std::endl;
            exit(-1);
          }
        }

      }

      // difference between two stepping schemes
      phi_whole -= phi_half;

      // estimate truncation error, step doubling dt
      phi_err = 0.;
      for( UInt m=0; m<_NIcomp; ++m )
      {
        phi_err = std::max( phi_err, std::abs(phi_whole[m].max().real()) );
      }

      phi_iter++;

      //next 2 lines commented out after the introduction of outer_loop_const 
      //if(!_VarDtFlag) // If not adaptive time stepping, exit the dt refinement loop
      //  break; 
      
      dt_new = update_timestep( phi_err, dt);

    }

    if(_EnableFluctuations)
    { // {{{
      //Noise calculation is extracted out of step_phi for adaptive time-stepping
      //so that it is not repeated during the dt refinement inner loop.
      //
      //KTD: The mobility and semi-implicit operators are recomputed here
      //      even though they were computed in step_phi.
      //      Review order of operations and optimize this away so that this task is not repeated.
      
      // Calculate Mobility: same as within step_phi
      SmartFieldMat Mobility( _NIcomp, _NIcomp );
      _MobilityModel->set_mobility( _Phi, Mobility );
      Vec2dFieldType M_max = _MobilityModel->get_mobility_max( Mobility );
      
      // Calculate lhs: same as within step_phi 
      SmartFieldOpMat lhs( _NIcomp, _NIcomp );
      lhs.zero();
      lhs.setflag_inrealspace( false );
      _Phi.fft(); //calc_lhs requires _Phi in k-space
      calc_lhs ( _Phi, M_max, mu_current_linimp, dt, lhs );
      _Phi.ifft(); 
      
      // Calculate rhs (noise only): same as within step_phi, for constant time-stepping
      SmartFieldVec noise_rhs(_NIcomp);
      calc_noise_phi ( Mobility, _ReductionFactor, dt, noise_rhs );

      // -----------------------------
      // Solve for noise(n+1)
      // lhs noise(n+1) = rhs
      //    lhs = implicit ops
      //    rhs = noise(n)
      // -----------------------------
      SmartFieldVec noise (_NIcomp);
      _PhiOp->Solve_im_matrix( lhs, noise_rhs, noise, t );

      // Add the noise in real space
      noise.ifft();
      phi_half += noise;
    
    } // }}}

    // --- time step ---
    _Phi = phi_half;
    t = t + dt;
    Nt++;

    bool enableSK(false);
    if(enableSK)
    { // {{{
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
          //tmp *= _CurrGrid->GetCell()->getVol();
          //SK[k] += tmp;
          SK[k].xpby_inplace(tmp,V);
        }
      // Revert phi fields to (r) space
      for(int i=0; i<_NIcomp; i++)
        _Phi[i].fft_ktor();
    } // }}}

    // --- write diagnostics and data ---
    // step doubling diagnostics
    if ( Nt % _NStepDisp == 0 || t >= t_check || t == _TMax || Nt == _NtMax )
    {
      stepfile << Nt << "\t" << t << "\t" << dt << "\t";
      stepfile << phi_err << "\t" << phi_iter << "\t" <<  _PhiErrTol << std::endl;
      printf("%8lu%12.3e%12.3e%12.3e%12lu%12.3e", Nt, t, dt, phi_err, phi_iter, _PhiErrTol);
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

      // Write S(k) on the display steps
      if(enableSK)
      {
        std::stringstream fname;
        int k=0;
        for(int i=0; i<_NIcomp; i++)
          for(int j=i; j<_NIcomp; j++, k++)
          {
            SKhist.mapfromfield(SK[k].GetField(0), numSKsamples);
            fname.str("");
            fname<<"SK"<<i<<j<<".dat";
            SKhist.writehistogram(fname.str().c_str());
          }
      }
    }

    // write data to file
    if ( t - t_check >= -1e-14 || t == _TMax || Nt == _NtMax )
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

void TimeInt_ModelB :: step_phi( SmartFieldVec &phi,
                                 SmartFieldVec const &mu,
                                 SmartFieldVec const &mu_lin,
                                 SmartFieldOpMat const &mu_lin_np1,
                                 RealType dt,
                                 RealType t )
{ // {{{

  // Solves the following (1 step of Model B):
  //
  // ( I - dt Mij_max del^2 mu_lin_op_j  ) phi_n+1 = 
  //     phi_n + dt del . (Mij del . mu_n) - dt*Mij_max del^2 mu_lin_j
  //

  // Get Mobility
  SmartFieldMat Mobility( _NIcomp, _NIcomp );
  _MobilityModel->set_mobility( phi, Mobility );

  // Get maximum mobility
  Vec2dFieldType M_max = _MobilityModel->get_mobility_max( Mobility );

  // -----------------------------
  // Find rhs (all explicit terms)
  // rhs = phi_n + dt del . (Mij del . mu_n) - dt Mij_max del^2 mu_lin_j
  //             + dt Rxn_terms - dt Rxn_terms_linear
  //             + noise(already integrated)
  // -----------------------------

  SmartFieldVec rhs( _NIcomp );
  SmartFieldVec tmpvec( _NIcomp );
  SmartFieldVec tmpvec2( _NIcomp );

  // phi_n
  phi.fft();
  rhs = phi;

  // get full mu term
  tmpvec.setflag_inrealspace( false );
  _PhiOp->Del_A_Del_f_ex( mu, Mobility, tmpvec ); // tmp <- del.(Mij del. mu_n)
  rhs.xpby_inplace(tmpvec, dt);
  
  // get full reactions term
  if (_ReactionsModel->GetReactionsFlag() != 0)
  {
    _ReactionsModel->calc_reactions( phi, tmpvec );
    rhs.xpby_inplace(tmpvec, dt);
  }

  // add full noise term only in constant time-stepping
  // for adaptive time-stepping, noise is stepped outside of step_phi for performance
  if(_EnableFluctuations && !(_VarDtFlag))
  {
    calc_noise_phi ( Mobility, _ReductionFactor, dt, tmpvec );// tmp <- noise
    rhs.xpby_inplace( tmpvec, 1.0);
  }

  // subtract off the explicit linear term of mu, which will be treated implicitly
  _PhiOp->Del2_f_ex( mu_lin, tmpvec ); // tmpvec <- del^2 mu_lin_n
  tmpvec2.dot(M_max, tmpvec); // mu_lin_n <- M_max . tmpvec
  rhs.xpby_inplace(tmpvec2, -dt);
  
  // subtract off the linear reactions term, which will be treated implicitly
  if (_ReactionsModel->GetReactionsFlag() != 0)
  {
    _ReactionsModel->calc_reactions_lin_exp( phi, tmpvec );
    rhs.xpby_inplace( tmpvec, -dt);
  }

  // -----------------------------
  // Find lhs (all implicit terms)
  // lhs = I - dt Mij_max del^2 mu_lin_op_j
  // -----------------------------
  SmartFieldOpMat lhs( _NIcomp, _NIcomp );
  lhs.zero();
  lhs.setflag_inrealspace( false );
  
  calc_lhs ( phi, M_max, mu_lin_np1, dt, lhs );

  // -----------------------------
  // Solve for phi_n+1
  // lhs (phi_n+1) = rhs
  // -----------------------------
  _PhiOp->Solve_im_matrix( lhs, rhs, phi, t );
  phi.ifft();
  phi.zeroimag();

  if ( (_CapPhi0Flag == 1) && (phi[0].maxsigned().real() > _Phi0Max) )
  {
    capField(phi[0], _Phi0Max);
  }

  return;

} // }}}

void TimeInt_ModelB :: calc_lhs ( SmartFieldVec const &phi, Vec2dFieldType const &M_max,
                                  SmartFieldOpMat const &mu_lin_np1, RealType const &dt, 
				                          SmartFieldOpMat &lhs )
{ // {{{
  // Calculate the lhs (i.e., the A matrix) in  Ax=b of step_phi
  // Note: phi comes in k-space (used for calculation of rxn terms only)

  // -----------------------------------------------------
  // lhs = I - dt Mij_max del^2 mu_lin_op_j - dt Rxn_lin
  // -----------------------------------------------------

  // I
  _PhiOp->Eye_im( lhs );

  // subtract the implicit, linear part of mu
  SmartFieldOpMat tmpop( _NIcomp, _NIcomp );
  tmpop.zero();
  tmpop.setflag_inrealspace( false );
  _PhiOp->A_Del2_F_im( M_max, mu_lin_np1, tmpop ); // tmpop <- M_max del^2 mu_lin_np1
  tmpop *= (-dt);
  lhs += tmpop;
  
  // subtract the implicit, linear part of the reaction terms
  if (_ReactionsModel->GetReactionsFlag() != 0)
  {
    _ReactionsModel->calc_reactions_lin_imp( phi, tmpop );
    tmpop *= (-dt);
    lhs += tmpop;
  }

} // }}}


RealType TimeInt_ModelB :: update_timestep( RealType phi_err, RealType dt )
{ // {{{

  static bool prev_step_accept = true;
  RealType theta = 0.95; // safety factor so we don't reject steps if possible
  RealType dt_new = 0.;   // new step size

  // only step bigger if (1) big enough and (2) didn't reject prev step; always step smaller
  if ( (prev_step_accept && ((dt_new - dt)/dt > 0.01)) || dt_new < dt )
  {
    // dt_step_dbl = dt * sqrt{ _PhiErrTol/ max[abs(phi_err)] }
    dt_new = dt*sqrt( theta*_PhiErrTol/phi_err );
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

