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

#include "TimeInt_ModelB_DelMu.h"

// ----------------- Constructor/Destructor -----------------

TimeInt_ModelB_DelMu :: TimeInt_ModelB_DelMu( int Ncomp, int NIcomp, Grid* CurrGrid ):
TimeInt_Base( Ncomp, NIcomp, CurrGrid ),
_EnableFluctuations(false),
_Phi( NIcomp )
// {{{
{

  std::cout << std::endl;
  std::cout << " * Time Integration derived class initialization\n";
  std::cout << "   - Model B (diffusion only) using Del_Mu terms" << std:: endl;

} // }}}

TimeInt_ModelB_DelMu :: ~TimeInt_ModelB_DelMu()
{ // {{{
} // }}}

// ----------------- Input/Output -----------------

void TimeInt_ModelB_DelMu :: set_params()
{//{{{

  std::string tmp_str;

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

void TimeInt_ModelB_DelMu::read_initial_state()
{ // {{{

  // read in initial conditions
  char filename[40];
  for (int i=0; i<_NIcomp; i++)
  {
    sprintf( filename, "phi%d.in", i+1 );
    _Phi[i].readfield(filename);
  }

} // }}}

void TimeInt_ModelB_DelMu::write_final_state()
{ // {{{

  // write out final state
  char filename[40];
  for (int i=0; i<_NIcomp; i++)
  {
    sprintf( filename, "phi%d.out", i+1 );
    _Phi[i].writerealfield(filename);
  }

} // }}}

void TimeInt_ModelB_DelMu :: write_curr_state( int n_disp )
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

void TimeInt_ModelB_DelMu :: outer_loop()
{ // {{{ 

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

  // Keep a copy of mu_n and mu_n+1/2 to avoid repeat evaluations in stepper routines
  SmartFieldMat del_mu_current(_NIcomp, _Dim);
  SmartFieldVec mu_current_linexp(_NIcomp);
  SmartFieldOpMat mu_current_linimp(_NIcomp,_NIcomp);
  SmartFieldMat del_mu_half(_NIcomp, _Dim);
  SmartFieldVec mu_half_linexp(_NIcomp);
  SmartFieldOpMat mu_half_linimp(_NIcomp,_NIcomp);
  del_mu_current.setflag_inrealspace(false);
  del_mu_half.setflag_inrealspace(false);

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

    // Compute and store current mu before any adaptive time stepper looping
    _Phi.fft();
    _EnergyModel->calc_del_mu(_Phi, del_mu_current);
    _EnergyModel->calc_mu_lin_imp(_Phi, mu_current_linimp);
    _EnergyModel->calc_mu_lin_exp(_Phi, mu_current_linexp);
    //mu_current_linimp.OpDot(_Phi, mu_current_linexp);
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
        step_phi( phi_half, del_mu_current, mu_current_linexp, mu_current_linimp, dt/2., t );
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
        _EnergyModel->calc_del_mu(phi_half, del_mu_half);
        _EnergyModel->calc_mu_lin_imp(phi_half, mu_half_linimp);
        _EnergyModel->calc_mu_lin_exp(phi_half, mu_half_linexp);
        //mu_half_linimp.OpDot(phi_half, mu_half_linexp);
        phi_half.ifft();

        // second half step for phi
        step_phi( phi_half, del_mu_half, mu_half_linexp, mu_half_linimp, dt/2., t ); 
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
        step_phi( phi_whole, del_mu_current, mu_current_linexp, mu_current_linimp, dt, t ); 
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
      if(!_VarDtFlag) // If not adaptive time stepping, exit the dt refinement loop
        break;
      dt_new = update_timestep( phi_err, dt);

    }

    // Generate noise for fluctuating simulations
    // TODO: this method only works for PS algorithms. Either generalize it, or test for FD and quit.
    if(_EnableFluctuations)
    { // {{{
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
      L.setflag_inrealspace(Mobility.getflag_inrealspace());
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

void TimeInt_ModelB_DelMu :: step_phi( SmartFieldVec &phi,
                                 SmartFieldMat const &del_mu,
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
  // -----------------------------

  SmartFieldVec rhs( _NIcomp );
  SmartFieldVec tmpvec( _NIcomp );
  SmartFieldVec tmpvec2( _NIcomp );

  // phi_n
  phi.fft();
  rhs = phi;

  // get full mu term
  tmpvec.setflag_inrealspace( false );
  _PhiOp->Del_A_Del_f_ex( del_mu, Mobility, tmpvec ); // tmp <- del.(Mij del. mu_n)
  rhs.xpby_inplace(tmpvec, dt);
  if (_ReactionsModel->GetReactionsFlag() != 0)
  {
    _ReactionsModel->calc_reactions( phi, tmpvec );
    rhs.xpby_inplace(tmpvec, dt);
  }

  // Subtract off the explicit linear term of mu, which will be treated implicitly
  _PhiOp->Del2_f_ex( mu_lin, tmpvec ); // tmpvec <- del^2 mu_lin_n
  tmpvec2.dot(M_max, tmpvec); // mu_lin_n <- M_max . tmpvec
  rhs.xpby_inplace(tmpvec2, -dt);
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

  // I
  _PhiOp->Eye_im( lhs );

  // implicit, linear part of mu
  SmartFieldOpMat tmpop( _NIcomp, _NIcomp );
  tmpop.zero();
  tmpop.setflag_inrealspace( false );
  _PhiOp->A_Del2_F_im( M_max, mu_lin_np1, tmpop ); // tmpop <- M_max del^2 mu_lin_np1
  tmpop *= (-dt);
  lhs += tmpop;
  if (_ReactionsModel->GetReactionsFlag() != 0)
  {
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

  if ( (_CapPhi0Flag == 1) && (phi[0].maxsigned().real() > _Phi0Max) )
  {
    capField(phi[0], _Phi0Max);
  }

  return;

} // }}}

RealType TimeInt_ModelB_DelMu :: update_timestep( RealType phi_err, RealType dt )
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

