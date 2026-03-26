/***************************************************************
*
* Base class for the time integration schemes
* 
* DRT -- Mon 16 May 2016
*
****************************************************************/

#include "TimeInt_Base.h"

// ----------------- Constructor/Destructor -----------------

TimeInt_Base :: TimeInt_Base( int Ncomp, int NIcomp, Grid* CurrGrid) :
_Ncomp(Ncomp),
_NIcomp(NIcomp),
_CurrGrid(CurrGrid),
_EnergyModel( NULL ),
_MobilityModel( NULL ),
_ViscModel( NULL ),
_ReactionsModel( NULL ),
_PhiOp( NULL ),
_PhiBCs( NULL )
{ // {{{

  _Dim = _CurrGrid->Dim();
  _OpFactory = new Operators_Factory();
  _MFactory = new Model_Factory();
  _BCFactory = new BCs_Factory();

} // }}}

TimeInt_Base :: ~TimeInt_Base()
{ // {{{

  delete _MFactory;
  delete _OpFactory;
  delete _BCFactory;

} // }}}

// -------------------------------- Thermal fluctuations -----------------------------------------------

void TimeInt_Base :: calc_noise_phi ( SmartFieldMat const &M, RealType const rFac,
                                      RealType const dt, SmartFieldVec &noise)
{ // {{{
  // Calculates noise fields using the fluctuation-disspiation theorem for Model B with timestep dt. 
  // Inputs:  M      = mobility matrix
  //          rFac   = reduction factor (scales noise amplitude std. deviation artificially)
  //          dt     = current time step size
  // Outputs: noise  =  noise fields, same index as the mobility matrix
  // 
  // We output the calculated noise in k-space as subsequent steps (solving Ax=b) often take place in k-space.

  // KTD: this method only works for PS algorithms.
  // Either generalize it, or test for FD and quit.
  if (_CurrGrid->GridFlag() != 0) 
  {
    std::cout << "Fluctuation-dissipation noise generation is not yet supported for non-pseudo-spectral grids." << std::endl;
    std::cout << "Terminating execution." << std::endl;
    exit(1);
  }

  // Cholesky decompose the mobility matrix
  SmartFieldMat L(_NIcomp, _NIcomp);
  L.setflag_inrealspace(M.getflag_inrealspace());
  // Begin by copying the upper triangle of M to L
  for(UInt i=0; i<_NIcomp; i++)
    for(UInt j=i; j<_NIcomp; j++)
      L(i,j) = M(i,j);
  // Cholesky decomposition: the lower triangle of L will be replaced
  L.cholesky();

  // NOISE - generate one noise field per non-implicit phi field
  // This must be done for every field before adding noise to any field,
  // because any off-diagonal contributions in
  // the mobility matrix will mix the noise fields together
  SmartFieldVec noise_tmp(_NIcomp);
  SmartField    noisecoef;
  double dV = _CurrGrid->GetCell()->getVol() / _CurrGrid->NPW();
  // Fill a field with sqrt(2*dt/dV)*|k|
  noisecoef.setflag_inrealspace(false);
  noisecoef = -2.*dt/dV/sqrt(_EnergyModel->mu_scale())/rFac/rFac;//divided 2x due to sqrt later.
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
  
  // Note also that the local approximation we used to treat the mobility
  // is not conserving.  Enforce conservation manually
  noise_tmp.fft();
  for(UInt i=0; i<_NIcomp; i++)
    noise_tmp[i].setaverage(0.);
  
  noise = noise_tmp;

} // }}}

void TimeInt_Base :: calc_noise_vel ( SmartField const &eta, RealType const rFac,
                                      SmartFieldVec &noise)
{ // {{{
  // Calculates noise fields using the fluctuation-disspiation theorem 
  // for Model H momenta with timestep dt. 
  // Inputs:  eta    = viscosity field
  //          rFac   = reduction factor (scales noise amplitude std. deviation artificially)
  // Outputs: noise  =  noise fields, noise(x,y,z)
  // 
  // We output the calculated noise in k-space as subsequent steps (solving Ax=b) often take place in k-space.

  // KTD: this method only works for PS algorithms.
  // Either generalize it, or test for FD and quit.
  if (_CurrGrid->GridFlag() != 0) 
  {
    std::cout << "Fluctuation-dissipation noise generation is not yet supported for non-pseudo-spectral grids." << std::endl;
    std::cout << "Terminating execution." << std::endl;
    exit(1);
  }

  // NOISE - generate one noise field per velocity field, vx, vy, vz
  SmartFieldVec noise_tmp(_Dim);
  SmartField    noisecoef;
  double dV = _CurrGrid->GetCell()->getVol() / _CurrGrid->NPW();
  // Fill a field with sqrt(2*dt/dV)*|k|
  noisecoef.setflag_inrealspace(false);
  noisecoef = -2./dV*sqrt(_EnergyModel->mu_scale())/rFac/rFac;//divided 2x due to sqrt later (variance)
  _PhiOp->Del2_f_ex_inplace(noisecoef); // * (-k^2)
  noisecoef.sqrt(noisecoef);//(std dev)
  
  // Generate _Dim noise fields with variance 2*k^2*dt/dV
  for(UInt i=0; i<_Dim; i++)
  {
    noise[i].setflag_inrealspace(true);
    noise[i].fillnoise(3+_CurrGrid->Dim());
    noise[i].fft_rtok();
    noise[i] *= noisecoef;
    noise[i].fft_ktor();
  }
 
  // Multiply the square root of viscosity in real-space
  SmartField eta_tmp;
  eta_tmp = eta;
  eta_tmp.sqrt(eta_tmp);
  for(UInt i=0; i<_Dim; i++)
  {
    noise[i] *= eta_tmp;
  }
  
  // Note also that the local approximation we used to treat the viscosity
  // is not conserving.  Enforce conservation manually
  noise.fft();
  for(UInt i=0; i<_NIcomp; i++)
    noise[i].setaverage(0.);

} // }}}

// ----------------- Field checks and manipulations -----------------
void TimeInt_Base :: capField ( SmartField &F, const RealType k )
{ // {{{
  // Caps the field to a maximum value of k. 
  F -= k;
  F.zeropos();
  F += k;
} // }}}{

void TimeInt_Base :: floorField ( SmartField &F, const RealType k )
{ // {{{
  // Floors the field to a minimum value of k. 
  F -= k;
  F.zeroneg();
  F += k;
} // }}}{

void TimeInt_Base :: floorField ( SmartFieldVec &F, const RealType k )
{ // {{{
  // Floors all fields to a minimum value of k. 
  for( UInt i=0; i<_NIcomp; ++i )
  {
    floorField( F[i], k );
  }
} // }}}{

bool TimeInt_Base :: check_phi( const SmartFieldVec &phi )
{ // {{{

  SmartField phi_tmp;
  bool phiphysical(true);

  for( UInt i=0; i<_Ncomp; ++i )
  {
    if (i < _NIcomp)
    { 
      phi_tmp = phi[i]; 
    }
    else
    {
      phi_tmp = 1.;
      for ( UInt j=0; j<_NIcomp; ++j ){ phi_tmp -= phi[j]; }
    }

    if ( phi_tmp.maxsigned().real() > 1. || phi_tmp.minsigned().real() < 0.)
    {
      std::cout << "*** Warning from TimeInt_Base::check_phi() ***\n";
      std::cout << "phi: " << i << " is unphysical:";
      std::cout << " max(phi) = " << phi_tmp.maxsigned().real();
      std::cout << ", min(phi) = " << phi_tmp.minsigned().real();
      std::cout << std::endl;
      return false; // unphysical value of phi
    }
  }
  return true; // phi is physical
} // }}}

bool TimeInt_Base :: check_phi( const SmartFieldVec &phi, const RealType &glassMax )
{ // {{{

  SmartField phi_tmp;
  bool phiphysical(true);

  for( UInt i=0; i<_Ncomp; ++i )
  {
    if (i < _NIcomp)
    { 
      phi_tmp = phi[i]; 
    }
    else
    {
      phi_tmp = 1.;
      for ( UInt j=0; j<_NIcomp; ++j ){ phi_tmp -= phi[j]; }
    }
    /*if (i==0)//check glassMax for polymer only --- make general later
    {
      if (phi_tmp.maxsigned().real() > glassMax)
      {
        std::cout << "*** Warning from TimeInt_Base::check_phi() ***\n";
        std::cout << "Calculated max(phi_polymer) = " << phi_tmp.maxsigned().real();
        std::cout << " is greater than glass transition check = " << glassMax << std::endl;
        return false; //unphysical value of phi 
      }  
    }*/
    if ( phi_tmp.maxsigned().real() > 1. || phi_tmp.minsigned().real() < 0.)
    {
      std::cout << "*** Warning from TimeInt_Base::check_phi() ***\n";
      std::cout << "phi: " << i << " is unphysical:";
      std::cout << " max(phi) = " << phi_tmp.maxsigned().real();
      std::cout << ", min(phi) = " << phi_tmp.minsigned().real();
      std::cout << std::endl;
      return false; // unphysical value of phi
    }
  }
  return true; // phi is physical
} // }}}

// ----------------- Setup/Cleanup -----------------

// 
void TimeInt_Base :: setup_operators( )
{ // {{{

  // set up the BCs and operators for phi
  std::string filename( "params_BCs_phi.in" );

  // Read in BC key from file
  if (jsonFile(filename.c_str()))
  {
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
      std::cout << "*** Error in TimeInt_Base :: setup_operators() ***\n";
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
      std::cout << "*** Error in TimeInt_Base :: setup_operators() ***\n";
      std::cout << "*** Cannot open " << filename << " ***\n";
      exit(1);
    }
  }

  // call factory to create new BC object based on the kind of BCs
  _PhiBCs=_BCFactory->make_BCs( _BCFlag, 0, filename, _NIcomp, _CurrGrid );

  // call factory to make new operators object for phi
  _PhiOp = _OpFactory->make_operators( _CurrGrid, _PhiBCs );

} // }}}

void TimeInt_Base :: cleanup_operators()
{ // {{{

  _OpFactory->recycle_operators( _PhiOp );
  delete _PhiBCs;

} // }}}

void TimeInt_Base :: setup_models( int EnergyKey, int MobilityKey, int ViscKey, int ReactionsKey )
{ // {{{

  std::cout << "\n";
  std::cout << " * Initializing Models\n";

  // Energy Model
  _EnergyModel = _MFactory->make_energy_model( EnergyKey, *_PhiOp, _Ncomp, _NIcomp );
  _EnergyModel->set_params();

  // Viscosity Model
  // (set before Mobility Model, because Mobility Model needs access to eta)
  _ViscModel = _MFactory->make_visc_model( ViscKey, _Ncomp, _NIcomp );
  _ViscModel->set_params();

  // Mobility Model
  _MobilityModel = _MFactory->make_mobility_model( MobilityKey, _Ncomp, _NIcomp, _ViscModel );
  _MobilityModel->set_params();

  // Reactions Model
  _ReactionsModel = _MFactory->make_reactions_model( ReactionsKey, _Ncomp, _NIcomp );
  _ReactionsModel->set_params();

} // }}}

void TimeInt_Base :: cleanup_models()
{ // {{{

  delete _ViscModel;
  delete _MobilityModel;
  delete _EnergyModel;

} // }}}

