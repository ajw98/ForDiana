#define MAIN

#include "version.h"
#include "global.h"
#include "pllhandler.h"
#include "Grid.h"
#include "FieldStack.h"
#include "SmartField.h"
#include "Operators_Base.h"
#include "Operators_Factory.h"
#include "Model_Energy_Base.h"
#include "Model_Factory.h"
#include "TimeInt_Base.h"
#include "TimeInt_Factory.h"

#include "cRandom.h"

void read_parallel(int &MPInnodes, int &OMPnthreads, int &CUDAdeviceID, 
                   int &CUDAnthreads, int &maxnumelements);
void read_keys( int & Ncomp, int & NIcomp, int & GridKey, int & EnergyKey, 
                int & MobilityKey, int & ViscKey, int & TimeIntKey,
                int & ReactionsKey );

int main(int argc, char **argv)
{

  // --- Initialize parallelism (type determined at run time) ---

  pllhandler pll(argc, argv);

  int MPInnodes, OMPnthreads, CUDAdeviceID, CUDAnthreads, maxnumelements;
  read_parallel( MPInnodes, OMPnthreads, CUDAdeviceID, 
                 CUDAnthreads, maxnumelements );

  pll.setup(OMPnthreads,CUDAnthreads,maxnumelements,CUDAdeviceID);

#ifdef __MPI__
  // if not root node, send cout/printf to null stream
  std::cout.flush();
  MPI_Barrier(MPI_COMM_WORLD);
  if ( !pll.isrootnode() ) 
  {
    std::cout.rdbuf(0); // redirect cout
    freopen("", "w", stdout); // redirect printf
  }
#endif

  // Now that the input parameters are known, and the pll strategy is fully initialized
  // set up the random number generator...random number seed differently on each MPI node
  // NOTE: it is crucial that we don't do this until after pll.setup!
  long rndseed = 0;
  // TODO: allow user input of random seed here. If the default of 0 is maintained then
  // it is set automatically below.
  if(rndseed == 0)
  {
    // If available on this system, we seed the random number generator using
    // /dev/urandom.
    // If not available, we seed from the system clock. The latter has the disadvantage of
    // potentially producing correlated sequences of random numbers between different runs started at the
    // same time.
    std::ifstream urnd;
    urnd.open("/dev/urandom",std::ios::in | std::ios::binary);
    if(urnd.fail())
      rndseed = -long(time(0));
    else
      urnd.read(reinterpret_cast<char*>(&rndseed),sizeof(rndseed)/sizeof(char));
    // Restriction on seed: it should be negative and < 2147483565
    // NOTE: these restrictions should really apply to inputted random seeds, but to avoid
    // breaking reproducibility of previous runs this is not enforced. Rather, test for
    // it below and warn.
    rndseed %= 2147483565;
    rndseed = -std::abs(rndseed); // Ensure negative
  }
  // Make the seed different on every MPI node
  rndseed *= (pll.myrank()+1); // Again make the seed different on every node
  std::cout << " Random number seed (node0)        = " << rndseed << std::endl;
  if(rndseed >= 0 || rndseed  < -2147483565)
    std::cout << " WARNING: random seeds should be in the range [-2147483565,0)" << std::endl;
  cRandom::Instance().initGenerator(rndseed,false);

  // --- Time Stamp ---
  struct timeval begtime, endtime;
  double totaltime;

  gettimeofday(&begtime, NULL);
  printf("\nBeginning Simulation: %s\n", ctime(&begtime.tv_sec) );
  printf(" * Code Branch: %s\n", GIT_BRANCH);
  printf(" * Code Date: %s\n", GIT_DATE);
  printf(" * Code Hash: %s\n", GIT_SHA);

  std::cout << std::endl;
  std::cout << " * Parallelization initiated\n" << std::scientific;
  std::cout << "   - Parallel Parameters:\n";
  std::cout << "      MPI Nodes       = " << MPInnodes << "\n";
  std::cout << "      OMP threads     = " << OMPnthreads << "\n";
  std::cout << "      Cuda Dev ID     = " << CUDAdeviceID << "\n";
  std::cout << "      Cuda threads    = " << CUDAnthreads << "\n";
  std::cout << "      Max GPU elem    = " << maxnumelements << "\n";
                
  // --- Get # components and keys for initializing classes ---
  int Ncomp, NIcomp, GridKey, EnergyKey;
  int MobilityKey, ViscKey, TimeIntKey;
  int ReactionsKey;

  read_keys( Ncomp, NIcomp, GridKey, EnergyKey, 
             MobilityKey, ViscKey, TimeIntKey,
             ReactionsKey );
      
  // --- initialize the grid/layout ---
                      
  bool IsPeriodic = (GridKey == 0 ? true : false);
  Grid* CurrGrid = new Grid( pll, IsPeriodic );
                      
  // --- initialize the SmartFields ---
  // Need to set up SmartFields *after* we initialize the cell 
  // and the layout, but *before* we initialize the gradients
                      
  init_FieldStack_statics( CurrGrid->GetLayout() );
  init_SmartField_statics( CurrGrid->FDSize() );

  // --- initialize the time integration class and run the outer loop ---

  TimeInt_Factory *TIFactory = new TimeInt_Factory();
  TimeInt_Base *TIObj = 
    TIFactory->make_timeint(TimeIntKey, Ncomp, NIcomp, CurrGrid );

  // initialize the operators and boundary conditions
  TIObj->setup_operators();

  // initialize the models
  TIObj->setup_models( EnergyKey, MobilityKey, ViscKey, ReactionsKey );

  // set integration parameters (read from file)
  TIObj->set_params();

  // read in initial conditions from file
  TIObj->read_initial_state();

  // ******************************
  //   main time integration loop
  //  
         TIObj->outer_loop(); 
  //  
  // ******************************

  //write out final conditions to file
  TIObj->write_final_state();

  // clean up models
  TIObj->cleanup_models();

  // clean up operators
  TIObj->cleanup_operators();

  delete TIObj;
  delete TIFactory;

  // --- clean up fields and grids ---

  // Need to delete the SmartFields before we delete
  // the Grid but after we delete the Operators
  delete_SmartField_statics();
  delete_FieldStack_statics();

  delete CurrGrid;

  // --- Time Stamp ---
  gettimeofday(&endtime, NULL); // time of day
  totaltime = (endtime.tv_sec - begtime.tv_sec) + 
              (endtime.tv_usec - begtime.tv_usec)/1000000.0;
  printf("\nEnding Simulation: %s", ctime(&endtime.tv_sec) );
  printf(" - Elapsed Wall Time: %10.3e (sec.)\n\n", totaltime );

  return 0;
}

void read_parallel(int & MPInnodes, int & OMPnthreads, int & CUDAdeviceID, 
                   int & CUDAnthreads, int & maxnumelements)
{ // {{{

  std::string tmp_str;
  std::ifstream fid;

  if (jsonFile("params_Parallel.in"))
  {
    fid.open("params_Parallel.in");
    if ( fid.is_open() )
    {
      rapidjson::IStreamWrapper isw(fid);
      rapidjson::Document d;
      d.ParseStream(isw);

      MPInnodes = d["MPInnodes"].GetInt();

      OMPnthreads = d["OMPnthreads"].GetInt();

      CUDAdeviceID = d["CUDAdeviceID"].GetInt();

      CUDAnthreads = d["CUDAnthreads"].GetInt();

      maxnumelements = d["maxnumelements"].GetInt();

    }
    else
    {
      std::cout << "(Error): No Parallel Parameters\n";
      exit(1);
    }
    fid.close();
  }
  else
  {
    
    fid.open("params_Parallel.in");
    if ( fid.is_open() )
    {
      std::getline(fid, tmp_str, '#');
      MPInnodes = atoi(tmp_str.c_str());
      std::getline(fid, tmp_str, '\n');

      std::getline(fid, tmp_str, '#');
      OMPnthreads = atoi(tmp_str.c_str());
      std::getline(fid, tmp_str, '\n');

      std::getline(fid, tmp_str, '#');
      CUDAdeviceID = atoi(tmp_str.c_str());
      std::getline(fid, tmp_str, '\n');

      std::getline(fid, tmp_str, '#');
      CUDAnthreads = atoi(tmp_str.c_str());
      std::getline(fid, tmp_str, '\n');

      std::getline(fid, tmp_str, '#');
      maxnumelements = atoi(tmp_str.c_str());
      std::getline(fid, tmp_str, '\n');
    }
    fid.close();
  }
  std::cout << "Preparing for the simulation\n"; 
  std::cout << "\n" << std::scientific;
  std::cout << " * Parameters from params_Parallel.in:\n" << std::scientific;
  std::cout << "      MPInnodes          = " << MPInnodes << "\n";
  std::cout << "      OMPnthreads        = " << OMPnthreads << "\n";
  std::cout << "      CUDAdeviceID       = " << CUDAdeviceID << "\n";
  std::cout << "      CUDAnthreads       = " << CUDAnthreads << "\n";
  std::cout << "      maxnumelements     = " << maxnumelements << "\n";
  std::cout << std::endl;

} // }}}

void read_keys( int & Ncomp, int & NIcomp, int & GridKey, 
                int & EnergyKey, int & MobilityKey, 
                int & ViscKey, int & TimeIntKey, int & ReactionsKey )
{ // {{{

  std::string tmp_str;
  std::ifstream f1;
  
  if (jsonFile("params_Keys.in"))
  {
    f1.open("params_Keys.in");
    if ( f1.is_open() )
    {
      rapidjson::IStreamWrapper isw(f1);
      rapidjson::Document d;
      d.ParseStream(isw);

      Ncomp = d["Ncomp"].GetInt();
      if (Ncomp < 2)
      {
        std::cout << "*** (Run-time Error) in main.cpp::read_keys() ***\n";
        std::cout << "*** Number of components (" << Ncomp;
        std::cout << ") must be >= 2. ***\n";
      }

      NIcomp = d["NIcomp"].GetInt();

      GridKey = d["GridFlag"].GetInt();

      EnergyKey = d["EnergyModelFlag"].GetInt();

      MobilityKey = d["MobilityModelFlag"].GetInt();

      ViscKey = d["ViscModelFlag"].GetInt();

      TimeIntKey = d["TimeIntFlag"].GetInt();

      ReactionsKey = d["ReactionsModelFlag"].GetInt();

    }
    else
    { 
      std::cout << "Error opening params_Keys.in" << std::endl;
      exit(1);
    }
    f1.close();
  }
  else
  {
    f1.open("params_Keys.in");
    if ( f1.is_open() )
    {
      std::getline(f1, tmp_str, '#');
      Ncomp = atoi(tmp_str.c_str());
      std::getline(f1, tmp_str, '\n');
      
      if (Ncomp < 2)
      {
        std::cout << "*** (Run-time Error) in main.cpp::read_keys() ***\n";
        std::cout << "*** Number of components (" << Ncomp;
        std::cout << ") must be >= 2. ***\n";
      }

      std::getline(f1, tmp_str, '#');
      NIcomp = atoi(tmp_str.c_str());
      std::getline(f1, tmp_str, '\n');

      std::getline(f1, tmp_str, '#');
      GridKey = atoi(tmp_str.c_str());
      std::getline(f1, tmp_str, '\n');

      std::getline(f1, tmp_str, '#');
      EnergyKey = atoi(tmp_str.c_str());
      std::getline(f1, tmp_str, '\n');

      std::getline(f1, tmp_str, '#');
      MobilityKey = atoi(tmp_str.c_str());
      std::getline(f1, tmp_str, '\n');

      std::getline(f1, tmp_str, '#');
      ViscKey = atoi(tmp_str.c_str());
      std::getline(f1, tmp_str, '\n');

      std::getline(f1, tmp_str, '#');
      TimeIntKey = atoi(tmp_str.c_str());
      std::getline(f1, tmp_str, '\n');

      std::getline(f1, tmp_str, '#');
      ReactionsKey = atoi(tmp_str.c_str());
      std::getline(f1, tmp_str, '\n');

    }
    else
    { 
      std::cout << "Error opening params_Keys.in" << std::endl;
      exit(1);
    }
    f1.close();
  }

  std::cout << "\n" << std::scientific;
  std::cout << " * Reading # of components and parameter keys:\n" << std::scientific;
  std::cout << "   - Components" << "\n";
  std::cout << "      Ncomp           = " << Ncomp << "\n";
  std::cout << "      NIcomp          = " << NIcomp << "\n";
  std::cout << "   - Keys" << "\n";
  std::cout << "      GridKey         = " << GridKey << "\n";
  std::cout << "      EnergyKey       = " << EnergyKey << "\n";
  std::cout << "      MobilityKey     = " << MobilityKey << "\n";
  std::cout << "      ViscKey         = " << ViscKey << "\n";
  std::cout << "      TimeIntKey      = " << TimeIntKey << "\n";
  std::cout << "      ReactionsKey    = " << ReactionsKey << "\n";

} // }}}

