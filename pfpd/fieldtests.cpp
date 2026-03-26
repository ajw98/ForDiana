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

void read_parallel(int &MPInnodes, int &OMPnthreads, int &CUDAdeviceID, 
                   int &CUDAnthreads, int &maxnumelements);
void read_keys( int & Ncomp, int & NIcomp, int & GridKey, int & EnergyKey, 
                int & MobilityKey, int & ViscKey, int & TimeIntKey,
                int & ICKey, int & VelICKey, int & BCKey );

void unit_tests();

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
  int ICKey, VelICKey, BCKey;
  read_keys( Ncomp, NIcomp, GridKey, EnergyKey, 
             MobilityKey, ViscKey, TimeIntKey,
             ICKey, VelICKey, BCKey);
      
  // --- initialize the grid/layout ---
                      
  bool IsPeriodic = (BCKey == 0 ? true : false);
  Grid* CurrGrid = new Grid( pll, IsPeriodic );

  // --- initialize the SmartFields ---
  // Need to set up SmartFields *after* we initialize the cell 
  // and the layout, but *before* we initialize the gradients
                      
  init_FieldStack_statics( CurrGrid->GetLayout() );
  init_SmartField_statics( CurrGrid->FDSize() );

  // --- initialize the time integration class and run the outer loop ---

  unit_tests();

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

void unit_tests()
{ // {{{

  std::cout << "\n *** SmartField Tests ***\n";

  // the very first field takes a lot of time (allocating the FFT plan?)
  {
    SmartField tmp;
  }

  struct timeval begtime, midtime, endtime;
  double firstalloc_time, otheralloc_time, totaltime;
  gettimeofday(&begtime, NULL);
  {
    SmartFieldMat tmp0(10, 10);
//    SmartFieldVec tmp0(100);

//    SmartField tmp0;
//    SmartField tmp1;
//    SmartField tmp2;
//    SmartField tmp3;
//    SmartField tmp4;
//    SmartField tmp5;
//    SmartField tmp6;
//    SmartField tmp7;
//    SmartField tmp8;
//    SmartField tmp9;
//    SmartField tmp10;
//    SmartField tmp11;
//    SmartField tmp12;
//    SmartField tmp13;
//    SmartField tmp14;
//    SmartField tmp15;
//    SmartField tmp16;
//    SmartField tmp17;
//    SmartField tmp18;
//    SmartField tmp19;
//    SmartField tmp20;
//    SmartField tmp21;
//    SmartField tmp22;
//    SmartField tmp23;
//    SmartField tmp24;
//    SmartField tmp25;
//    SmartField tmp26;
//    SmartField tmp27;
//    SmartField tmp28;
//    SmartField tmp29;
//    SmartField tmp30;
//    SmartField tmp31;
//    SmartField tmp32;
//    SmartField tmp33;
//    SmartField tmp34;
//    SmartField tmp35;
//    SmartField tmp36;
//    SmartField tmp37;
//    SmartField tmp38;
//    SmartField tmp39;
//    SmartField tmp40;
//    SmartField tmp41;
//    SmartField tmp42;
//    SmartField tmp43;
//    SmartField tmp44;
//    SmartField tmp45;
//    SmartField tmp46;
//    SmartField tmp47;
//    SmartField tmp48;
//    SmartField tmp49;
//    SmartField tmp50;
//    SmartField tmp51;
//    SmartField tmp52;
//    SmartField tmp53;
//    SmartField tmp54;
//    SmartField tmp55;
//    SmartField tmp56;
//    SmartField tmp57;
//    SmartField tmp58;
//    SmartField tmp59;
//    SmartField tmp60;
//    SmartField tmp61;
//    SmartField tmp62;
//    SmartField tmp63;
//    SmartField tmp64;
//    SmartField tmp65;
//    SmartField tmp66;
//    SmartField tmp67;
//    SmartField tmp68;
//    SmartField tmp69;
//    SmartField tmp70;
//    SmartField tmp71;
//    SmartField tmp72;
//    SmartField tmp73;
//    SmartField tmp74;
//    SmartField tmp75;
//    SmartField tmp76;
//    SmartField tmp77;
//    SmartField tmp78;
//    SmartField tmp79;
//    SmartField tmp80;
//    SmartField tmp81;
//    SmartField tmp82;
//    SmartField tmp83;
//    SmartField tmp84;
//    SmartField tmp85;
//    SmartField tmp86;
//    SmartField tmp87;
//    SmartField tmp88;
//    SmartField tmp89;
//    SmartField tmp90;
//    SmartField tmp91;
//    SmartField tmp92;
//    SmartField tmp93;
//    SmartField tmp94;
//    SmartField tmp95;
//    SmartField tmp96;
//    SmartField tmp97;
//    SmartField tmp98;
//    SmartField tmp99;
  }

  gettimeofday(&midtime, NULL); // time of day

  for (int i=0; i<1000; i++)
  { // {{{

    SmartFieldMat tmp0(10, 10);
//    SmartFieldVec tmp0(100);

//    SmartField tmp0;
//    SmartField tmp1;
//    SmartField tmp2;
//    SmartField tmp3;
//    SmartField tmp4;
//    SmartField tmp5;
//    SmartField tmp6;
//    SmartField tmp7;
//    SmartField tmp8;
//    SmartField tmp9;
//    SmartField tmp10;
//    SmartField tmp11;
//    SmartField tmp12;
//    SmartField tmp13;
//    SmartField tmp14;
//    SmartField tmp15;
//    SmartField tmp16;
//    SmartField tmp17;
//    SmartField tmp18;
//    SmartField tmp19;
//    SmartField tmp20;
//    SmartField tmp21;
//    SmartField tmp22;
//    SmartField tmp23;
//    SmartField tmp24;
//    SmartField tmp25;
//    SmartField tmp26;
//    SmartField tmp27;
//    SmartField tmp28;
//    SmartField tmp29;
//    SmartField tmp30;
//    SmartField tmp31;
//    SmartField tmp32;
//    SmartField tmp33;
//    SmartField tmp34;
//    SmartField tmp35;
//    SmartField tmp36;
//    SmartField tmp37;
//    SmartField tmp38;
//    SmartField tmp39;
//    SmartField tmp40;
//    SmartField tmp41;
//    SmartField tmp42;
//    SmartField tmp43;
//    SmartField tmp44;
//    SmartField tmp45;
//    SmartField tmp46;
//    SmartField tmp47;
//    SmartField tmp48;
//    SmartField tmp49;
//    SmartField tmp50;
//    SmartField tmp51;
//    SmartField tmp52;
//    SmartField tmp53;
//    SmartField tmp54;
//    SmartField tmp55;
//    SmartField tmp56;
//    SmartField tmp57;
//    SmartField tmp58;
//    SmartField tmp59;
//    SmartField tmp60;
//    SmartField tmp61;
//    SmartField tmp62;
//    SmartField tmp63;
//    SmartField tmp64;
//    SmartField tmp65;
//    SmartField tmp66;
//    SmartField tmp67;
//    SmartField tmp68;
//    SmartField tmp69;
//    SmartField tmp70;
//    SmartField tmp71;
//    SmartField tmp72;
//    SmartField tmp73;
//    SmartField tmp74;
//    SmartField tmp75;
//    SmartField tmp76;
//    SmartField tmp77;
//    SmartField tmp78;
//    SmartField tmp79;
//    SmartField tmp80;
//    SmartField tmp81;
//    SmartField tmp82;
//    SmartField tmp83;
//    SmartField tmp84;
//    SmartField tmp85;
//    SmartField tmp86;
//    SmartField tmp87;
//    SmartField tmp88;
//    SmartField tmp89;
//    SmartField tmp90;
//    SmartField tmp91;
//    SmartField tmp92;
//    SmartField tmp93;
//    SmartField tmp94;
//    SmartField tmp95;
//    SmartField tmp96;
//    SmartField tmp97;
//    SmartField tmp98;
//    SmartField tmp99;

  } // }}}

  gettimeofday(&endtime, NULL); // time of day
  firstalloc_time = (midtime.tv_sec - begtime.tv_sec) + (midtime.tv_usec - begtime.tv_usec)/1000000.0;
  otheralloc_time = (endtime.tv_sec - midtime.tv_sec) + (endtime.tv_usec - midtime.tv_usec)/1000000.0;
  totaltime = (endtime.tv_sec - begtime.tv_sec) + (endtime.tv_usec - begtime.tv_usec)/1000000.0;

  printf(" First Allocation: %22.15e (sec.)\n", firstalloc_time );
  printf(" Next 1K Allocations: %22.15e (sec.)\n", otheralloc_time );
  printf(" Total Time: %22.15e (sec.)\n\n", totaltime );
  std::cout << "\n *** End SmartField Tests ***\n";

} // }}}

