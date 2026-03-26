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

void unit_tests( Grid* _CurrGrid );
void time_tests( Grid* _CurrGrid );
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

  unit_tests( CurrGrid );
  //time_tests(CurrGrid );
  
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


void unit_tests( Grid* _CurrGrid )
{ // {{{

  std::cout << "\n *** SmartField Unit Tests ***\n";
  std::cout << "   (All results should be zero)\n\n";
  { // {{{

    FieldType ftnum = 1.;

    FieldStack initstack;
    Field<FieldType> myfield(*(initstack.CurrLayout), true, "");
    myfield = FieldType(2.);

    std::vector< FieldType > myvec( SmartField::NPW );
    for (UInt m =0; m<SmartField::NPW; m++)
    {
      myvec[m] = FieldType(3*(m+1));
    }
   
    std::vector< std::vector<FieldType> > myvecvec( SmartField::FDSize, myvec );
    for (UInt m =0; m<SmartField::FDSize; m++)
    {
      myvecvec[m] = std::vector<FieldType>( SmartField::NPW, FieldType(4*(m+1)) );
    }

    std::vector< Field<FieldType> > myfieldvec( SmartField::FDSize, myfield );
    for (UInt m =0; m<SmartField::FDSize; m++)
    {
      myfieldvec[m] = FieldType(5*(m+1));
    }

    // assignment operators
    std::cout << " * Assignment operators\n";
    {
      SmartField tmp1, tmp2;

      tmp1 = ftnum;
      std::cout << std::setw(45) << std::left << "  - SmartField = constant";
      std::cout << tmp1[0].getelement(1) - ftnum << std::endl;
      
      tmp1 = myvec;
      std::cout << std::setw(45) << std::left << "  - SmartField = stl vector";
      std::cout << tmp1[0].getelement(1) - myvec[1] << std::endl;
      
      tmp1 = myfield;
      std::cout << std::setw(45) << std::left << "  - SmartField = Field";
      std::cout << tmp1[0].getelement(1) - myfield.getelement(1) << std::endl;

      tmp1 = myvecvec;
      std::cout << std::setw(45) << std::left << "  - SmartField = stl vector of vectors";
      std::cout << tmp1[0].getelement(1) - myvecvec[0][1] << std::endl;

      tmp1 = myfieldvec;
      std::cout << std::setw(45) << std::left << "  - SmartField = stl vector of Field";
      std::cout << tmp1[0].getelement(1) - myfieldvec[0].getelement(1) << std::endl;
      
      tmp2 = tmp1;
      std::cout << std::setw(45) << std::left << "  - SmartField = SmartField";
      std::cout << tmp2[0].getelement(1) - tmp1[0].getelement(1) << std::endl;
    }

    // compound subtraction
    std::cout << " * Compound subtraction operators\n";
    {
      SmartField tmp1 , tmp2;
      
      tmp1 = ftnum;
      tmp1 -= ftnum;
      std::cout << std::setw(45) << std::left << "  - SmartField = constant";
      std::cout << tmp1[0].getelement(1) << std::endl;
      
      tmp1 = myfield;
      tmp1 -= myfield;
      std::cout << std::setw(45) << std::left << "  - SmartField = Field";
      std::cout << tmp1[0].getelement(1) << std::endl;
      
      tmp1 = myfieldvec;
      tmp1 -= myfieldvec;
      std::cout << std::setw(45) << std::left << "  - SmartField = stl vector of Field";
      std::cout << tmp1[0].getelement(1) << std::endl;
      
      tmp2 = myfieldvec;
      tmp2 = tmp1;
      tmp1 -= tmp2;
      std::cout << std::setw(45) << std::left << "  - SmartField = SmartField";
      std::cout << tmp1[0].getelement(1) << std::endl;
    }
    
    // compound addition
    std::cout << " * Compound addition operators\n";
    {
      SmartField tmp1 , tmp2;
      
      tmp1 = FieldType(0.);
      tmp1 -= ftnum;
      tmp1 += ftnum;
      std::cout << std::setw(45) << std::left << "  - SmartField = constant";
      std::cout << tmp1[0].getelement(1) << std::endl;
      
      tmp1 = FieldType(0.);
      tmp1 -= myfield;
      tmp1 += myfield;
      std::cout << std::setw(45) << std::left << "  - SmartField = Field";
      std::cout << tmp1[0].getelement(1) << std::endl;
      
      tmp1 = FieldType(0.);
      tmp1 -= myfieldvec;
      tmp1 += myfieldvec;
      std::cout << std::setw(45) << std::left << "  - SmartField = stl vector of Field";
      std::cout << tmp1[0].getelement(1) << std::endl;
      
      tmp2 = FieldType(3.);
      tmp1 = FieldType(0.);
      tmp1 -= tmp2;
      tmp1 += tmp2;
      std::cout << std::setw(45) << std::left << "  - SmartField = SmartField";
      std::cout << tmp1[0].getelement(1) << std::endl;
    }
    
    // compound multiplication
    std::cout << " * Compound multiplication operators\n";
    {
      SmartField tmp1 , tmp2;
      
      tmp1 = FieldType(3.);
      tmp1 *= ftnum;
      tmp1 -= FieldType(3.);
      std::cout << std::setw(45) << std::left << "  - SmartField = constant";
      std::cout << tmp1[0].getelement(1) << std::endl;
      
      tmp1 = FieldType(3.);
      tmp1 *= myfield;
      tmp1 -= FieldType(6.);
      std::cout << std::setw(45) << std::left << "  - SmartField = Field";
      std::cout << tmp1[0].getelement(1) << std::endl;
      
      tmp1 = FieldType(3.);
      tmp1 *= myfieldvec;
      for (UInt m=0; m<myfieldvec.size(); m++) tmp1[m] += FieldType( -3.*5.*(m+1.) );
      std::cout << std::setw(45) << std::left << "  - SmartField = stl vector of Field";
      std::cout << tmp1[0].getelement(1) << std::endl;
      
      tmp2 = FieldType(3.);
      tmp1 = FieldType(3.);
      tmp1 *= tmp2;
      tmp1 -= FieldType(9.);
      std::cout << std::setw(45) << std::left << "  - SmartField = SmartField";
      std::cout << tmp1[0].getelement(1) << std::endl;
    }
    
    // compound division
    std::cout << " * Compound division operators\n";
    {
      SmartField tmp1 , tmp2;
      
      tmp1 = FieldType(3.);
      tmp1 /= ftnum;
      tmp1 -= FieldType(3.);
      std::cout << std::setw(45) << std::left << "  - SmartField = constant";
      std::cout << tmp1[0].getelement(1) << std::endl;
      
      tmp1 = FieldType(3.);
      tmp1 /= myfield;
      tmp1 -= FieldType(1.5);
      std::cout << std::setw(45) << std::left << "  - SmartField = Field";
      std::cout << tmp1[0].getelement(1) << std::endl;

      tmp1 = FieldType(3.);
      tmp1 /= myfieldvec;
      for (UInt m=0; m<myfieldvec.size(); m++) tmp1[m] += FieldType( -3./(5.*(m+1)) );
      std::cout << std::setw(45) << std::left << "  - SmartField = stl vector of Field";
      std::cout << tmp1[0].getelement(1) << std::endl;
      
      tmp2 = FieldType(4.);
      tmp1 = FieldType(3.);
      tmp1 /= tmp2;
      tmp1 -= FieldType(0.75);
      std::cout << std::setw(45) << std::left << "  - SmartField = SmartField";
      std::cout << tmp1[0].getelement(1) << std::endl;     
    }

    std::cout << " * File IO\n";
    {
      SmartField tmp1, tmp2;
      tmp1 = myfieldvec;

      std::string fname = "tmp1.dat";
      tmp1.writefield(fname, true);

      tmp2.readfield("phi1.in");
      tmp1 += tmp2;

      fname = "tmp1_v2.dat";
      tmp1.writecplxfield(fname);
    }

    std::cout << " * Getters and Setters\n";
    {
      SmartField tmp1;
      tmp1 = myfieldvec;

      UInt idx1 = SmartField::NPW/3+1;
      UInt idx2 = idx1+1;
      UInt idx3 = idx2+1;

      tmp1.setelement(FieldType(1024.), idx1);
      tmp1.setelement(FieldType(-8.), idx2);
      tmp1.setelement(FieldType(-2048.), idx3);
      std::cout << "  - element " << idx1 << ": " << tmp1.getelement( idx1 ) - FieldType(1024.) << std::endl;
      std::cout << "  - element " << idx2 << ": " << tmp1.getelement( idx2 ) + FieldType(8.0) << std::endl;
      std::cout << "  - l1norm (nonzero): " << tmp1.l1norm() << std::endl;
      std::cout << "  - max: " << tmp1.max() + FieldType(2048.) << std::endl;
      std::cout << "  - maxsigned: " << tmp1.maxsigned() - FieldType(1024) << std::endl;
      std::cout << "  - minsigned: " << tmp1.minsigned() + FieldType(2048) << std::endl;
      std::cout << "  - maxidx: " << tmp1.maxidx() - idx1 << std::endl;
      std::cout << "  - sum (nonzero): " << tmp1.sumelem() << std::endl;
      std::cout << "  - integrate (nonzero): " << tmp1.integrate() << std::endl;
      std::cout << "  - integrate square (nonzero): " << tmp1.integrate_square() << std::endl;

      tmp1.log(tmp1);
      std::cout << "  - log: " << tmp1[0].getelement(1) << std::endl;
      tmp1.exponentiate(tmp1);
      std::cout << "  - exp: " << tmp1[0].getelement(1) << std::endl;
      tmp1.sqrt(tmp1);
      std::cout << "  - sqrt: " << tmp1[0].getelement(1) << std::endl;

      std::cout << "  - flag_inrealspace: " << tmp1.getflag_inrealspace() << std::endl;
      tmp1.setflag_inrealspace(false);
      std::cout << "  - flag_inrealspace: " << tmp1.getflag_inrealspace() << std::endl;
      tmp1.setflag_inrealspace(true);

      SmartField tmp2 = FieldType(3.0, 4.0);
      std::cout << "  - l1norm: " << tmp2.l1norm() << std::endl;
      std::cout << "  - l2norm: " << tmp2.l2norm() << std::endl;

      std::cout << "  - array copy (before): " << myvecvec[0][1] << std::endl;
      tmp2.copytoarray( myvecvec );
      std::cout << "  - array copy (after): " << myvecvec[0][1] << std::endl;

    }

  } // }}}
  std::cout << "\n *** End SmartField Unit Tests ***\n";

  std::cout << "\n\n ***Testing SmartFieldMat Matrix Operations***\n";
  std::cout << "   (All results should be zero)\n\n";
  { // {{{

    SmartFieldMat A( 3, 3 );
    SmartFieldMat B( 3, 3 );
    SmartFieldMat C( 3, 3 );
    SmartFieldVec b( 3 );
    SmartFieldVec x( 3 );
    FieldType ftnum(0.);
    SmartField myfield(0.);
    Vec2dFieldType Simple_A(3, Vec1dFieldType(3, 0.));
    Vec1dFieldType Simple_x(3, 0.);

    A.zero();
    B.zero();
    C.zero();
    b.zero();
    x.zero();

    A.setflag_inrealspace(true);
    B.setflag_inrealspace(true);
    C.setflag_inrealspace(true);
    b.setflag_inrealspace(true);
    x.setflag_inrealspace(true);

    // assignment
    std::cout << " * Assignment operators\n";
    { 

      A(0,0) = FieldType(1.);
      A(0,1) = FieldType(2.);
      A(0,2) = FieldType(3.);
      A(1,0) = FieldType(4.);
      A(1,1) = FieldType(5.);
      A(1,2) = FieldType(4.);
      A(2,0) = FieldType(3.);
      A(2,1) = FieldType(2.);
      A(2,2) = FieldType(1.);

      B(0,0) = FieldType(2.);
      B(0,1) = FieldType(5.);
      B(0,2) = FieldType(-3.);
      B(1,0) = FieldType(7.);
      B(1,1) = FieldType(-4.);
      B(1,2) = FieldType(8.);
      B(2,0) = FieldType(0.);
      B(2,1) = FieldType(-1.3);
      B(2,2) = FieldType(-4);

      Simple_A[0][0] = FieldType(1.);
      Simple_A[0][1] = FieldType(2.);
      Simple_A[0][2] = FieldType(3.);
      Simple_A[1][0] = FieldType(4.);
      Simple_A[1][1] = FieldType(5.);
      Simple_A[1][2] = FieldType(4.);
      Simple_A[2][0] = FieldType(3.);
      Simple_A[2][1] = FieldType(2.);
      Simple_A[2][2] = FieldType(1.);

      myfield = FieldType(2.);
      ftnum = 1.;

      A = B;
      std::cout << std::setw(45) << std::left << "  - SmartFieldMat = SmartFieldMat";
      std::cout << A(0,0).getelement(1) - B(0,0).getelement(1) << std::endl;

      A = Simple_A;
      std::cout << std::setw(45) << std::left << "  - SmartFieldMat = 2D constant vector";
      std::cout << A(0,0).getelement(1) - Simple_A[0][0] << std::endl;

      A = myfield;
      std::cout << std::setw(45) << std::left << "  - SmartFieldMat = SmartField scalar";
      std::cout << A(0,0).getelement(1) - myfield.getelement(1) << std::endl;

      A = ftnum;
      std::cout << std::setw(45) << std::left << "  - SmartFieldMat = constant scalar";
      std::cout << A(0,0).getelement(1) - ftnum << std::endl;

      A.zero();
      std::cout << std::setw(45) << std::left << "  - SmartFieldMat.zero()" << std::endl << "\t";
      std::cout << A(0,0).getelement(0) << ", ";
      std::cout << A(0,1).getelement(0) << ", ";
      std::cout << A(0,2).getelement(0) << std::endl << "\t"; 
      std::cout << A(1,0).getelement(0) << ", ";
      std::cout << A(1,1).getelement(0) << ", ";
      std::cout << A(1,2).getelement(0) << std::endl << "\t"; 
      std::cout << A(2,0).getelement(0) << ", ";
      std::cout << A(2,1).getelement(0) << ", ";
      std::cout << A(2,2).getelement(0) << std::endl << "\t";

      A.eye(FieldType(1.));
      std::cout << std::setw(45) << std::left << "  - SmartFieldMat.eye()" << std::endl << "\t";
      std::cout << A(0,0).getelement(0) - FieldType(1.) << ", ";
      std::cout << A(0,1).getelement(0) << ", ";
      std::cout << A(0,2).getelement(0) << std::endl << "\t"; 
      std::cout << A(1,0).getelement(0) << ", ";
      std::cout << A(1,1).getelement(0) - FieldType(1.) << ", ";
      std::cout << A(1,2).getelement(0) << std::endl << "\t"; 
      std::cout << A(2,0).getelement(0) << ", ";
      std::cout << A(2,1).getelement(0) << ", ";
      std::cout << A(2,2).getelement(0) - FieldType(1.) << std::endl << "\t";

    }

    // addition
    std::cout << " * Addition operators\n";
    { 

      B(0,0) = FieldType(2.);
      B(0,1) = FieldType(5.);
      B(0,2) = FieldType(-3.);
      B(1,0) = FieldType(7.);
      B(1,1) = FieldType(-4.);
      B(1,2) = FieldType(8.);
      B(2,0) = FieldType(0.);
      B(2,1) = FieldType(-1.3);
      B(2,2) = FieldType(-4);

      Simple_A[0][0] = FieldType(1.);
      Simple_A[0][1] = FieldType(2.);
      Simple_A[0][2] = FieldType(3.);
      Simple_A[1][0] = FieldType(4.);
      Simple_A[1][1] = FieldType(5.);
      Simple_A[1][2] = FieldType(4.);
      Simple_A[2][0] = FieldType(3.);
      Simple_A[2][1] = FieldType(2.);
      Simple_A[2][2] = FieldType(1.);

      myfield = FieldType(2.);
      ftnum = 1.;

      A.zero();
      A += B;
      std::cout << std::setw(45) << std::left << "  - SmartFieldMat += SmartFieldMat";
      std::cout << A(0,0).getelement(1) - B(0,0).getelement(1) << std::endl;

      A.zero();
      A += Simple_A;
      std::cout << std::setw(45) << std::left << "  - SmartFieldMat += 2D constant vector";
      std::cout << A(0,0).getelement(1) - Simple_A[0][0] << std::endl;

      A.zero();
      A += myfield;
      std::cout << std::setw(45) << std::left << "  - SmartFieldMat += SmartField scalar";
      std::cout << A(0,0).getelement(1) - myfield.getelement(1) << std::endl;

      A.zero();
      A += ftnum;
      std::cout << std::setw(45) << std::left << "  - SmartFieldMat += constant scalar";
      std::cout << A(0,0).getelement(1) - ftnum << std::endl;
    }

    // subtraction 
    std::cout << " * Subtraction operators\n";
    { 

      B(0,0) = FieldType(2.);
      B(0,1) = FieldType(5.);
      B(0,2) = FieldType(-3.);
      B(1,0) = FieldType(7.);
      B(1,1) = FieldType(-4.);
      B(1,2) = FieldType(8.);
      B(2,0) = FieldType(0.);
      B(2,1) = FieldType(-1.3);
      B(2,2) = FieldType(-4);

      Simple_A[0][0] = FieldType(1.);
      Simple_A[0][1] = FieldType(2.);
      Simple_A[0][2] = FieldType(3.);
      Simple_A[1][0] = FieldType(4.);
      Simple_A[1][1] = FieldType(5.);
      Simple_A[1][2] = FieldType(4.);
      Simple_A[2][0] = FieldType(3.);
      Simple_A[2][1] = FieldType(2.);
      Simple_A[2][2] = FieldType(1.);

      myfield = FieldType(2.);
      ftnum = 1.;

      A.zero();
      A -= B;
      std::cout << std::setw(45) << std::left << "  - SmartFieldMat -= SmartFieldMat";
      std::cout << A(0,0).getelement(1) + B(0,0).getelement(1) << std::endl;

      A.zero();
      A -= Simple_A;
      std::cout << std::setw(45) << std::left << "  - SmartFieldMat -= 2D constant vector";
      std::cout << A(0,0).getelement(1) + Simple_A[0][0] << std::endl;

      A.zero();
      A -= myfield;
      std::cout << std::setw(45) << std::left << "  - SmartFieldMat -= SmartField scalar";
      std::cout << A(0,0).getelement(1) + myfield.getelement(1) << std::endl;

      A.zero();
      A -= ftnum;
      std::cout << std::setw(45) << std::left << "  - SmartFieldMat -= constant scalar";
      std::cout << A(0,0).getelement(1) + ftnum << std::endl;
    }

    // multiplication
    std::cout << " * Multiplication operators\n";
    { 

      A(0,0) = FieldType(1.);
      A(0,1) = FieldType(2.);
      A(0,2) = FieldType(3.);
      A(1,0) = FieldType(4.);
      A(1,1) = FieldType(5.);
      A(1,2) = FieldType(4.);
      A(2,0) = FieldType(3.);
      A(2,1) = FieldType(2.);
      A(2,2) = FieldType(1.);

      B(0,0) = FieldType(2.);
      B(0,1) = FieldType(5.);
      B(0,2) = FieldType(-3.);
      B(1,0) = FieldType(7.);
      B(1,1) = FieldType(-4.);
      B(1,2) = FieldType(8.);
      B(2,0) = FieldType(0.);
      B(2,1) = FieldType(-1.3);
      B(2,2) = FieldType(-4);

      Simple_A[0][0] = FieldType(1.);
      Simple_A[0][1] = FieldType(2.);
      Simple_A[0][2] = FieldType(3.);
      Simple_A[1][0] = FieldType(4.);
      Simple_A[1][1] = FieldType(5.);
      Simple_A[1][2] = FieldType(4.);
      Simple_A[2][0] = FieldType(3.);
      Simple_A[2][1] = FieldType(2.);
      Simple_A[2][2] = FieldType(1.);

      myfield = FieldType(2.);
      ftnum = 3.;

      A *= B;
      std::cout << std::setw(45) << std::left << "  - SmartFieldMat *= SmartFieldMat";
      std::cout << A(0,0).getelement(1) - Simple_A[0][0]*B(0,0).getelement(1) << std::endl;

      A = Simple_A;
      A *= Simple_A;
      std::cout << std::setw(45) << std::left << "  - SmartFieldMat *= 2D constant vector";
      std::cout << A(0,0).getelement(1) - Simple_A[0][0]*Simple_A[0][0] << std::endl;

      A = Simple_A;
      A *= myfield;
      std::cout << std::setw(45) << std::left << "  - SmartFieldMat *= SmartField scalar";
      std::cout << A(0,0).getelement(1) - myfield.getelement(1)*Simple_A[0][0] << std::endl;

      A = Simple_A;
      A *= ftnum;
      std::cout << std::setw(45) << std::left << "  - SmartFieldMat *= constant scalar";
      std::cout << A(0,0).getelement(1) - Simple_A[0][0]*ftnum << std::endl;
    }

    // dot product 
    std::cout << " * Dot Product\n";
    { 

      A(0,0) = FieldType(1.);
      A(0,1) = FieldType(2.);
      A(0,2) = FieldType(3.);
      A(1,0) = FieldType(4.);
      A(1,1) = FieldType(5.);
      A(1,2) = FieldType(4.);
      A(2,0) = FieldType(3.);
      A(2,1) = FieldType(2.);
      A(2,2) = FieldType(1.);

      B(0,0) = FieldType(2.);
      B(0,1) = FieldType(5.);
      B(0,2) = FieldType(-3.);
      B(1,0) = FieldType(7.);
      B(1,1) = FieldType(-4.);
      B(1,2) = FieldType(8.);
      B(2,0) = FieldType(0.);
      B(2,1) = FieldType(-1.3);
      B(2,2) = FieldType(-4);

      x[0] = FieldType(1.);
      x[1] = FieldType(2.);
      x[2] = FieldType(3.);

      C.dot(A, B); // b[0] = 14; b[1] = 26; b[2] = 10;
      std::cout << std::setw(45) << std::left << "  - SmartFieldMat .dot. SmartFieldMat" << std::endl << "\t";
      std::cout << C(0,0).getelement(0) - FieldType(16.) << ", ";
      std::cout << C(0,1).getelement(0) - FieldType(-6.9) << ", ";
      std::cout << C(0,2).getelement(0) - FieldType(1.) << std::endl << "\t"; 
      std::cout << C(1,0).getelement(0) - FieldType(43.) << ", ";
      std::cout << C(1,1).getelement(0) - FieldType(-5.2) << ", ";
      std::cout << C(1,2).getelement(0) - FieldType(12.) << std::endl << "\t"; 
      std::cout << C(2,0).getelement(0) - FieldType(20.) << ", ";
      std::cout << C(2,1).getelement(0) - FieldType(5.7) << ", ";
      std::cout << C(2,2).getelement(0) - FieldType(3.) << std::endl;

      b.dot(A, x); // b[0] = 14; b[1] = 26; b[2] = 10;
      std::cout << std::setw(45) << std::left << "  - SmartFieldMat .dot. SmartFieldVec" << std::endl << "\t";
      std::cout << b[0].getelement(0) - FieldType(14.) << ", ";
      std::cout << b[1].getelement(0) - FieldType(26.) << ", ";
      std::cout << b[2].getelement(0) - FieldType(10.) << std::endl;

      b.dot(x, A); // b[0] = 18; b[1] = 18; b[2] = 14;
      std::cout << std::setw(45) << std::left << "  - SmartFieldVec .dot. SmartFieldMat"<<std::endl<<"\t";
      std::cout << b[0].getelement(0) - FieldType(18.) << ", ";
      std::cout << b[1].getelement(0) - FieldType(18.) << ", ";
      std::cout << b[2].getelement(0) - FieldType(14.) << std::endl;

      // ** scalar product isn't written ** 
      // b[0] = FieldType(18.);
      // b[1] = FieldType(18.);
      // b[2] = FieldType(14.);
      // myfield.dot(x, b); // 
      // std::cout << std::setw(45) << std::left << "  - SmartFieldVec .dot. SmartFieldVec";
      // std::cout << myfield.getelement(0) - FieldType(96.) << std::endl;

    }

    std::cout << " * Matrix Ops\n";
    {

      A(0,0) = FieldType(1.);
      A(0,1) = FieldType(2.);
      A(0,2) = FieldType(3.);
      A(1,0) = FieldType(4.);
      A(1,1) = FieldType(5.);
      A(1,2) = FieldType(4.);
      A(2,0) = FieldType(3.);
      A(2,1) = FieldType(2.);
      A(2,2) = FieldType(1.);

      x[0] = FieldType(1.);
      x[1] = FieldType(2.);
      x[2] = FieldType(3.);

      b[0] = FieldType(14.);
      b[1] = FieldType(26.);
      b[2] = FieldType(10.);

      myfield.zero();

      // matrix norm
      A.norm(myfield);
      std::cout << std::setw(45) << std::left << "  - SmartFieldMat.norm() ";
      std::cout << myfield.getelement(0) - FieldType(9.2195444572928871) << std::endl;

      // outer product
      B.outer(x, b);
      std::cout << std::setw(45) << std::left << "  - SmartFieldVec .outer_prod. SmartFieldVec" << std::endl << "\t";
      std::cout << B(0,0).getelement(0) - FieldType(14.) << ", ";
      std::cout << B(0,1).getelement(0) - FieldType(26.) << ", ";
      std::cout << B(0,2).getelement(0) - FieldType(10.) << std::endl << "\t"; 
      std::cout << B(1,0).getelement(0) - FieldType(28.) << ", ";
      std::cout << B(1,1).getelement(0) - FieldType(52) << ", ";
      std::cout << B(1,2).getelement(0) - FieldType(20.) << std::endl << "\t"; 
      std::cout << B(2,0).getelement(0) - FieldType(42.) << ", ";
      std::cout << B(2,1).getelement(0) - FieldType(78.) << ", ";
      std::cout << B(2,2).getelement(0) - FieldType(30.) << std::endl;

      // in-place transpose
      B.transpose();
      std::cout << std::setw(45) << std::left << "  - SmartFieldMat.transpose()" << std::endl << "\t";
      std::cout << B(0,0).getelement(0) - FieldType(14.) << ", ";
      std::cout << B(0,1).getelement(0) - FieldType(28.) << ", ";
      std::cout << B(0,2).getelement(0) - FieldType(42.) << std::endl << "\t";
      std::cout << B(1,0).getelement(0) - FieldType(26.) << ", ";
      std::cout << B(1,1).getelement(0) - FieldType(52.) << ", ";
      std::cout << B(1,2).getelement(0) - FieldType(78.) << std::endl << "\t";
      std::cout << B(2,0).getelement(0) - FieldType(10.) <<  ", ";
      std::cout << B(2,1).getelement(0) - FieldType(20.) <<  ", ";
      std::cout << B(2,2).getelement(0) - FieldType(30.) << std::endl;

    }

    // linear solve
    std::cout << " * Linear Solve\n";
    { 

      A(0,0) = FieldType(1.);
      A(0,1) = FieldType(2.);
      A(0,2) = FieldType(3.);
      A(1,0) = FieldType(4.);
      A(1,1) = FieldType(5.);
      A(1,2) = FieldType(4.);
      A(2,0) = FieldType(3.);
      A(2,1) = FieldType(2.);
      A(2,2) = FieldType(1.);

      b[0] = FieldType(14.);
      b[1] = FieldType(26.);
      b[2] = FieldType(10.);

      x.linsolve(A, b); // a[0] = 1; a[1] = 2; a[2] = 3;
      std::cout << std::setw(45) << std::left << "  - linsolve( SmartFieldMat, SmartFieldVec) " << std::endl << "\t";
      std::cout << x[0].getelement(0) - FieldType(1.) << ", ";
      std::cout << x[1].getelement(0) - FieldType(2.) << ", ";
      std::cout << x[2].getelement(0) - FieldType(3.) << std::endl;

      A(0,0) = FieldType(1.);
      A(0,1) = FieldType(2.);
      A(0,2) = FieldType(3.);
      A(1,0) = FieldType(4.);
      A(1,1) = FieldType(5.);
      A(1,2) = FieldType(4.);
      A(2,0) = FieldType(3.);
      A(2,1) = FieldType(2.);
      A(2,2) = FieldType(1.);

      B.invert( A ); // A^{-1} = 1/8 * [[3, -4, 7], [-8, 8, -8], [7, -4, 3]]
      std::cout << std::setw(45) << std::left << "  - matrix invert ( SmartFieldMat )" << std::endl << "\t";
      std::cout << B(0,0).getelement(0) - FieldType(3./8.) << ", ";
      std::cout << B(0,1).getelement(0) - FieldType(-1./2.) << ", ";
      std::cout << B(0,2).getelement(0) - FieldType(7./8.) << std::endl << "\t"; 
      std::cout << B(1,0).getelement(0) - FieldType(-1.) << ", ";
      std::cout << B(1,1).getelement(0) - FieldType(1.) << ", ";
      std::cout << B(1,2).getelement(0) - FieldType(-1.) << std::endl << "\t"; 
      std::cout << B(2,0).getelement(0) - FieldType(7./8.) << ", ";
      std::cout << B(2,1).getelement(0) - FieldType(-1./2.) << ", ";
      std::cout << B(2,2).getelement(0) - FieldType(3./8.) << std::endl;

    }

    // getters and setters
    std::cout << " * Getters, Setters and FFTs\n";
    {

      A(0,0) = FieldType(1.);
      A(0,1) = FieldType(2.);
      A(0,2) = FieldType(3.);
      A(1,0) = FieldType(4.);
      A(1,1) = FieldType(5.);
      A(1,2) = FieldType(4.);
      A(2,0) = FieldType(3.);
      A(2,1) = FieldType(2.);
      A(2,2) = FieldType(1.);
    
      bool realspace;
  
      realspace = A.getflag_inrealspace();
      std::cout << std::setw(45) << std::left << "  - SmartFieldMat.getflag_inrealspace()";
      std::cout << int(realspace) - 1 << std::endl;

      A.setflag_inrealspace(false);
      realspace = A.getflag_inrealspace();
      std::cout << std::setw(45) << std::left << "  - SmartFieldMat.setflag_inrealspace()";
      std::cout << int(realspace) << std::endl;

      A.setflag_inrealspace(true);
      Simple_A = A.getelement(1);
      std::cout << std::setw(45) << std::left << "  - SmartFieldMat.getelement(idx)" << std::endl << "\t";
      std::cout << Simple_A[0][0] - FieldType(1.) << ", ";
      std::cout << Simple_A[0][1] - FieldType(2.) << ", ";
      std::cout << Simple_A[0][2] - FieldType(3.) << std::endl << "\t"; 
      std::cout << Simple_A[1][0] - FieldType(4.) << ", ";
      std::cout << Simple_A[1][1] - FieldType(5.) << ", ";
      std::cout << Simple_A[1][2] - FieldType(4.) << std::endl << "\t"; 
      std::cout << Simple_A[2][0] - FieldType(3.) << ", ";
      std::cout << Simple_A[2][1] - FieldType(2.) << ", ";
      std::cout << Simple_A[2][2] - FieldType(1.) << std::endl;

      A.fft();
      // element(0) == constant value
      // elements(1-end) == 0
      std::cout << std::setw(45) << std::left << "  - SmartFieldMat.fft()" << std::endl << "\t";
      std::cout << A(0,0).getelement(1) << ", ";
      std::cout << A(0,1).getelement(1) << ", ";
      std::cout << A(0,2).getelement(1) << std::endl << "\t"; 
      std::cout << A(1,0).getelement(1) << ", ";
      std::cout << A(1,1).getelement(1) << ", ";
      std::cout << A(1,2).getelement(1) << std::endl << "\t"; 
      std::cout << A(2,0).getelement(1) << ", ";
      std::cout << A(2,1).getelement(1) << ", ";
      std::cout << A(2,2).getelement(1) << std::endl;

      A.ifft();
      std::cout << std::setw(45) << std::left << "  - SmartFieldMat.ifft()" << std::endl << "\t";
      std::cout << A(0,0).getelement(1) - FieldType(1.) << ", ";
      std::cout << A(0,1).getelement(1) - FieldType(2.) << ", ";
      std::cout << A(0,2).getelement(1) - FieldType(3.) << std::endl << "\t"; 
      std::cout << A(1,0).getelement(1) - FieldType(4.) << ", ";
      std::cout << A(1,1).getelement(1) - FieldType(5.) << ", ";
      std::cout << A(1,2).getelement(1) - FieldType(4.) << std::endl << "\t"; 
      std::cout << A(2,0).getelement(1) - FieldType(3.) << ", ";
      std::cout << A(2,1).getelement(1) - FieldType(2.) << ", ";
      std::cout << A(2,2).getelement(1) - FieldType(1.) << std::endl;

    }

    // subscript and slice (only do if FDSize > 1)
    std::cout << " * Subscript and Slice (for Finite Differences if FDSize > 1)\n";
    std::cout << "  - FDSize = " << _CurrGrid->FDSize() << std::endl;
    if (_CurrGrid->FDSize() > 1)
    {

      A.zero();
      A(0,0)[1] = FieldType(1.);
      A(0,1)[1] = FieldType(2.);
      A(0,2)[1] = FieldType(3.);
      A(1,0)[1] = FieldType(4.);
      A(1,1)[1] = FieldType(5.);
      A(1,2)[1] = FieldType(4.);
      A(2,0)[1] = FieldType(3.);
      A(2,1)[1] = FieldType(2.);
      A(2,2)[1] = FieldType(1.);

      A(0,0)[3] = FieldType(1.);
      A(0,1)[3] = FieldType(2.);
      A(0,2)[3] = FieldType(3.);
      A(1,0)[3] = FieldType(4.);
      A(1,1)[3] = FieldType(5.);
      A(1,2)[3] = FieldType(4.);
      A(2,0)[3] = FieldType(3.);
      A(2,1)[3] = FieldType(2.);
      A(2,2)[3] = FieldType(1.);


      B.zero();

      int Nx = _CurrGrid->FDSize(); // need an int to do modulus

      std::vector<UInt> idx(Nx,0);
      for (int i=0; i<Nx; i++)
      {
        idx[i] = ((i+1+Nx)%Nx);
      }
      
      // Shifts elements of the FD array by -1.
      B.subscript(A, idx);
      std::cout << std::setw(45) << std::left << "  - subscript( SmartFieldMat, idx ) " << std::endl << "\t";
      std::cout << B(0,0)[0].getelement(0) - FieldType(1.) << ", ";
      std::cout << B(0,1)[0].getelement(0) - FieldType(2.) << ", ";
      std::cout << B(0,2)[0].getelement(0) - FieldType(3.) << std::endl << "\t"; 
      std::cout << B(1,0)[0].getelement(0) - FieldType(4.) << ", ";
      std::cout << B(1,1)[0].getelement(0) - FieldType(5.) << ", ";
      std::cout << B(1,2)[0].getelement(0) - FieldType(4.) << std::endl << "\t"; 
      std::cout << B(2,0)[0].getelement(0) - FieldType(3.) << ", ";
      std::cout << B(2,1)[0].getelement(0) - FieldType(2.) << ", ";
      std::cout << B(2,2)[0].getelement(0) - FieldType(1.) << std::endl;

      std::vector<UInt> idx2(3, 0);
      idx2[0] = 1;
      idx2[1] = 3;
      idx2[2] = 5;
      
      B.zero();
      B.slice(A, idx);
      std::cout << std::setw(45) << std::left << "  - slice( SmartFieldMat, idx ) " << std::endl << "\t";
      std::cout << B(0,0)[3].getelement(0) - FieldType(1.) << ", ";
      std::cout << B(0,1)[3].getelement(0) - FieldType(2.) << ", ";
      std::cout << B(0,2)[3].getelement(0) - FieldType(3.) << std::endl << "\t"; 
      std::cout << B(1,0)[3].getelement(0) - FieldType(4.) << ", ";
      std::cout << B(1,1)[3].getelement(0) - FieldType(5.) << ", ";
      std::cout << B(1,2)[3].getelement(0) - FieldType(4.) << std::endl << "\t"; 
      std::cout << B(2,0)[3].getelement(0) - FieldType(3.) << ", ";
      std::cout << B(2,1)[3].getelement(0) - FieldType(2.) << ", ";
      std::cout << B(2,2)[3].getelement(0) - FieldType(1.) << std::endl;

    }

  } // }}}
  std::cout << "\n ***Finished SmartFieldMat Matrix Operations Test***\n";

  std::cout << "\n\n ***Testing Differential Operators***\n";
  std::cout << "   --2D Periodic Tests only!!--\n";
  std::cout << "   (All results should be zero)\n\n";
  { // {{{

    // --- set up operators (model B only) ---
    // {{{
    int _Ncomp, _NIcomp;
    int _BCFlag; // Which BCs are we using?
    BCs_Factory         *_BCFactory; // for making operators
    Operators_Factory   *_OpFactory; // for making operators
    Operators_Base *_PhiOp; // Operators for Phi fields
    BCs_Base *_PhiBCs; // BCs object

    // set up the BCs and operators for phi
    std::string filename( "params_BCs_phi.in" );
    
    if (jsonFile(filename.c_str()))
    {
      // Read in BC key from file
      std::string tmp_str;
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
    {//legacy file format
      // Read in BC key from file
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
    // }}}

    // --- Create a test field and a bunch of derivatives ---
    // {{{
    RealType x, y;

    int FDSize = _CurrGrid->FDSize();
    int NPW = _CurrGrid->NPW();
    double Lx = _CurrGrid->Lx();
    double Ly = _CurrGrid->Ly();
    double dx = _CurrGrid->Dx();
    double dy = _CurrGrid->Ly()/_CurrGrid->Ny();

    Vec2dFieldType tmp(FDSize, Vec1dFieldType(NPW, FieldType(0.)));

    SmartField a;
    a.setflag_inrealspace(true);
    for( int i = 0; i < FDSize; i++) {
      for( int j = 0; j < NPW; j++) {
        x = dx*(i + floor(j/_CurrGrid->Ny()));
        y = dy*(j%_CurrGrid->Ny());
        tmp[i][j] = sin(2*PI*x/Lx)*sin(2*PI*y/Ly);
      }
    }
    a = tmp;

    SmartField da_dx;
    da_dx.setflag_inrealspace(true);
    for( int i = 0; i < FDSize; i++) {
      for( int j = 0; j < NPW; j++) {
        x = dx*(i + floor(j/_CurrGrid->Ny()));
        y = dy*(j%_CurrGrid->Ny());
        tmp[i][j] = 2*PI/Lx*cos(2*PI*x/Lx)*sin(2*PI*y/Ly);
      }
    }
    da_dx = tmp;

    SmartField da_dy;
    da_dy.setflag_inrealspace(true);
    for( int i = 0; i < FDSize; i++) {
      for( int j = 0; j < NPW; j++) {
        x = dx*(i + floor(j/_CurrGrid->Ny()));
        y = dy*(j%_CurrGrid->Ny());
        tmp[i][j] = 2*PI/Ly*sin(2*PI*x/Lx)*cos(2*PI*y/Ly);
      }
    }
    da_dy = tmp;

    SmartField d2a_dxdy;
    d2a_dxdy.setflag_inrealspace(true);
    for( int i = 0; i < FDSize; i++) {
      for( int j = 0; j < NPW; j++) {
        x = dx*(i + floor(j/_CurrGrid->Ny()));
        y = dy*(j%_CurrGrid->Ny());
        tmp[i][j] = 4*PI*PI/Lx/Ly*cos(2*PI*x/Lx)*cos(2*PI*y/Ly);
      }
    }
    d2a_dxdy = tmp;

    SmartField d2a_dx2;
    d2a_dx2.setflag_inrealspace(true);
    for( int i = 0; i < FDSize; i++) {
      for( int j = 0; j < NPW; j++) {
        x = dx*(i + floor(j/_CurrGrid->Ny()));
        y = dy*(j%_CurrGrid->Ny());
        tmp[i][j] = -4*PI*PI/Lx/Lx*sin(2*PI*x/Lx)*sin(2*PI*y/Ly);
      }
    }
    d2a_dx2 = tmp;

    SmartField d2a_dy2;
    d2a_dy2.setflag_inrealspace(true);
    for( int i = 0; i < FDSize; i++) {
      for( int j = 0; j < NPW; j++) {
        x = dx*(i + floor(j/_CurrGrid->Ny()));
        y = dy*(j%_CurrGrid->Ny());
        tmp[i][j] = -4*PI*PI/Ly/Ly*sin(2*PI*x/Lx)*sin(2*PI*y/Ly);
      }
    }
    d2a_dy2 = tmp;

    SmartField d4a_dx4;
    d4a_dx4.setflag_inrealspace(true);
    for( int i = 0; i < FDSize; i++) {
      for( int j = 0; j < NPW; j++) {
        x = dx*(i + floor(j/_CurrGrid->Ny()));
        y = dy*(j%_CurrGrid->Ny());
        tmp[i][j] = pow(2*PI/Lx, 4)*sin(2*PI*x/Lx)*sin(2*PI*y/Ly);
      }
    }
    d4a_dx4 = tmp;

    SmartField d4a_dy4;
    d4a_dy4.setflag_inrealspace(true);
    for( int i = 0; i < FDSize; i++) {
      for( int j = 0; j < NPW; j++) {
        x = dx*(i + floor(j/_CurrGrid->Ny()));
        y = dy*(j%_CurrGrid->Ny());
        tmp[i][j] = pow(2*PI/Ly, 4)*sin(2*PI*x/Lx)*sin(2*PI*y/Ly);
      }
    }
    d4a_dy4 = tmp;

    SmartField d4a_dx2dy2;
    d4a_dx2dy2.setflag_inrealspace(true);
    for( int i = 0; i < FDSize; i++) {
      for( int j = 0; j < NPW; j++) {
        x = dx*(i + floor(j/_CurrGrid->Ny()));
        y = dy*(j%_CurrGrid->Ny());
        tmp[i][j] = pow(2*PI/Lx, 2)*pow(2*PI/Ly, 2)*sin(2*PI*x/Lx)*sin(2*PI*y/Ly);
      }
    }
    d4a_dx2dy2 = tmp;

    SmartField b;
    b.setflag_inrealspace(true);
    for( int i = 0; i < FDSize; i++) {
      for( int j = 0; j < NPW; j++) {
        x = dx*(i + floor(j/_CurrGrid->Ny()));
        y = dy*(j%_CurrGrid->Ny());
        tmp[i][j] = cos(2*PI*x/Lx)*cos(2*PI*y/Ly);
      }
    }
    b = tmp;

    SmartField db_dx;
    db_dx.setflag_inrealspace(true);
    for( int i = 0; i < FDSize; i++) {
      for( int j = 0; j < NPW; j++) {
        x = dx*(i + floor(j/_CurrGrid->Ny()));
        y = dy*(j%_CurrGrid->Ny());
        tmp[i][j] = -2*PI/Lx*sin(2*PI*x/Lx)*cos(2*PI*y/Ly);
      }
    }
    db_dx = tmp;

    SmartField db_dy;
    db_dy.setflag_inrealspace(true);
    for( int i = 0; i < FDSize; i++) {
      for( int j = 0; j < NPW; j++) {
        x = dx*(i + floor(j/_CurrGrid->Ny()));
        y = dy*(j%_CurrGrid->Ny());
        tmp[i][j] = -2*PI/Ly*cos(2*PI*x/Lx)*sin(2*PI*y/Ly);
      }
    }
    db_dy = tmp;

    SmartField d2b_dxdy;
    d2b_dxdy.setflag_inrealspace(true);
    for( int i = 0; i < FDSize; i++) {
      for( int j = 0; j < NPW; j++) {
        x = dx*(i + floor(j/_CurrGrid->Ny()));
        y = dy*(j%_CurrGrid->Ny());
        tmp[i][j] = 4*PI*PI/Lx/Ly*sin(2*PI*x/Lx)*sin(2*PI*y/Ly);
      }
    }
    d2b_dxdy = tmp;

    SmartField d2b_dx2;
    d2b_dx2.setflag_inrealspace(true);
    for( int i = 0; i < FDSize; i++) {
      for( int j = 0; j < NPW; j++) {
        x = dx*(i + floor(j/_CurrGrid->Ny()));
        y = dy*(j%_CurrGrid->Ny());
        tmp[i][j] = -4*PI*PI/Lx/Lx*cos(2*PI*x/Lx)*cos(2*PI*y/Ly);
      }
    }
    d2b_dx2 = tmp;

    SmartField d2b_dy2;
    d2b_dy2.setflag_inrealspace(true);
    for( int i = 0; i < FDSize; i++) {
      for( int j = 0; j < NPW; j++) {
        x = dx*(i + floor(j/_CurrGrid->Ny()));
        y = dy*(j%_CurrGrid->Ny());
        tmp[i][j] = -4*PI*PI/Ly/Ly*cos(2*PI*x/Lx)*cos(2*PI*y/Ly);
      }
    }
    d2b_dy2 = tmp;

    SmartField d4b_dx4;
    d4b_dx4.setflag_inrealspace(true);
    for( int i = 0; i < FDSize; i++) {
      for( int j = 0; j < NPW; j++) {
        x = dx*(i + floor(j/_CurrGrid->Ny()));
        y = dy*(j%_CurrGrid->Ny());
        tmp[i][j] = pow(2*PI/Lx, 4)*cos(2*PI*x/Lx)*cos(2*PI*y/Ly);
      }
    }
    d4b_dx4 = tmp;

    SmartField d4b_dy4;
    d4b_dy4.setflag_inrealspace(true);
    for( int i = 0; i < FDSize; i++) {
      for( int j = 0; j < NPW; j++) {
        x = dx*(i + floor(j/_CurrGrid->Ny()));
        y = dy*(j%_CurrGrid->Ny());
        tmp[i][j] = pow(2*PI/Ly, 4)*cos(2*PI*x/Lx)*cos(2*PI*y/Ly);
      }
    }
    d4b_dy4 = tmp;

    SmartField d4b_dx2dy2;
    d4b_dx2dy2.setflag_inrealspace(true);
    for( int i = 0; i < FDSize; i++) {
      for( int j = 0; j < NPW; j++) {
        x = dx*(i + floor(j/_CurrGrid->Ny()));
        y = dy*(j%_CurrGrid->Ny());
        tmp[i][j] = pow(2*PI/Lx, 2)*pow(2*PI/Ly, 2)*cos(2*PI*x/Lx)*cos(2*PI*y/Ly);
      }
    }
    d4b_dx2dy2 = tmp;

    Vec2dFieldType m(2, Vec1dFieldType(2, FieldType(0.)));
    m[0][0] = 1; m[0][1] = 0;
    m[1][0] = 0; m[1][1] = 1;
    //m[0][0] = 0.7; m[0][1] = -2;
    //m[1][0] = 5;   m[1][1] = -0.33;

    SmartField del_M_del_ab_1;
    del_M_del_ab_1.setflag_inrealspace(true);
    for( int i = 0; i < FDSize; i++) {
      for( int j = 0; j < NPW; j++) {
        x = dx*(i + floor(j/_CurrGrid->Ny()));
        y = dy*(j%_CurrGrid->Ny());
        tmp[i][j] = m[0][0]*(d2a_dx2[i].getelement(j) + d2a_dy2[i].getelement(j)) 
                  + m[0][1]*(d2b_dx2[i].getelement(j) + d2b_dy2[i].getelement(j));
      }
    }
    del_M_del_ab_1 = tmp;

    SmartField del_M_del_ab_2;
    del_M_del_ab_2.setflag_inrealspace(true);
    for( int i = 0; i < FDSize; i++) {
      for( int j = 0; j < NPW; j++) {
        x = dx*(i + floor(j/_CurrGrid->Ny()));
        y = dy*(j%_CurrGrid->Ny());
        tmp[i][j] = m[1][0]*(d2a_dx2[i].getelement(j) + d2a_dy2[i].getelement(j)) 
                  + m[1][1]*(d2b_dx2[i].getelement(j) + d2b_dy2[i].getelement(j));
      }
    }
    del_M_del_ab_2 = tmp;

    SmartField del_M_del3_ab_1;
    del_M_del3_ab_1.setflag_inrealspace(true);
    for( int i = 0; i < FDSize; i++) {
      for( int j = 0; j < NPW; j++) {
        x = dx*(i + floor(j/_CurrGrid->Ny()));
        y = dy*(j%_CurrGrid->Ny());
        tmp[i][j] = m[0][0]*(d4a_dx4[i].getelement(j) + 2.*d4a_dx2dy2[i].getelement(j) + d4a_dy4[i].getelement(j))
                  + m[0][1]*(d4b_dx4[i].getelement(j) + 2.*d4b_dx2dy2[i].getelement(j) + d4b_dy4[i].getelement(j));
      }
    }
    del_M_del3_ab_1 = tmp;

    SmartField del_M_del3_ab_2;
    del_M_del3_ab_2.setflag_inrealspace(true);
    for( int i = 0; i < FDSize; i++) {
      for( int j = 0; j < NPW; j++) {
        x = dx*(i + floor(j/_CurrGrid->Ny()));
        y = dy*(j%_CurrGrid->Ny());
        tmp[i][j] = m[1][0]*(d4a_dx4[i].getelement(j) + 2.*d4a_dx2dy2[i].getelement(j) + d4a_dy4[i].getelement(j)) 
                  + m[1][1]*(d4b_dx4[i].getelement(j) + 2.*d4b_dx2dy2[i].getelement(j) + d4b_dy4[i].getelement(j));
      }
    }
    del_M_del3_ab_2 = tmp;


    SmartField A;
    SmartField B;
    SmartFieldVec A_vec(2);
    SmartFieldVec B_vec(2);
    SmartFieldVec C_vec(2);
    SmartFieldMat A_mat(2,2);

    SmartFieldMat M(2,2);
    M.setflag_inrealspace(true);
    M = m;

    Vec2dFieldType eye(2, Vec1dFieldType(2, FieldType(0.)));
    eye[0][0] = 1.; eye[0][1] = 0.;
    eye[1][0] = 0.; eye[1][1] = 1.;


    // }}}

    // --- Explicit derivatives --- 
    { // {{{

      std::cout << " * Explicit derivative operators\n";

      // scalar -> vector gradient operator
      A.setflag_inrealspace(true);
      A = a;

      A_vec.setflag_inrealspace(false);
      A.fft();
      _PhiOp->Del_f_ex(A, A_vec);
      A_vec.ifft();

      A_vec[0] -= da_dx;
      A_vec[1] -= da_dy;
      std::cout << std::setw(45) << std::left << "  - (Gradient) _PhiOp->Del_f_ex(scalar, vector)" << std::endl << "\t";
      std::cout << "x-derivative: " << A_vec[0].l2norm() << "\t";
      std::cout << "y-derivative: " << A_vec[1].l2norm() << std::endl;

      // vector -> tensor gradient operator
      A_vec.setflag_inrealspace(true);
      A_vec[0] = a;
      A_vec[1] = b;

      A_mat.setflag_inrealspace(false);
      A_vec.fft();
      _PhiOp->Del_f_ex(A_vec, A_mat);
      A_mat.ifft();

      A_mat(0,0) -= da_dx;
      A_mat(0,1) -= da_dy;
      A_mat(1,0) -= db_dx;
      A_mat(1,1) -= db_dy;
      std::cout << std::setw(45) << std::left << "  - (Gradient) _PhiOp->Del_f_ex(vector, tensor)" << std::endl << "\t";
      std::cout << "fn 1, x-derivative: " << A_mat(0,0).l2norm() << "\t";
      std::cout << "fn 1, y-derivative: " << A_mat(0,1).l2norm() << std::endl << "\t";
      std::cout << "fn 2, x-derivative: " << A_mat(1,0).l2norm() << "\t";
      std::cout << "fn 2, y-derivative: " << A_mat(1,1).l2norm() << std::endl;

      // divergence operator
      A_vec.setflag_inrealspace(true);
      A_vec[0] = a;
      A_vec[1] = b;
      A_vec.fft();
      A.setflag_inrealspace(false);
      _PhiOp->Div_f_ex(A_vec, A);
      A.ifft();

      A -= da_dx;
      A -= db_dy;
      std::cout << std::setw(45) << std::left << "  - (Divergence) _PhiOp->Div_f_ex(vector, scalar)" << "\t";
      std::cout << A.l2norm() << std::endl;

      // Laplacian operator
      A_vec.setflag_inrealspace(true);
      A_vec[0] = a;
      A_vec[1] = b;

      B_vec.setflag_inrealspace(false);
      A_vec.fft();
      _PhiOp->Del2_f_ex(A_vec, B_vec);
      B_vec.ifft();

      B_vec[0] -= d2a_dx2;
      B_vec[0] -= d2a_dy2;
      B_vec[1] -= d2b_dx2;
      B_vec[1] -= d2b_dy2;
      std::cout << std::setw(45) << std::left << "  - (Laplacian) _PhiOp->Del2_f_ex(vector, vector)" << std::endl << "\t";
      std::cout << "fn 1: " << B_vec[0].l2norm() << "\t";
      std::cout << "fn 2: " << B_vec[1].l2norm() << std::endl;

      // Biharmonic operator
      A_vec.setflag_inrealspace(true);
      A_vec[0] = a;
      A_vec[1] = b;

      B_vec.setflag_inrealspace(false);
      A_vec.fft();
      _PhiOp->Del4_f_ex(A_vec, B_vec);
      B_vec.ifft();

      B_vec[0] -= d4a_dx4;
      B_vec[0] -= d4a_dy4;
      B_vec[0] -= d4a_dx2dy2;
      B_vec[0] -= d4a_dx2dy2;
      B_vec[1] -= d4b_dx4;
      B_vec[1] -= d4b_dy4;
      B_vec[1] -= d4b_dx2dy2;
      B_vec[1] -= d4b_dx2dy2;
      std::cout << std::setw(45) << std::left << "  - (Biharmonic/Del^4) _PhiOp->Del4_f_ex(vector, vector)" << std::endl << "\t";
      std::cout << "fn 1: " << B_vec[0].l2norm() << "\t";
      std::cout << "fn 2: " << B_vec[1].l2norm() << std::endl;

      // del . (M . del f)
      A_vec.setflag_inrealspace(true);
      A_vec[0] = a;
      A_vec[1] = b;

      B_vec.setflag_inrealspace(false);
      A_vec.fft();
      _PhiOp->Del_A_Del_f_ex(A_vec, M, B_vec);
      B_vec.ifft();

      B_vec[0] -= del_M_del_ab_1;
      B_vec[1] -= del_M_del_ab_2;
      std::cout << std::setw(45) << std::left << "  - (Del . A . Del f) _PhiOp->Del_A_Del_f_ex(vector, vector)" << std::endl << "\t";
      std::cout << "fn 1: " << B_vec[0].l2norm() << "\t";
      std::cout << "fn 2: " << B_vec[1].l2norm() << std::endl;

      // del . (M . del del^{2} f)
      A_vec.setflag_inrealspace(true);
      A_vec[0] = a;
      A_vec[1] = b;

      B_vec.setflag_inrealspace(false);
      A_vec.fft();
      _PhiOp->Del_A_Del3_f_ex(A_vec, M, B_vec);
      B_vec.ifft();

      B_vec[0] -= del_M_del3_ab_1;
      B_vec[1] -= del_M_del3_ab_2;
      std::cout << std::setw(45) << std::left << "  - (Del . A . Del Del^2 f) _PhiOp->Del_A_Del3_f_ex(vector, vector)" << std::endl << "\t";
      std::cout << "fn 1: " << B_vec[0].l2norm() << "\t";
      std::cout << "fn 2: " << B_vec[1].l2norm() << std::endl;

    } // }}}

    // --- Implicit derivatives ---
    { // {{{

      std::cout << "\n";
      std::cout << " * Implicit derivative operators\n";
  
      SmartFieldOp FieldOp;
      SmartFieldOpMat FieldOpMat(2,2);
      SmartFieldOpMat Del2_f(2,2);

      // Identity operator
      A_vec.setflag_inrealspace(true);
      B_vec.setflag_inrealspace(true);
      A_vec[0] = a;
      A_vec[1] = b;
      A_vec.fft();

      FieldOpMat.zero();
      FieldOpMat.setflag_inrealspace( false );
      _PhiOp->Eye_im(FieldOpMat);
      FieldOpMat.OpDot( A_vec, B_vec );
      B_vec.ifft();

      B_vec[0] -= a;
      B_vec[1] -= b;
      std::cout << std::setw(45) << std::left << "  - (Identity) _PhiOp->Eye_im(vector)" << std::endl << "\t";
      std::cout << "fn 1: " << B_vec[0].l2norm() << "\t";
      std::cout << "fn 2: " << B_vec[1].l2norm() << std::endl;

      // Scalar Laplacian operator
      A.setflag_inrealspace(true);
      B.setflag_inrealspace(false);
      A = a;
      A.fft();

      FieldOp.zero();
      FieldOp.setflag_inrealspace( false );
      _PhiOp->Del2_im(FieldOp);
      FieldOp.OpDot( A, B );
      B.ifft();

      B -= d2a_dx2;
      B -= d2a_dy2;
      std::cout << std::setw(45) << std::left << "  - (Laplacian) _PhiOp->Del2_im(scalar)" << "\t";
      std::cout << B.l2norm() << "\n";

      // Laplacian operation on a vector
      A_vec.setflag_inrealspace(true);
      B_vec.setflag_inrealspace(false);
      A_vec[0] = a;
      A_vec[1] = b;
      A_vec.fft();

      FieldOpMat.zero();
      FieldOpMat.setflag_inrealspace( false );
      _PhiOp->A_Del2_im(m, FieldOpMat);
      FieldOpMat.OpDot( A_vec, B_vec );
      B_vec.ifft();

      B_vec[0] -= del_M_del_ab_1;
      B_vec[1] -= del_M_del_ab_2;
      std::cout << std::setw(45) << std::left << "  - (Laplacian on a vector) _PhiOp->A_Del2_im(vector)" << std::endl << "\t";
      std::cout << "fn 1: " << B_vec[0].l2norm() << "\t";
      std::cout << "fn 2: " << B_vec[1].l2norm() << std::endl;

      // Biharmonic operator, i.e. del^4
      A_vec.setflag_inrealspace(true);
      B_vec.setflag_inrealspace(false);
      A_vec[0] = a;
      A_vec[1] = b;
      A_vec.fft();

      FieldOpMat.zero();
      FieldOpMat.setflag_inrealspace( false );
      _PhiOp->A_Del4_im(m, FieldOpMat);
      FieldOpMat.OpDot( A_vec, B_vec );
      B_vec.ifft();

      B_vec[0] -= del_M_del3_ab_1;
      B_vec[1] -= del_M_del3_ab_2;
      std::cout << std::setw(45) << std::left << "  - (Biharmonic on a vector) _PhiOp->A_Del4_im(vector)" << std::endl << "\t";
      std::cout << "fn 1: " << B_vec[0].l2norm() << "\t";
      std::cout << "fn 2: " << B_vec[1].l2norm() << std::endl;

      // Nested Laplacian --> Biharmonic, i.e. del^2( Mij del^2 fj)
      A_vec.setflag_inrealspace(true);
      B_vec.setflag_inrealspace(false);
      A_vec[0] = a;
      A_vec[1] = b;
      A_vec.fft();

      Del2_f.zero();
      Del2_f.setflag_inrealspace( false );
      FieldOpMat.zero();
      FieldOpMat.setflag_inrealspace( false );
      _PhiOp->A_Del2_im(m, Del2_f);
      _PhiOp->A_Del2_F_im(eye, Del2_f, FieldOpMat);
      FieldOpMat.OpDot( A_vec, B_vec );
      B_vec.ifft();

      B_vec[0] -= del_M_del3_ab_1;
      B_vec[1] -= del_M_del3_ab_2;
      std::cout << std::setw(45) << std::left << "  - (Nested Laplacian: Del^2 [A . Del^2]) _PhiOp->A_Del2_F_im(vector)" << std::endl << "\t";
      std::cout << "fn 1: " << B_vec[0].l2norm() << "\t";
      std::cout << "fn 2: " << B_vec[1].l2norm() << std::endl;

    } // }}}

    // --- Test Implicit Solver ---
    { // {{{

      std::cout << "\n";
      std::cout << " * Implicit solver\n";

      // Total Solver (Del^2 B_vec = A_vec)
      SmartFieldOpMat T1(2, 2);

      A_vec.setflag_inrealspace(true);
      A_vec[0] = a;
      A_vec[1] = b;
      A_vec.fft();

      // T1 = I + del^2
      T1.zero();
      T1.setflag_inrealspace( false );
      _PhiOp->Eye_im(T1);
      _PhiOp->A_Del2_im(eye, T1);

      // solve del^2 B_vec = A_vec
      B_vec.zero();
      B_vec.setflag_inrealspace( false );
      _PhiOp->Solve_im_matrix(T1, A_vec, B_vec, 0);

      // Now do the reverse to test it
      C_vec.zero();
      C_vec.setflag_inrealspace( false );
      T1.OpDot(B_vec, C_vec);
      C_vec[0] -= A_vec[0];
      C_vec[1] -= A_vec[1];

      std::cout << std::setw(45) << std::left << "  - (Solver) _PhiOp->Sovler_im_matrix(Op, rhs_vec, x_vec)" << std::endl << "\t";
      std::cout << "fn 1: " << C_vec[0].l2norm() << "\t";
      std::cout << "fn 2: " << C_vec[1].l2norm() << std::endl;

      std::cout << "\n";
      std::cout << " * Timing tests of parts of the solver (no reported values)\n";

      // Sub-operations (for timing purposes, FD only)
      if (_CurrGrid->GridFlag() == 1) {

        // Convert SmartFields to Intercalated SmartFields
        UInt NIcomp = T1.Nrow();
        UInt BW = T1.BW();
        UInt Nx = T1.FDSize();
        UInt Nbands = 2*BW+1;

        MatInterFieldMat T1_int( Nbands, Nx, NIcomp, NIcomp );
        VecInterFieldVec rhs_int( Nx, NIcomp );
        VecInterFieldVec x_int( Nx, NIcomp );

        T1_int.setflag_inrealspace(false);
        rhs_int.setflag_inrealspace(false);
        x_int.setflag_inrealspace(false);

        rhs_int.smart_to_inter( A_vec );
        std::cout << std::setw(45) << std::left;
        std::cout << "  - (Smart->Inter Conversion) VecInterFieldVec.smart_to_inter(SmartFieldVec)" << std::endl;

        T1_int.smart_to_inter( T1 );
        std::cout << std::setw(45) << std::left;
        std::cout << "  - (Smart->Inter Conversion) MatInterFieldMat.smart_to_inter(SmartFieldMat)" << std::endl;

        // All BCs spend some time here calculating boundary conditions
        // I'm not going to reproduce this code

        // Doing the banded LU decomposition
        T1_int.band_LU();
        std::cout << std::setw(45) << std::left;
        std::cout << "  - (Banded LU Decomp) InterFieldMat.band_LU()" << std::endl;

        // Back Substitution to solve the system
        T1_int.band_solve( rhs_int, x_int );
        std::cout << std::setw(45) << std::left;
        std::cout << "  - (Banded Back Substitution) InterFieldMat.band_solve()" << std::endl;

        // Convert Intercalated SmartFields back to SmartFields
        x_int.inter_to_smart( B_vec );
        std::cout << std::setw(45) << std::left;
        std::cout << "  - (Inter->Smart Conversion) InterFieldMat.inter_to_smart(SmartFieldMat)" << std::endl;

      }
      else {
        std::cout << std::setw(45) << std::left << "  - (Parts of FD solver) Skipping ..." << std::endl;
      }
    } // }}}

  } // }}}
  std::cout << "\n ***Finished Differential Operators Test***\n";

} // }}}

//void time_tests( Grid* _CurrGrid )
//{ //{{{
//std::cout << "\n *** Smartfield Time Tests ***\n";    
//  { // {{{
//    
//    struct timeval start, stop;
//    double timetot;
//    int num_runs = 100000;
//    std::cout << "\n*** SmartField Time Tests ***\n\n";
//    std::cout << "Time results are from ";
//    std::cout << num_runs;
//    std::cout << "of runs\n\n";  
//    FieldType ftnum = 1.;
//  
//    FieldStack initstack;
//    Field<FieldType> myfield(*(initstack.CurrLayout), true, "");
//    myfield = FieldType(2.);
//  
//    std::vector< FieldType > myvec( SmartField::NPW );
//    for (UInt m =0; m<SmartField::NPW; m++)
//    {
//      myvec[m] = FieldType(3*(m+1));
//    }
//   
//    std::vector< std::vector<FieldType> > myvecvec( SmartField::FDSize, myvec );
//    for (UInt m =0; m<SmartField::FDSize; m++)
//    {
//      myvecvec[m] = std::vector<FieldType>( SmartField::NPW, FieldType(4*(m+1)) );
//    }
//  
//    std::vector< Field<FieldType> > myfieldvec( SmartField::FDSize, myfield );
//    for (UInt m =0; m<SmartField::FDSize; m++)
//    {
//      myfieldvec[m] = FieldType(5*(m+1));
//    }
//  
//    //assignment test
//    std::cout << "Assignment: \n";
//    {
//      SmartField tmp1, tmp2;
//      
//      gettimeofday(&start, NULL);
//      for (int i = 0; i <= num_runs; ++i) {
//        tmp1 = ftnum;
//      }
//      gettimeofday(&stop, NULL);
//      timetot = (stop.tv_sec - start.tv_sec) + 
//                  (stop.tv_usec - start.tv_usec)/1000000.0;
//      timetot = timetot/num_runs;
//      printf(" - Elapsed Time %10.3e (sec.)\n\n", timetot );
//    }
//    // Compound Subtraction test
//    std::cout << "Compound Subtraction: \n";
//    {
//      
//      SmartField tmp1, tmp2;
//      tmp1 = myfield;
//  
//      gettimeofday(&start, NULL);
//      for (int i = 0; i <= num_runs; ++i) {
//        tmp1 -= myfield;
//      }
//      gettimeofday(&stop, NULL);
//      timetot = (stop.tv_sec - start.tv_sec) + 
//                  (stop.tv_usec - start.tv_usec)/1000000.0;
//      timetot = timetot/num_runs;
//      printf(" - Elapsed Time %10.3e (sec.)\n\n", timetot );
//    }
//  
//    //Compound addition test
//    std::cout << "Compound Addition: \n";
//    {
//      SmartField tmp1, tmp2;
//      tmp1 = FieldType(0.);
//  
//      gettimeofday(&start, NULL);
//      for (int i = 0; i <= num_runs; ++i) {
//        tmp1 += myfieldvec;
//      }
//      gettimeofday(&stop, NULL);
//      timetot = (stop.tv_sec - start.tv_sec) + 
//                  (stop.tv_usec - start.tv_usec)/1000000.0;
//      timetot = timetot/num_runs;
//      printf(" - Elapsed Time %10.3e (sec.)\n\n", timetot );
//    }
//  
//    //Compound multiplication test
//    std::cout << "Compound multiplication: \n";
//    {
//      SmartField tmp1, tmp2;
//      tmp1 = FieldType(3.);
//      tmp2 = FieldType(3.);
//  
//      gettimeofday(&start, NULL);
//      for (int i = 0; i <= num_runs; ++i) {
//        tmp1 *= tmp2;
//      }
//      gettimeofday(&stop, NULL);
//      timetot = (stop.tv_sec - start.tv_sec) + 
//                  (stop.tv_usec - start.tv_usec)/1000000.0;
//      timetot = timetot/num_runs;
//      printf(" - Elapsed Time %10.3e (sec.)\n\n", timetot );
//    }
//  
//  
//    //Compound division test
//    std::cout << "Compound division: \n";
//    {
//      SmartField tmp1, tmp2;
//      tmp1 = FieldType(3.);
//      
//      gettimeofday(&start, NULL);
//      for (int i = 0; i <= num_runs; ++i) {
//        tmp1 /= ftnum;
//      }
//      gettimeofday(&stop, NULL);
//      timetot = (stop.tv_sec - start.tv_sec) + 
//                  (stop.tv_usec - start.tv_usec)/1000000.0;
//      timetot = timetot/num_runs;
//      printf(" - Elapsed Time %10.3e (sec.)\n\n", timetot );
//    }
//    //skipped file I/O and copying
//  
//    //Getting and setting stuff
//    std::cout << "Getters and setters and stuff:\n";
//    SmartField tmp1;
//    tmp1 = myfieldvec;
//    SmartField tmp2 = FieldType(3.0, 4.0);
//   
//    UInt idx1 = SmartField::NPW/3+1;
//  
//    for (int j = 1; j <= 15; ++j) {
//      std::cout << j;
//      gettimeofday(&start, NULL);
//      for (int i = 0; i <= num_runs; ++i) {
//        switch (j)
//        {
//          case 1: tmp1.l1norm();
//            break;
//          case 2: tmp1.max();
//            break;
//          case 3: tmp1.maxsigned();
//            break;
//          case 4: tmp1.minsigned();
//            break;
//          case 5: tmp1.maxidx();
//            break;
//          case 6: tmp1.sumelem();
//            break;
//          case 7: tmp1.integrate();
//            break;
//          case 8: tmp1.integrate_square();
//            break;
//          case 9: tmp1.log(tmp1);
//            break;
//          case 10: tmp1.exponentiate(tmp1);
//            break;
//          case 11: tmp1.sqrt(tmp1);
//            break;
//          case 12: tmp2.l1norm();
//            break;
//          case 13: tmp2.l2norm();
//            break;
//          case 14: tmp1.setelement(FieldType(1024.), idx1);
//            break;
//          case 15: tmp1.getelement( idx1 );
//            break;
//        }
//      }
//      gettimeofday(&stop, NULL);
//      timetot = (stop.tv_sec - start.tv_sec) + 
//                (stop.tv_usec - start.tv_usec)/1000000.0;
//      timetot = timetot/num_runs;
//      printf(" - Elapsed Time %10.3e (sec.)\n\n", timetot );
//    }
//  
//  //copy test
//    std::cout << "Copying test: \n";
//    gettimeofday(&start, NULL);
//    for (int i = 0; i <= num_runs/1000; ++i) {
//      tmp2.copytoarray( myvecvec );
//    }
//    gettimeofday(&stop, NULL);
//    timetot = (stop.tv_sec - start.tv_sec) + 
//              (stop.tv_usec - start.tv_usec)/1000000.0;
//    timetot = timetot/(num_runs/1000);
//    printf(" - Elapsed Time %10.3e (sec.)\n\n", timetot );
//  
//  } // }}}
//std::cout << "\n *** End Smartfield Time Tests ***\n";
//
//std::cout << "\n *** SmartfieldMat Time Tests ***\n";    
//  { // {{{
//    struct timeval start, stop;
//    double timetot;
//    int num_runs = 100000;
//  
//    std::cout << "Smartfield operation time tests\n\n";
//    std::cout << "Values averaged over " << num_runs << "runs\n\n";
//  
//  
//    SmartFieldMat A( 3, 3 );
//    SmartFieldMat B( 3, 3 );
//    SmartFieldMat C( 3, 3 );
//    SmartFieldVec b( 3 );
//    SmartFieldVec x( 3 );
//    FieldType ftnum(0.);
//    SmartField myfield(0.);
//    Vec2dFieldType Simple_A(3, Vec1dFieldType(3, 0.));
//    Vec1dFieldType Simple_x(3, 0.);
//  
//    A.zero();
//    B.zero();
//    C.zero();
//    b.zero();
//    x.zero();
//  
//    //assignment
//    std::cout << "Assignment test:\n";
//    {
//    
//      A[0][0] = FieldType(1.);
//      A[0][1] = FieldType(2.);
//      A[0][2] = FieldType(3.);
//      A[1][0] = FieldType(4.);
//      A[1][1] = FieldType(5.);
//      A[1][2] = FieldType(4.);
//      A[2][0] = FieldType(3.);
//      A[2][1] = FieldType(2.);
//      A[2][2] = FieldType(1.);
//  
//      B[0][0] = FieldType(2.);
//      B[0][1] = FieldType(5.);
//      B[0][2] = FieldType(-3.);
//      B[1][0] = FieldType(7.);
//      B[1][1] = FieldType(-4.);
//      B[1][2] = FieldType(8.);
//      B[2][0] = FieldType(0.);
//      B[2][1] = FieldType(-1.3);
//      B[2][2] = FieldType(-4);
//  
//      Simple_A[0][0] = FieldType(1.);
//      Simple_A[0][1] = FieldType(2.);
//      Simple_A[0][2] = FieldType(3.);
//      Simple_A[1][0] = FieldType(4.);
//      Simple_A[1][1] = FieldType(5.);
//      Simple_A[1][2] = FieldType(4.);
//      Simple_A[2][0] = FieldType(3.);
//      Simple_A[2][1] = FieldType(2.);
//      Simple_A[2][2] = FieldType(1.);
//  
//      myfield = FieldType(2.);
//      ftnum = 1.;
//    
//      gettimeofday(&start, NULL); 
//      for (int i = 0; i <= num_runs; ++i) {
//        A = B;
//      }
//      gettimeofday(&stop, NULL);
//      timetot = (stop.tv_sec - start.tv_sec) +
//                (stop.tv_usec - start.tv_usec)/1000000.0;
//      timetot = timetot/num_runs;
//      printf(" - Elapsed Time %10.3e (sec.)\n\n", timetot );
//    }
//  
//    //Addition
//    std::cout << "Addition Test:\n";
//    {
//  
//      B[0][0] = FieldType(2.);
//      B[0][1] = FieldType(5.);
//      B[0][2] = FieldType(-3.);
//      B[1][0] = FieldType(7.);
//      B[1][1] = FieldType(-4.);
//      B[1][2] = FieldType(8.);
//      B[2][0] = FieldType(0.);
//      B[2][1] = FieldType(-1.3);
//      B[2][2] = FieldType(-4);
//  
//      Simple_A[0][0] = FieldType(1.);
//      Simple_A[0][1] = FieldType(2.);
//      Simple_A[0][2] = FieldType(3.);
//      Simple_A[1][0] = FieldType(4.);
//      Simple_A[1][1] = FieldType(5.);
//      Simple_A[1][2] = FieldType(4.);
//      Simple_A[2][0] = FieldType(3.);
//      Simple_A[2][1] = FieldType(2.);
//      Simple_A[2][2] = FieldType(1.);
//  
//      myfield = FieldType(2.);
//      ftnum = 1.;
//  
//      A.zero();
//      
//      gettimeofday(&start, NULL); 
//      for (int i = 0; i <= num_runs; ++i) {
//        A += Simple_A;
//      }
//      gettimeofday(&stop, NULL);
//      timetot = (stop.tv_sec - start.tv_sec) +
//                (stop.tv_usec - start.tv_usec)/1000000.0;
//      timetot = timetot/num_runs;
//      printf(" - Elapsed Time %10.3e (sec.)\n\n", timetot );
//    }
//  
//    //Subtraction
//    std::cout << "Subtraction Test:\n";
//    {
//  
//      B[0][0] = FieldType(2.);
//      B[0][1] = FieldType(5.);
//      B[0][2] = FieldType(-3.);
//      B[1][0] = FieldType(7.);
//      B[1][1] = FieldType(-4.);
//      B[1][2] = FieldType(8.);
//      B[2][0] = FieldType(0.);
//      B[2][1] = FieldType(-1.3);
//      B[2][2] = FieldType(-4);
//  
//      Simple_A[0][0] = FieldType(1.);
//      Simple_A[0][1] = FieldType(2.);
//      Simple_A[0][2] = FieldType(3.);
//      Simple_A[1][0] = FieldType(4.);
//      Simple_A[1][1] = FieldType(5.);
//      Simple_A[1][2] = FieldType(4.);
//      Simple_A[2][0] = FieldType(3.);
//      Simple_A[2][1] = FieldType(2.);
//      Simple_A[2][2] = FieldType(1.);
//  
//      myfield = FieldType(2.);
//      ftnum = 1.;
//  
//      A.zero();
//      
//      gettimeofday(&start, NULL); 
//      for (int i = 0; i <= num_runs; ++i) {
//        A -= myfield;
//      }
//      gettimeofday(&stop, NULL);
//      timetot = (stop.tv_sec - start.tv_sec) +
//                (stop.tv_usec - start.tv_usec)/1000000.0;
//      timetot = timetot/num_runs;
//      printf(" - Elapsed Time %10.3e (sec.)\n\n", timetot );
//    }
//  
//    // multiplication
//    std::cout << "Multiplication test\n";
//    { 
//  
//      A[0][0] = FieldType(1.);
//      A[0][1] = FieldType(2.);
//      A[0][2] = FieldType(3.);
//      A[1][0] = FieldType(4.);
//      A[1][1] = FieldType(5.);
//      A[1][2] = FieldType(4.);
//      A[2][0] = FieldType(3.);
//      A[2][1] = FieldType(2.);
//      A[2][2] = FieldType(1.);
//  
//      B[0][0] = FieldType(2.);
//      B[0][1] = FieldType(5.);
//      B[0][2] = FieldType(-3.);
//      B[1][0] = FieldType(7.);
//      B[1][1] = FieldType(-4.);
//      B[1][2] = FieldType(8.);
//      B[2][0] = FieldType(0.);
//      B[2][1] = FieldType(-1.3);
//      B[2][2] = FieldType(-4);
//  
//      Simple_A[0][0] = FieldType(1.);
//      Simple_A[0][1] = FieldType(2.);
//      Simple_A[0][2] = FieldType(3.);
//      Simple_A[1][0] = FieldType(4.);
//      Simple_A[1][1] = FieldType(5.);
//      Simple_A[1][2] = FieldType(4.);
//      Simple_A[2][0] = FieldType(3.);
//      Simple_A[2][1] = FieldType(2.);
//      Simple_A[2][2] = FieldType(1.);
//  
//      myfield = FieldType(2.);
//      ftnum = 3.;
//  
//      A = Simple_A;
//  
//      gettimeofday(&start, NULL); 
//      for (int i = 0; i <= num_runs; ++i) {
//        A *= ftnum;
//      }
//      gettimeofday(&stop, NULL);
//      timetot = (stop.tv_sec - start.tv_sec) +
//                (stop.tv_usec - start.tv_usec)/1000000.0;
//      timetot = timetot/num_runs;
//      printf(" - Elapsed Time %10.3e (sec.)\n\n", timetot );
//    }
//  
//    // dot product 
//    std::cout << "Dot Product test\n";
//    { 
//  
//      A[0][0] = FieldType(1.);
//      A[0][1] = FieldType(2.);
//      A[0][2] = FieldType(3.);
//      A[1][0] = FieldType(4.);
//      A[1][1] = FieldType(5.);
//      A[1][2] = FieldType(4.);
//      A[2][0] = FieldType(3.);
//      A[2][1] = FieldType(2.);
//      A[2][2] = FieldType(1.);
//  
//      B[0][0] = FieldType(2.);
//      B[0][1] = FieldType(5.);
//      B[0][2] = FieldType(-3.);
//      B[1][0] = FieldType(7.);
//      B[1][1] = FieldType(-4.);
//      B[1][2] = FieldType(8.);
//      B[2][0] = FieldType(0.);
//      B[2][1] = FieldType(-1.3);
//      B[2][2] = FieldType(-4);
//  
//      x[0] = FieldType(1.);
//      x[1] = FieldType(2.);
//      x[2] = FieldType(3.);
//  
//      gettimeofday(&start, NULL); 
//      for (int i = 0; i <= num_runs; ++i) {
//        C.dot(A, B); // b[0] = 14; b[1] = 26; b[2] = 10;
//      }
//      gettimeofday(&stop, NULL);
//      timetot = (stop.tv_sec - start.tv_sec) +
//                (stop.tv_usec - start.tv_usec)/1000000.0;
//      timetot = timetot/num_runs;
//      printf(" - Elapsed Time %10.3e (sec.)\n\n", timetot );
//    }
//  
//    //Matrix operations test
//    std::cout << "Matrix Ops test\n";
//    {
//      bool realspace;
//  
//      for (int j = 1; j <= 9; ++j) {
//        A[0][0] = FieldType(1.);
//        A[0][1] = FieldType(2.);
//        A[0][2] = FieldType(3.);
//        A[1][0] = FieldType(4.);
//        A[1][1] = FieldType(5.);
//        A[1][2] = FieldType(4.);
//        A[2][0] = FieldType(3.);
//        A[2][1] = FieldType(2.);
//        A[2][2] = FieldType(1.);
//  
//        x[0] = FieldType(1.);
//        x[1] = FieldType(2.);
//        x[2] = FieldType(3.);
//  
//        b[0] = FieldType(14.);
//        b[1] = FieldType(26.);
//        b[2] = FieldType(10.);
//  
//        myfield.zero();
//  
//        std::cout << j;
//        gettimeofday(&start, NULL);
//        for (int i = 0; i <= num_runs; ++i) {
//          switch (j)
//          {
//            case 1: A.norm(myfield);
//              break;
//            case 2: B.outer(x, b);
//              break;
//            case 3: B.transpose();
//              break;
//            case 4: x.linsolve(A, b); 
//              break;
//            case 5: B.invert( A );
//              break;
//            case 6: realspace = A.getflag_inrealspace();
//              break;
//            case 7: A.setflag_inrealspace(true);
//              break;
//            case 8: A.getelement(1);
//              break;
//            case 9: A.fft();
//                    A.ifft();
//                  //  A.setflag_inrealspace(true);
//              break;
//     //       case 10: A.ifft();
//       //              A.setflag_inrealspace(false);
//         //     break;
//          }
//        }
//      gettimeofday(&stop, NULL);
//      timetot = (stop.tv_sec - start.tv_sec) +
//                (stop.tv_usec - start.tv_usec)/1000000.0;
//      timetot = timetot/num_runs;
//      printf(" - Elapsed Time %10.3e (sec.)\n\n", timetot );
//  
//      }
//    }
//  
//    // subscript and slice (only do if FDSize > 1)
//    std::cout << " * Subscript and Slice (for Finite Differences if FDSize > 1)\n";
//    std::cout << "  - FDSize = " << _CurrGrid->FDSize() << std::endl;
//    if (_CurrGrid->FDSize() > 1)
//    {
//  
//      A.zero();
//      A[0][0][1] = FieldType(1.);
//      A[0][1][1] = FieldType(2.);
//      A[0][2][1] = FieldType(3.);
//      A[1][0][1] = FieldType(4.);
//      A[1][1][1] = FieldType(5.);
//      A[1][2][1] = FieldType(4.);
//      A[2][0][1] = FieldType(3.);
//      A[2][1][1] = FieldType(2.);
//      A[2][2][1] = FieldType(1.);
//  
//      A[0][0][3] = FieldType(1.);
//      A[0][1][3] = FieldType(2.);
//      A[0][2][3] = FieldType(3.);
//      A[1][0][3] = FieldType(4.);
//      A[1][1][3] = FieldType(5.);
//      A[1][2][3] = FieldType(4.);
//      A[2][0][3] = FieldType(3.);
//      A[2][1][3] = FieldType(2.);
//      A[2][2][3] = FieldType(1.);
//  
//  
//      B.zero();
//  
//      int Nx = _CurrGrid->FDSize(); // need an int to do modulus
//      std::vector<UInt> idx(Nx,0);
//      for (int i=0; i<Nx; i++)
//      {
//        idx[i] = ((i+1+Nx)%Nx);
//      }
//      
//      std::cout << "subscript\n";
//      gettimeofday(&start, NULL); 
//      for (int i = 0; i <= num_runs; ++i) {
//        B.subscript(A, idx);
//      }
//      gettimeofday(&stop, NULL);
//      timetot = (stop.tv_sec - start.tv_sec) +
//                (stop.tv_usec - start.tv_usec)/1000000.0;
//      timetot = timetot/num_runs;
//      printf(" - Elapsed Time %10.3e (sec.)\n\n", timetot );
//        
//      std::vector<UInt> idx2(3, 0);
//      idx2[0] = 1;
//      idx2[1] = 3;
//      idx2[2] = 5;
//        
//      B.zero();
//      
//      std::cout << "slice\n";
//      gettimeofday(&start, NULL); 
//      for (int i = 0; i <= num_runs; ++i) {
//        B.slice(A, idx);
//      }
//      gettimeofday(&stop, NULL);
//      timetot = (stop.tv_sec - start.tv_sec) +
//                (stop.tv_usec - start.tv_usec)/1000000.0;
//      timetot = timetot/num_runs;
//      printf(" - Elapsed Time %10.3e (sec.)\n\n", timetot );
//  
//    }
//  
//  } // }}}
//std::cout << "\n *** End SmartfieldiMat Time Tests ***\n";
//
//std::cout << "\n *** Differential Operators Time Tests ***\n";
//  { //{{{
//    struct timeval start, stop;
//    double timetot;
//    int num_runs = 100000;
//
//    { // {{{
//
//      // --- set up operators (model B only) ---
//      // {{{
//      int _Ncomp, _NIcomp;
//      int _BCFlag; // Which BCs are we using?
//      BCs_Factory         *_BCFactory; // for making operators
//      Operators_Factory   *_OpFactory; // for making operators
//      Operators_Base *_PhiOp; // Operators for Phi fields
//      BCs_Base *_PhiBCs; // BCs object
//
//      // set up the BCs and operators for phi
//      std::string filename( "params_BCs_phi.in" );
//
//      // Read in BC key from file
//      std::string tmp_str;
//      std::ifstream f2(filename.c_str());
//      if ( f2.is_open() )
//      {
//        rapidjson::IStreamWrapper isw(f2);
//        rapidjson::Document d;
//        d.ParseStream(isw);
//        _BCFlag = d["BCFlag"].GetInt();
//      }
//      else
//      {
//        std::cout << "*** Error in TimeInt_Base :: setup_operators() ***\n";
//        std::cout << "*** Cannot open " << filename << " ***\n";
//        exit(1);
//      }
//      
//      // call factory to create new BC object based on the kind of BCs
//      _PhiBCs=_BCFactory->make_BCs( _BCFlag, 0, filename, _NIcomp, _CurrGrid );
//
//      // call factory to make new operators object for phi
//      _PhiOp = _OpFactory->make_operators( _CurrGrid, _PhiBCs );
//      // }}}
//
//      // --- Create a test field and a bunch of derivatives ---
//      // {{{
//      RealType x, y;
//
//      int FDSize = _CurrGrid->FDSize();
//      int NPW = _CurrGrid->NPW();
//      double Lx = _CurrGrid->Lx();
//      double Ly = _CurrGrid->Ly();
//      double dx = _CurrGrid->Dx();
//      double dy = _CurrGrid->Ly()/_CurrGrid->Ny();
//
//      Vec2dFieldType a(FDSize, Vec1dFieldType(NPW, FieldType(0.)));
//      Vec2dFieldType da_dx(FDSize, Vec1dFieldType(NPW, FieldType(0.)));
//      Vec2dFieldType da_dy(FDSize, Vec1dFieldType(NPW, FieldType(0.)));
//      Vec2dFieldType d2a_dxdy(FDSize, Vec1dFieldType(NPW, FieldType(0.)));
//      Vec2dFieldType d2a_dx2(FDSize, Vec1dFieldType(NPW, FieldType(0.)));
//      Vec2dFieldType d2a_dy2(FDSize, Vec1dFieldType(NPW, FieldType(0.)));
//      Vec2dFieldType d4a_dx4(FDSize, Vec1dFieldType(NPW, FieldType(0.)));
//      Vec2dFieldType d4a_dy4(FDSize, Vec1dFieldType(NPW, FieldType(0.)));
//      Vec2dFieldType d4a_dx2dy2(FDSize, Vec1dFieldType(NPW, FieldType(0.)));
//
//      Vec2dFieldType b(FDSize, Vec1dFieldType(NPW, FieldType(0.)));
//      Vec2dFieldType db_dx(FDSize, Vec1dFieldType(NPW, FieldType(0.)));
//      Vec2dFieldType db_dy(FDSize, Vec1dFieldType(NPW, FieldType(0.)));
//      Vec2dFieldType d2b_dxdy(FDSize, Vec1dFieldType(NPW, FieldType(0.)));
//      Vec2dFieldType d2b_dx2(FDSize, Vec1dFieldType(NPW, FieldType(0.)));
//      Vec2dFieldType d2b_dy2(FDSize, Vec1dFieldType(NPW, FieldType(0.)));
//      Vec2dFieldType d4b_dx4(FDSize, Vec1dFieldType(NPW, FieldType(0.)));
//      Vec2dFieldType d4b_dy4(FDSize, Vec1dFieldType(NPW, FieldType(0.)));
//      Vec2dFieldType d4b_dx2dy2(FDSize, Vec1dFieldType(NPW, FieldType(0.)));
//
//      Vec2dFieldType del_M_del_ab_1(FDSize, Vec1dFieldType(NPW, FieldType(0.)));
//      Vec2dFieldType del_M_del_ab_2(FDSize, Vec1dFieldType(NPW, FieldType(0.)));
//      Vec2dFieldType del_M_del3_ab_1(FDSize, Vec1dFieldType(NPW, FieldType(0.)));
//      Vec2dFieldType del_M_del3_ab_2(FDSize, Vec1dFieldType(NPW, FieldType(0.)));
//      Vec2dFieldType m(2, Vec1dFieldType(2, FieldType(0.)));
//      m[0][0] = 0.7;
//      m[0][1] = -2;
//      m[1][0] = 5;
//      m[1][1] = -0.33;
//      Vec2dFieldType eye(2, Vec1dFieldType(2, FieldType(0.)));
//      eye[0][0] = 1.;
//      eye[0][1] = 0.;
//      eye[1][0] = 0.;
//      eye[1][1] = 1.;
//
//      for( int i = 0; i < FDSize; i++) {
//        for( int j = 0; j < NPW; j++) {
//      
//          x = dx*(i + floor(j/_CurrGrid->Ny()));
//          y = dy*(j%_CurrGrid->Ny());
//    
//          a[i][j] = sin(2*PI*x/Lx)*sin(2*PI*y/Ly);
//          da_dx[i][j] = 2*PI/Lx*cos(2*PI*x/Lx)*sin(2*PI*y/Ly);
//          da_dy[i][j] = 2*PI/Ly*sin(2*PI*x/Lx)*cos(2*PI*y/Ly);
//          d2a_dx2[i][j] = -4*PI*PI/Lx/Lx*sin(2*PI*x/Lx)*sin(2*PI*y/Ly);
//          d2a_dy2[i][j] = -4*PI*PI/Ly/Ly*sin(2*PI*x/Lx)*sin(2*PI*y/Ly);
//          d4a_dx4[i][j] = pow(2*PI/Lx, 4)*sin(2*PI*x/Lx)*sin(2*PI*y/Ly);
//          d4a_dy4[i][j] = pow(2*PI/Ly, 4)*sin(2*PI*x/Lx)*sin(2*PI*y/Ly);
//          d4a_dx2dy2[i][j] = pow(2*PI/Lx, 2)*pow(2*PI/Ly, 2)*sin(2*PI*x/Lx)*sin(2*PI*y/Ly);
//
//          b[i][j] = cos(2*PI*x/Lx)*cos(2*PI*y/Ly);
//          db_dx[i][j] = -2*PI/Lx*sin(2*PI*x/Lx)*cos(2*PI*y/Ly);
//          db_dy[i][j] = -2*PI/Ly*cos(2*PI*x/Lx)*sin(2*PI*y/Ly);
//          d2b_dx2[i][j] = -4*PI*PI/Lx/Lx*cos(2*PI*x/Lx)*cos(2*PI*y/Ly);
//          d2b_dy2[i][j] = -4*PI*PI/Ly/Ly*cos(2*PI*x/Lx)*cos(2*PI*y/Ly);
//          d4b_dx4[i][j] = pow(2*PI/Lx, 4)*cos(2*PI*x/Lx)*cos(2*PI*y/Ly);
//          d4b_dy4[i][j] = pow(2*PI/Ly, 4)*cos(2*PI*x/Lx)*cos(2*PI*y/Ly);
//          d4b_dx2dy2[i][j] = pow(2*PI/Lx, 2)*pow(2*PI/Ly, 2)*cos(2*PI*x/Lx)*cos(2*PI*y/Ly);
//
//          del_M_del_ab_1[i][j] = m[0][0]*(d2a_dx2[i][j] + d2a_dy2[i][j]) + m[0][1]*(d2b_dx2[i][j] + d2b_dy2[i][j]);
//          del_M_del_ab_2[i][j] = m[1][0]*(d2a_dx2[i][j] + d2a_dy2[i][j]) + m[1][1]*(d2b_dx2[i][j] + d2b_dy2[i][j]);
//          del_M_del3_ab_1[i][j] = m[0][0]*(d4a_dx4[i][j] + 2.*d4a_dx2dy2[i][j] + d4a_dy4[i][j])
//                                + m[0][1]*(d4b_dx4[i][j] + 2.*d4b_dx2dy2[i][j] + d4b_dy4[i][j]);
//          del_M_del3_ab_2[i][j] = m[1][0]*(d4a_dx4[i][j] + 2.*d4a_dx2dy2[i][j] + d4a_dy4[i][j]) 
//                                + m[1][1]*(d4b_dx4[i][j] + 2.*d4b_dx2dy2[i][j] + d4b_dy4[i][j]);
//
//        }
//      }
//      
//      SmartField A;
//      SmartField B;
//      SmartFieldVec A_vec(2);
//      SmartFieldVec B_vec(2);
//      SmartFieldVec C_vec(2);
//      SmartFieldMat A_mat(2,2);
//      SmartFieldMat M(2,2);
//      M = m;
//      SmartFieldVec tmpVec(2);
//      SmartFieldMat tmpMat(2,2);
//      SmartField tmp;
//      // }}}
//
//      // --- Explicit derivatives --- 
//      { // {{{
//
//        std::cout << " * Explicit derivative operators\n";
//
//        // scalar -> vector gradient operator
//        std::cout << "Vector gradient operator\n";
//        A.setflag_inrealspace(true);
//        A = a;
//
//        A_vec.setflag_inrealspace(false);
//        A.fft();
//        tmpVec = A_vec;
//
//        gettimeofday(&start, NULL);
//        for (int i = 0; i <= num_runs; ++i) {
//          _PhiOp->Del_f_ex(A, A_vec);
//          A_vec = tmpVec;
//        }
//        gettimeofday(&stop, NULL);
//        timetot = (stop.tv_sec - start.tv_sec) + 
//                    (stop.tv_usec - start.tv_usec)/1000000.0;
//        timetot = timetot/num_runs;
//        printf(" - Elapsed Time %10.3e (sec.)\n\n", timetot );
//        A_vec.ifft();
//
//        // vector -> tensor gradient operator
//        std::cout << "tensor gradient operator\n";
//        A_vec.setflag_inrealspace(true);
//        A_vec[0] = a;
//        A_vec[1] = b;
//
//        A_mat.setflag_inrealspace(false);
//        A_vec.fft();
//        tmpMat = A_mat;
//
//        gettimeofday(&start, NULL);
//        for (int i = 0; i <= num_runs; ++i) {
//          _PhiOp->Del_f_ex(A_vec, A_mat);
//          A_mat = tmpMat;
//        }
//        gettimeofday(&stop, NULL);
//        timetot = (stop.tv_sec - start.tv_sec) + 
//                    (stop.tv_usec - start.tv_usec)/1000000.0;
//        timetot = timetot/num_runs;
//        printf(" - Elapsed Time %10.3e (sec.)\n\n", timetot );
//        A_mat.ifft();
//
//        // divergence operator
//        std::cout << "divergence operator\n";
//        A_vec.setflag_inrealspace(true);
//        A_vec[0] = a;
//        A_vec[1] = b;
//        A_vec.fft();
//        A.setflag_inrealspace(false);
//        tmp = A;
//        gettimeofday(&start, NULL);
//        for (int i = 0; i <= num_runs; ++i) {
//          _PhiOp->Div_f_ex(A_vec, A);
//          A = tmp;
//        }
//        gettimeofday(&stop, NULL);
//        timetot = (stop.tv_sec - start.tv_sec) + 
//                    (stop.tv_usec - start.tv_usec)/1000000.0;
//        timetot = timetot/num_runs;
//        printf(" - Elapsed Time %10.3e (sec.)\n\n", timetot );
//        A.ifft();
//
//        // Laplacian operator
//        std::cout << "Laplacian operator\n";
//        A_vec.setflag_inrealspace(true);
//        A_vec[0] = a;
//        A_vec[1] = b;
//
//        B_vec.setflag_inrealspace(false);
//        A_vec.fft();
//        tmpVec = B_vec;
//        gettimeofday(&start, NULL);
//        for (int i = 0; i <= num_runs; ++i) {
//          _PhiOp->Del2_f_ex(A_vec, B_vec);
//          B_vec = tmpVec;
//        }
//        gettimeofday(&stop, NULL);
//        timetot = (stop.tv_sec - start.tv_sec) + 
//                    (stop.tv_usec - start.tv_usec)/1000000.0;
//        timetot = timetot/num_runs;
//        printf(" - Elapsed Time %10.3e (sec.)\n\n", timetot );
//        B_vec.ifft();
//
//        // Biharmonic operator
//        std::cout << "Biharmonic operator\n";
//        A_vec.setflag_inrealspace(true);
//        A_vec[0] = a;
//        A_vec[1] = b;
//
//        B_vec.setflag_inrealspace(false);
//        A_vec.fft();
//        gettimeofday(&start, NULL);
//        for (int i = 0; i <= num_runs; ++i) {
//          _PhiOp->Del4_f_ex(A_vec, B_vec);
//          B_vec = tmpVec;
//        }
//        gettimeofday(&stop, NULL);
//        timetot = (stop.tv_sec - start.tv_sec) + 
//                    (stop.tv_usec - start.tv_usec)/1000000.0;
//        timetot = timetot/num_runs;
//        printf(" - Elapsed Time %10.3e (sec.)\n\n", timetot );
//        B_vec.ifft();
//
//        // del . (M . del f)
//        std::cout << "del . (M. del f)\n";
//        A_vec.setflag_inrealspace(true);
//        A_vec[0] = a;
//        A_vec[1] = b;
//
//        B_vec.setflag_inrealspace(false);
//        A_vec.fft();
//        gettimeofday(&start, NULL);
//        for (int i = 0; i <= num_runs; ++i) {
//          _PhiOp->Del_A_Del_f_ex(A_vec, M, B_vec);
//          B_vec = tmpVec;
//        }
//        gettimeofday(&stop, NULL);
//        timetot = (stop.tv_sec - start.tv_sec) + 
//                    (stop.tv_usec - start.tv_usec)/1000000.0;
//        timetot = timetot/num_runs;
//        printf(" - Elapsed Time %10.3e (sec.)\n\n", timetot );
//        B_vec.ifft();
//
//        // del . (M . del del^{2} f)
//        std::cout << "del . (M . del del^{2} f)\n";
//        A_vec.setflag_inrealspace(true);
//        A_vec[0] = a;
//        A_vec[1] = b;
//
//        B_vec.setflag_inrealspace(false);
//        A_vec.fft();
//        gettimeofday(&start, NULL);
//        for (int i = 0; i <= num_runs; ++i) {
//          _PhiOp->Del_A_Del3_f_ex(A_vec, M, B_vec);
//          B_vec = tmpVec;
//        }
//        gettimeofday(&stop, NULL);
//        timetot = (stop.tv_sec - start.tv_sec) + 
//                    (stop.tv_usec - start.tv_usec)/1000000.0;
//        timetot = timetot/num_runs;
//        printf(" - Elapsed Time %10.3e (sec.)\n\n", timetot );
//        B_vec.ifft();
//
//      } // }}}
//
//      // --- Implicit derivatives ---
//      { // {{{
//  
//        std::cout << "\n";
//        std::cout << " * Implicit derivative operators\n";
//    
//        SmartFieldOp FieldOp;
//        SmartFieldOpMat FieldOpMat(2, 2);
//        SmartFieldOpMat Del2_f(2, 2);
//    
//        // Identity operator
//        std::cout << "Identity operator\n";
//        A_vec.setflag_inrealspace(true);
//        B_vec.setflag_inrealspace(true);
//        A_vec[0]= a;
//        A_vec[1]= b;
//        A_vec.fft();
//  
//        FieldOpMat.zero();
//        FieldOpMat.setflag_inrealspace( false );
//        gettimeofday(&start, NULL);
//        for (int i = 0; i <= num_runs; ++i) {
//          _PhiOp->Eye_im(FieldOpMat);
//          FieldOpMat.OpDot( A_vec, B_vec );
//        }
//        gettimeofday(&stop, NULL);
//        timetot = (stop.tv_sec - start.tv_sec) + 
//                    (stop.tv_usec - start.tv_usec)/1000000.0;
//        timetot = timetot/num_runs;
//        printf(" - Elapsed Time %10.3e (sec.)\n\n", timetot );
//        B_vec.ifft();
//  
//        // Scalar Laplacian operator
//        std::cout << "Scalar Laplacian operator\n";
//        A.setflag_inrealspace(true);
//        B.setflag_inrealspace(false);
//        A = a;
//        A.fft();
//  
//        FieldOp.zero();
//        FieldOp.setflag_inrealspace( false );
//        gettimeofday(&start, NULL);
//        for (int i = 0; i <= num_runs; ++i) {
//          _PhiOp->Del2_im(FieldOp);
//          FieldOp.OpDot( A, B );
//        }
//        gettimeofday(&stop, NULL);
//        timetot = (stop.tv_sec - start.tv_sec) + 
//                    (stop.tv_usec - start.tv_usec)/1000000.0;
//        timetot = timetot/num_runs;
//        printf(" - Elapsed Time %10.3e (sec.)\n\n", timetot );
//        B.ifft();
//  
//        // Laplacian operation on a vector
//        std::cout << "Laplacian operation on a vector\n";
//        A_vec.setflag_inrealspace(true);
//        B_vec.setflag_inrealspace(false);
//        A_vec[0] = a;
//        A_vec[1] = b;
//        A_vec.fft();
//  
//        FieldOpMat.zero();
//        FieldOpMat.setflag_inrealspace( false );
//        gettimeofday(&start, NULL);
//        for (int i = 0; i <= num_runs; ++i) {
//          _PhiOp->A_Del2_im(m, FieldOpMat);
//          FieldOpMat.OpDot( A_vec, B_vec );
//        }
//        gettimeofday(&stop, NULL);
//        timetot = (stop.tv_sec - start.tv_sec) + 
//                    (stop.tv_usec - start.tv_usec)/1000000.0;
//        timetot = timetot/num_runs;
//        printf(" - Elapsed Time %10.3e (sec.)\n\n", timetot );
//        B_vec.ifft();
//  
//        // Biharmonic operator, i.e. del^4
//        std::cout << "Biharmonic operator\n";
//        A_vec.setflag_inrealspace(true);
//        B_vec.setflag_inrealspace(false);
//        A_vec[0] = a;
//        A_vec[1] = b;
//        A_vec.fft();
//  
//        FieldOpMat.zero();
//        FieldOpMat.setflag_inrealspace( false );
//        gettimeofday(&start, NULL);
//        for (int i = 0; i <= num_runs; ++i) {
//          _PhiOp->A_Del4_im(m, FieldOpMat);
//          FieldOpMat.OpDot( A_vec, B_vec );
//        }
//        gettimeofday(&stop, NULL);
//        timetot = (stop.tv_sec - start.tv_sec) + 
//                    (stop.tv_usec - start.tv_usec)/1000000.0;
//        timetot = timetot/num_runs;
//        printf(" - Elapsed Time %10.3e (sec.)\n\n", timetot );
//        B_vec.ifft();
//  
//        // Nested Laplacian --> Biharmonic, i.e. del^2( Mij del^2 fj)
//        std::cout << "Nested Laplacian -->Biharmonic\n";
//        A_vec.setflag_inrealspace(true);
//        B_vec.setflag_inrealspace(false);
//        A_vec[0] = a;
//        A_vec[1] = b;
//        A_vec.fft();
//  
//        Del2_f.zero();
//        Del2_f.setflag_inrealspace( false );
//        FieldOpMat.zero();
//        FieldOpMat.setflag_inrealspace( false );
//        gettimeofday(&start, NULL);
//        for (int i = 0; i <= num_runs; ++i) {
//          _PhiOp->A_Del2_im(m, Del2_f);
//          _PhiOp->A_Del2_F_im(eye, Del2_f, FieldOpMat);
//          FieldOpMat.OpDot( A_vec, B_vec );
//        }
//        gettimeofday(&stop, NULL);
//        timetot = (stop.tv_sec - start.tv_sec) + 
//                    (stop.tv_usec - start.tv_usec)/1000000.0;
//        timetot = timetot/num_runs;
//        printf(" - Elapsed Time %10.3e (sec.)\n\n", timetot );
//        B_vec.ifft();
//  
//  
//      } // }}}
//
//      // --- Test Implicit Solver ---
//      { // {{{
//  
//        std::cout << "\n";
//        std::cout << " * Implicit solver\n";
//  
//        // Total Solver (Del^2 B_vec = A_vec)
//        SmartFieldOpMat T1(2, 2);
//  
//        A_vec.setflag_inrealspace(true);
//        A_vec[0] = a;
//        A_vec[1] = b;
//        A_vec.fft();
//  
//        // T1 = del^2
//        T1.zero();
//        T1.setflag_inrealspace( false );
//        _PhiOp->Eye_im(T1);
//        _PhiOp->A_Del2_im(eye, T1);
//  
//        // solve del^2 B_vec = A_vec
//        B_vec.zero();
//        B_vec.setflag_inrealspace( false );
//        gettimeofday(&start, NULL);
//        for (int i = 0; i <= num_runs; ++i) {
//          _PhiOp->Solve_im_matrix(T1, A_vec, B_vec, 0);
//        }
//        gettimeofday(&stop, NULL);
//        timetot = (stop.tv_sec - start.tv_sec) + 
//                    (stop.tv_usec - start.tv_usec)/1000000.0;
//        timetot = timetot/num_runs;
//        printf(" - Elapsed Time %10.3e (sec.)\n\n", timetot );
//  
//        // Now do the reverse to test it
//        C_vec.zero();
//        C_vec.setflag_inrealspace( false );
//        T1.OpDot(B_vec, C_vec);
//        C_vec[0] -= A_vec[0];
//        C_vec[1] -= A_vec[1];
//  
//        std::cout << std::setw(45) << std::left << "  - (Solver) _PhiOp->Solve_im_matrix(Op, rhs_vec, x_vec)" << std::endl << "\t";
//  
//        std::cout << "\n";
//        std::cout << " * Timing tests of parts of the solver (no reported values)\n";
//  
//        // Sub-operations (for timing purposes, FD only)
//        if (_CurrGrid->GridFlag() == 1) {
//  
//          // Convert SmartFields to Intercalated SmartFields
//          UInt NIcomp = T1.Nrow();
//          UInt BW = T1.BW();
//          UInt Nx = T1.FDSize();
//          UInt Nbands = 2*BW+1;
//  
//          MatInterFieldMat T1_int( Nbands, Nx, NIcomp, NIcomp );
//          VecInterFieldVec rhs_int( Nx, NIcomp );
//          VecInterFieldVec x_int( Nx, NIcomp );
//  
//          T1_int.setflag_inrealspace(false);
//          rhs_int.setflag_inrealspace(false);
//          x_int.setflag_inrealspace(false);
//  
//          std::cout << std::setw(45) << std::left;
//          std::cout << "  - (Smart->Inter Conversion) VecInterFieldVec.smart_to_inter(SmartFieldVec)" << std::endl;
//          gettimeofday(&start, NULL);
//          for (int i = 0; i <= num_runs; ++i) {
//            rhs_int.smart_to_inter( A_vec );
//          }
//          gettimeofday(&stop, NULL);
//          timetot = (stop.tv_sec - start.tv_sec) + 
//                      (stop.tv_usec - start.tv_usec)/1000000.0;
//          timetot = timetot/num_runs;
//          printf(" - Elapsed Time %10.3e (sec.)\n\n", timetot );
//  
//          std::cout << std::setw(45) << std::left;
//          std::cout << "  - (Smart->Inter Conversion) MatInterFieldMat.smart_to_inter(SmartFieldMat)" << std::endl;
//          gettimeofday(&start, NULL);
//          for (int i = 0; i <= num_runs; ++i) {
//            T1_int.smart_to_inter( T1 );
//          }
//          gettimeofday(&stop, NULL);
//          timetot = (stop.tv_sec - start.tv_sec) + 
//                      (stop.tv_usec - start.tv_usec)/1000000.0;
//          timetot = timetot/num_runs;
//          printf(" - Elapsed Time %10.3e (sec.)\n\n", timetot );
//  
//          // All BCs spend some time here calculating boundary conditions
//          // I'm not going to reproduce this code
//  
//          // Doing the banded LU decomposition
//          std::cout << std::setw(45) << std::left;
//          std::cout << "  - (Banded LU Decomp) InterFieldMat.band_LU()" << std::endl;
//          gettimeofday(&start, NULL);
//          for (int i = 0; i <= num_runs; ++i) {
//            T1_int.band_LU();
//          }
//          gettimeofday(&stop, NULL);
//          timetot = (stop.tv_sec - start.tv_sec) + 
//                      (stop.tv_usec - start.tv_usec)/1000000.0;
//          timetot = timetot/num_runs;
//          printf(" - Elapsed Time %10.3e (sec.)\n\n", timetot );
//  
//          // Back Substitution to solve the system
//          std::cout << std::setw(45) << std::left;
//          std::cout << "  - (Banded Back Substitution) InterFieldMat.band_solve()" << std::endl;
//          gettimeofday(&start, NULL);
//          for (int i = 0; i <= num_runs; ++i) {
//            T1_int.band_solve( rhs_int, x_int );
//          }
//          gettimeofday(&stop, NULL);
//          timetot = (stop.tv_sec - start.tv_sec) + 
//                      (stop.tv_usec - start.tv_usec)/1000000.0;
//          timetot = timetot/num_runs;
//          printf(" - Elapsed Time %10.3e (sec.)\n\n", timetot );
//  
//          // Convert Intercalated SmartFields back to SmartFields
//          std::cout << std::setw(45) << std::left;
//          std::cout << "  - (Inter->Smart Conversion) InterFieldMat.inter_to_smart(SmartFieldMat)" << std::endl;
//          gettimeofday(&start, NULL);
//          for (int i = 0; i <= num_runs; ++i) {
//            x_int.inter_to_smart( B_vec );
//          }
//          gettimeofday(&stop, NULL);
//          timetot = (stop.tv_sec - start.tv_sec) + 
//                      (stop.tv_usec - start.tv_usec)/1000000.0;
//          timetot = timetot/num_runs;
//          printf(" - Elapsed Time %10.3e (sec.)\n\n", timetot );
//  
//        }
//        else {
//          std::cout << std::setw(45) << std::left << "  - (Parts of FD solver) Skipping ..." << std::endl;
//        }
//      } // }}}
//  
//    } //}}}
//
//  } //}}}
//} //}}}

