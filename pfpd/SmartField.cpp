/***************************************************************
*
* Class defining a wrapper around the field class for better
* use of scratch space
* 
* DRT -- Wed, 19 Aug 2015
*
****************************************************************/

#include "SmartField.h"

UInt SmartField::ICFlag = 0;
UInt SmartField::BinaryInFlag = 0;
UInt SmartField::BinaryOutFlag = 0;
UInt SmartField::FDSize = 0;
UInt SmartField::NPW = 0;

// ----------------- Non-member routines for the statics -----------------

void init_SmartField_statics( UInt fd_size )
{ // {{{

  // Read params_ICs.in
  std::string tmp_str;

  if (jsonFile("params_ICs.in"))
  {
    std::ifstream f1("params_ICs.in");
    if ( f1.is_open() )
    {
      rapidjson::IStreamWrapper isw(f1);
      rapidjson::Document d;
      d.ParseStream(isw);
      SmartField::BinaryInFlag = d["BinaryIn_flag"].GetInt();

      SmartField::BinaryOutFlag = d["BinaryOut_flag"].GetInt();

      SmartField::ICFlag = d["IC_flag"].GetInt();
    }
    else
    {
      std::cout << "cannot open params_ICs.in" << std::endl;
      exit(1);
    }
  }
  else
  {//legacy
    std::ifstream f1("params_ICs.in");
    if ( f1.is_open() )
    {
      std::getline(f1, tmp_str, '#');
      SmartField::BinaryInFlag = atoi(tmp_str.c_str());
      std::getline(f1, tmp_str, '\n');

      std::getline(f1, tmp_str, '#');
      SmartField::BinaryOutFlag = atoi(tmp_str.c_str());
      std::getline(f1, tmp_str, '\n');

      std::getline(f1, tmp_str, '#');
      SmartField::ICFlag = atoi(tmp_str.c_str());
      std::getline(f1, tmp_str, '\n');
    }
    else
    {
      std::cout << "cannot open params_ICs.in" << std::endl;
      exit(1);
    }
  }

  std::cout << "   - IC Parameters:\n" << std::scientific;
  std::cout << "      BinaryIntFlag      = " << SmartField::BinaryInFlag << std::endl;
  std::cout << "      BinaryOutFlag      = " << SmartField::BinaryOutFlag << std::endl;
  std::cout << "      ICFlag             = " << SmartField::ICFlag << std::endl;

  FieldStack initstack;
  SmartField::FDSize = fd_size;
  SmartField::NPW = initstack.CurrLayout->getNPWglobal();

} // }}}

void delete_SmartField_statics()
{ // {{{
} // }}}

void test_time()
{ // {{{
  
  struct timeval start, stop;
  double totaltime;
  int num_runs = 100000;
  std::cout << "\n*** SmartField Time Tests ***\n\n";
  std::cout << "Time results are from ";
  std::cout << num_runs;
  std::cout << "of runs\n\n";  
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

  //assignment test
  std::cout << "Assignment: \n";
  {
    SmartField tmp1, tmp2;
    
    gettimeofday(&start, NULL);
    for (int i = 0; i <= num_runs; ++i) {
      tmp1 = ftnum;
    }
    gettimeofday(&stop, NULL);
    totaltime = (stop.tv_sec - start.tv_sec) + 
                (stop.tv_usec - start.tv_usec)/1000000.0;
    totaltime = totaltime/num_runs;
    printf(" - Elapsed Time %10.3e (sec.)\n\n", totaltime );
  }
  // Compound Subtraction test
  std::cout << "Compound Subtraction: \n";
  {
    
    SmartField tmp1, tmp2;
    tmp1 = myfield;

    gettimeofday(&start, NULL);
    for (int i = 0; i <= num_runs; ++i) {
      tmp1 -= myfield;
    }
    gettimeofday(&stop, NULL);
    totaltime = (stop.tv_sec - start.tv_sec) + 
                (stop.tv_usec - start.tv_usec)/1000000.0;
    totaltime = totaltime/num_runs;
    printf(" - Elapsed Time %10.3e (sec.)\n\n", totaltime );
  }

  //Compound addition test
  std::cout << "Compound Addition: \n";
  {
    SmartField tmp1, tmp2;
    tmp1 = FieldType(0.);

    gettimeofday(&start, NULL);
    for (int i = 0; i <= num_runs; ++i) {
      tmp1 += myfieldvec;
    }
    gettimeofday(&stop, NULL);
    totaltime = (stop.tv_sec - start.tv_sec) + 
                (stop.tv_usec - start.tv_usec)/1000000.0;
    totaltime = totaltime/num_runs;
    printf(" - Elapsed Time %10.3e (sec.)\n\n", totaltime );
  }

  //Compound multiplication test
  std::cout << "Compound multiplication: \n";
  {
    SmartField tmp1, tmp2;
    tmp1 = FieldType(3.);
    tmp2 = FieldType(3.);

    gettimeofday(&start, NULL);
    for (int i = 0; i <= num_runs; ++i) {
      tmp1 *= tmp2;
    }
    gettimeofday(&stop, NULL);
    totaltime = (stop.tv_sec - start.tv_sec) + 
                (stop.tv_usec - start.tv_usec)/1000000.0;
    totaltime = totaltime/num_runs;
    printf(" - Elapsed Time %10.3e (sec.)\n\n", totaltime );
  }


  //Compound division test
  std::cout << "Compound division: \n";
  {
    SmartField tmp1, tmp2;
    tmp1 = FieldType(3.);
    
    gettimeofday(&start, NULL);
    for (int i = 0; i <= num_runs; ++i) {
      tmp1 /= ftnum;
    }
    gettimeofday(&stop, NULL);
    totaltime = (stop.tv_sec - start.tv_sec) + 
                (stop.tv_usec - start.tv_usec)/1000000.0;
    totaltime = totaltime/num_runs;
    printf(" - Elapsed Time %10.3e (sec.)\n\n", totaltime );
  }
  //skipped file I/O and copying

  //Getting and setting stuff
  std::cout << "Getters and setters and stuff:\n";
  SmartField tmp1;
  tmp1 = myfieldvec;
  SmartField tmp2 = FieldType(3.0, 4.0);
 
  UInt idx1 = SmartField::NPW/3+1;

  for (int j = 1; j <= 15; ++j) {
    std::cout << j;
    gettimeofday(&start, NULL);
    for (int i = 0; i <= num_runs; ++i) {
      switch (j)
      {
        case 1: tmp1.l1norm();
          break;
        case 2: tmp1.max();
          break;
        case 3: tmp1.maxsigned();
          break;
        case 4: tmp1.minsigned();
          break;
        case 5: tmp1.maxidx();
          break;
        case 6: tmp1.sumelem();
          break;
        case 7: tmp1.integrate();
          break;
        case 8: tmp1.integrate_square();
          break;
        case 9: tmp1.log(tmp1);
          break;
        case 10: tmp1.exponentiate(tmp1);
          break;
        case 11: tmp1.sqrt(tmp1);
          break;
        case 12: tmp2.l1norm();
          break;
        case 13: tmp2.l2norm();
          break;
        case 14: tmp1.setelement(FieldType(1024.), idx1);
          break;
        case 15: tmp1.getelement( idx1 );
          break;
      }
    }
    gettimeofday(&stop, NULL);
    totaltime = (stop.tv_sec - start.tv_sec) + 
              (stop.tv_usec - start.tv_usec)/1000000.0;
    totaltime = totaltime/num_runs;
    printf(" - Elapsed Time %10.3e (sec.)\n\n", totaltime );
  }

//copy test
  std::cout << "Copying test: \n";
  gettimeofday(&start, NULL);
  for (int i = 0; i <= num_runs/1000; ++i) {
    tmp2.copytoarray( myvecvec );
  }
  gettimeofday(&stop, NULL);
  totaltime = (stop.tv_sec - start.tv_sec) + 
            (stop.tv_usec - start.tv_usec)/1000000.0;
  totaltime = totaltime/(num_runs/1000);
  printf(" - Elapsed Time %10.3e (sec.)\n\n", totaltime );

} // }}}

void SmartField_unit_tests()
{ // {{{

  std::cout << "\n*** SmartField Unit Tests ***\n";
  std::cout << "   (All results should be zero)\n\n";

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

  std::cout << " * Getters and Setters (values are non-zero here)\n";
  {
    SmartField tmp1;
    tmp1 = myfieldvec;

    UInt idx1 = SmartField::NPW/3+1;
    UInt idx2 = idx1+1;
    UInt idx3 = idx2+1;

    tmp1.setelement(FieldType(1024.), idx1);
    tmp1.setelement(FieldType(-8.), idx2);
    tmp1.setelement(FieldType(-2048.), idx3);
    std::cout << "  - element " << idx1 << ": " << tmp1.getelement( idx1 ) << std::endl;
    std::cout << "  - element " << idx2 << ": " << tmp1.getelement( idx2 ) << std::endl;
    std::cout << "  - l1norm: " << tmp1.l1norm() << std::endl;
    std::cout << "  - max: " << tmp1.max() << std::endl;
    std::cout << "  - maxsigned: " << tmp1.maxsigned() << std::endl;
    std::cout << "  - minsigned: " << tmp1.minsigned() << std::endl;
    std::cout << "  - maxidx: " << tmp1.maxidx() << std::endl;
    std::cout << "  - sum: " << tmp1.sumelem() << std::endl;
    std::cout << "  - integrate: " << tmp1.integrate() << std::endl;
    std::cout << "  - integrate square: " << tmp1.integrate_square() << std::endl;

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

  std::cout << "\n*** End Unit Tests ***\n";
  
} // }}}

// ----------------- Constructor/Destructor -----------------

SmartField :: SmartField() // get some scratch space off of the stack
: data_ptr(FDSize, NULL)
{ // {{{

  CreateField();

} // }}}

SmartField :: ~SmartField() // put the scratch space back in the stack
{ // {{{

  for (UInt m=0; m<FDSize; m++)
  {
    FieldStack::RetireField( data_ptr[m] );
  }

} // }}}

SmartField :: SmartField( const SmartField & rhs ) // copy constructor
: data_ptr(FDSize, NULL)
{ // {{{

  CreateField();
  setflag_inrealspace(rhs.getflag_inrealspace());
  for (UInt m = 0; m<FDSize; m++)
  {
    *(data_ptr[m])=*(rhs.data_ptr[m]); 
  }

} // }}}

SmartField :: SmartField( const std::vector< Field<FieldType> > &rhs )
: data_ptr(FDSize, NULL)
{ // {{{

  CreateField();

  if ( rhs.size() != FDSize )
  {
    std::cout << "Error. The size of the LHS (SmartField) and the RHS (std::vector) do not match." << std::endl;
    exit(1);
  }
  else
  {
    for (UInt m = 0; m<FDSize; m++)
    {
      data_ptr[m]->setflag_inrealspace(rhs[m].getflag_inrealspace());
      *(data_ptr[m])=rhs[m]; 
    }
  }

} // }}}

SmartField :: SmartField( const std::vector< std::vector<FieldType> > &rhs )
: data_ptr(FDSize, NULL)
{ // {{{

  CreateField();

  if ( rhs.size() != FDSize )
  {
    std::cout << "Error. The size of the LHS (SmartField) and the RHS (std::vector) do not match." << std::endl;
    exit(1);
  }
  else
  {
    for (UInt m = 0; m<FDSize; m++)
    {
      *(data_ptr[m])=rhs[m]; 
    }
  }

} // }}}

SmartField :: SmartField( const Field<FieldType> &rhs )
: data_ptr(FDSize, NULL)
{ // {{{

  CreateField();

  for (UInt m = 0; m<FDSize; m++)
  {
    *(data_ptr[m])=rhs; 
  }

} // }}}

SmartField :: SmartField( const std::vector<FieldType> &rhs )
: data_ptr(FDSize, NULL)
{ // {{{

  CreateField();

  for (UInt m = 0; m<FDSize; m++)
  {
    *(data_ptr[m])=rhs; 
  }

} // }}}

SmartField :: SmartField( const FieldType rhs)
: data_ptr(FDSize, NULL)
{ // {{{

  CreateField();
  for (UInt m = 0; m<FDSize; m++)
  {
    *(data_ptr[m])=rhs; 
  }

} // }}}

void SmartField::CreateField()
{ // {{{

  for (UInt m=0; m<FDSize; m++)
  {
    FieldStack::FetchField( data_ptr[m] );
  }

} // }}}

// ----------------- Operators -----------------

// assignment
SmartField & SmartField::operator= (const SmartField &rhs)
{ // {{{

  if (this != &rhs) 
  { 
    setflag_inrealspace(rhs.getflag_inrealspace());
    for (UInt m = 0; m<FDSize; m++)
    {
      *(data_ptr[m])=*(rhs.data_ptr[m]); 
    }
  }
  return *this;

} // }}}

SmartField & SmartField::operator= (const std::vector< Field<FieldType> > &rhs)
{ // {{{
  
  if ( rhs.size() != FDSize )
  {
    std::cout << "Error. The size of the LHS (SmartField) and the RHS (std::vector) do not match." << std::endl;
    exit(1);
  }
  else
  {
    for (UInt m = 0; m<FDSize; m++)
    {
      data_ptr[m]->setflag_inrealspace(rhs[m].getflag_inrealspace());
      *(data_ptr[m])=rhs[m]; 
    }
  }
  return *this;

} // }}}

SmartField & SmartField::operator= (const std::vector< Field<RealType> > &rhs)
{ // {{{
  
  if ( rhs.size() != FDSize )
  {
    std::cout << "Error. The size of the LHS (SmartField) and the RHS (std::vector) do not match." << std::endl;
    exit(1);
  }
  else
  {
    for (UInt m = 0; m<FDSize; m++)
    {
      data_ptr[m]->setflag_inrealspace(rhs[m].getflag_inrealspace());
      *(data_ptr[m])=rhs[m];
    }
  }
  return *this;

} // }}}

SmartField & SmartField::operator= (const Field<FieldType> &rhs)
{ // {{{

  this->setflag_inrealspace(rhs.getflag_inrealspace());
  for (UInt m = 0; m<FDSize; m++)
  {
    *(data_ptr[m])=rhs;
  }
  return *this;

} // }}}

SmartField & SmartField::operator= (const std::vector< std::vector<FieldType> > &rhs)
{ // {{{
  
  if ( rhs.size() != FDSize )
  {
    std::cout << "Error. The size of the LHS (SmartField) and the RHS (std::vector) do not match." << std::endl;
    exit(1);
  }
  else
  {
    for (UInt m = 0; m<FDSize; m++)
    {
      *(data_ptr[m])=rhs[m]; 
    }
  }
  return *this;

} // }}}

SmartField & SmartField::operator= (const std::vector<FieldType> &rhs)
{ // {{{

  for (UInt m = 0; m<FDSize; m++)
  {
    *(data_ptr[m])=rhs; 
  }
  return *this;

} // }}}

SmartField & SmartField::operator= (const FieldType rhs)
{ // {{{

  for (UInt m = 0; m<FDSize; m++)
  {
    *(data_ptr[m])=rhs; 
  }
  return *this;

} // }}}


// compound addition
SmartField & SmartField::operator+=(const SmartField &rhs)
{ // {{{

  for (UInt m = 0; m<FDSize; m++)
  {
    *(data_ptr[m]) += *(rhs.data_ptr[m]); 
  }
  return *this;

} // }}}

SmartField & SmartField::operator+=(const std::vector< Field<FieldType> > &rhs)
{ // {{{
  
  if ( rhs.size() != FDSize )
  {
    std::cout << "Error. The size of the LHS (SmartField) and the RHS (std::vector) do not match." << std::endl;
    exit(1);
  }
  else
  {
    for (UInt m = 0; m<FDSize; m++)
    {
      *(data_ptr[m])+=rhs[m]; 
    }
  }
  return *this;

} // }}}

SmartField & SmartField::operator+=(const Field<FieldType> &rhs)
{ // {{{

  for (UInt m = 0; m<FDSize; m++)
  {
    *(data_ptr[m]) += rhs; 
  }
  return *this;

} // }}}

SmartField & SmartField::operator+=(FieldType rhs)
{ // {{{

  for (UInt m = 0; m<FDSize; m++)
  {
    *(data_ptr[m]) += rhs; 
  }
  return *this;

} // }}}

// compound subtraction
SmartField & SmartField::operator-=(const SmartField &rhs)
{ // {{{

  for (UInt m = 0; m<FDSize; m++)
  {
    *(data_ptr[m]) -= *(rhs.data_ptr[m]); 
  }
  return *this;

} // }}}

SmartField & SmartField::operator-=(const std::vector< Field<FieldType> > &rhs)
{ // {{{
  
  if ( rhs.size() != FDSize )
  {
    std::cout << "Error. The size of the LHS (SmartField) and the RHS (std::vector) do not match." << std::endl;
    exit(1);
  }
  else
  {
    for (UInt m = 0; m<FDSize; m++)
    {
      *(data_ptr[m])-=rhs[m]; 
    }
  }
  return *this;

} // }}}

SmartField & SmartField::operator-=(const Field<FieldType> &rhs)
{ // {{{

  for (UInt m = 0; m<FDSize; m++)
  {
    *(data_ptr[m]) -= rhs;
  }
  return *this;

} // }}}

SmartField & SmartField::operator-=(const FieldType rhs)
{ // {{{

  for (UInt m = 0; m<FDSize; m++)
  {
    *(data_ptr[m]) += (-rhs); 
  }
  return *this;

} // }}}

// compound multiplication
SmartField & SmartField::operator*=(const SmartField &rhs)
{ // {{{

  for (UInt m = 0; m<FDSize; m++)
  {
    *(data_ptr[m]) *= *(rhs.data_ptr[m]); 
  }
  return *this;

} // }}}

SmartField & SmartField::operator*=(const std::vector< Field<FieldType> > &rhs)
{ // {{{
  
  if ( rhs.size() != FDSize )
  {
    std::cout << "Error. The size of the LHS (SmartField) and the RHS (std::vector) do not match." << std::endl;
    exit(1);
  }
  else
  {
    for (UInt m = 0; m<FDSize; m++)
    {
      *(data_ptr[m])*=rhs[m]; 
    }
  }
  return *this;

} // }}}

SmartField & SmartField::operator*=(const Field<FieldType> &rhs)
{ // {{{

  for (UInt m = 0; m<FDSize; m++)
  {
    *(data_ptr[m]) *= rhs;
  }
  return *this;

} // }}}

SmartField & SmartField::operator*=(FieldType rhs)
{ // {{{

  for (UInt m = 0; m<FDSize; m++)
  {
    *(data_ptr[m]) *= rhs; 
  }
  return *this;

} // }}}

SmartField & SmartField::operator*=(RealType rhs)
{ // {{{

  for (UInt m = 0; m<FDSize; m++)
  {
    *(data_ptr[m]) *= rhs;
  }
  return *this;

} // }}}

SmartField & SmartField::prodconjg(const SmartField &rhs)
{ // {{{

  for (UInt m = 0; m<FDSize; m++)
  {
    data_ptr[m]->prodconjg(*(rhs.data_ptr[m]));
  }
  return *this;

} // }}}

SmartField & SmartField::prodconjg(const Field<FieldType> &rhs)
{ // {{{

  for (UInt m = 0; m<FDSize; m++)
  {
    data_ptr[m]->prodconjg(rhs);
  }
  return *this;

} // }}}

// compound division
SmartField & SmartField::operator/=(const SmartField &rhs)
{ // {{{

  for (UInt m = 0; m<FDSize; m++)
  {
    *(data_ptr[m]) /= *(rhs.data_ptr[m]); 
  }
  return *this;

} // }}}

SmartField & SmartField::operator/=(const std::vector< Field<FieldType> > &rhs)
{ // {{{
  
  if ( rhs.size() != FDSize )
  {
    std::cout << "Error. The size of the LHS (SmartField) and the RHS (std::vector) do not match." << std::endl;
    exit(1);
  }
  else
  {
    for (UInt m = 0; m<FDSize; m++)
    {
      *(data_ptr[m])/=rhs[m]; 
    }
  }
  return *this;

} // }}}

SmartField & SmartField::operator/=(const Field<FieldType> &rhs)
{ // {{{

  for (UInt m = 0; m<FDSize; m++)
  {
    *(data_ptr[m]) /= rhs;
  }
  return *this;

} // }}}

SmartField & SmartField::operator/=(const FieldType rhs)
{ // {{{

  for (UInt m = 0; m<FDSize; m++)
  {
    *(data_ptr[m]) *= (1./rhs); 
  }
  return *this;

} // }}}

// ----------------- Member Functions -----------------

// like subscript operator, [], to get a whole array (for FD)
void SmartField::subscript( const SmartField & rhs, const std::vector<UInt> idx )
{ // {{{

  if (idx.size() != FDSize )
  {
    std::cout << "***(Error in SmartField.subscript)***";
    std::cout << "SmartField and idx sizes are not compatible";
    exit(1);
  }

  if (&rhs == this)
  {
    std::cout << "***(Error in SmartField.subscript)***";
    std::cout << "Cannot do in place operation.";
    exit(1);
  }
  
  for (UInt m=0; m<FDSize; m++)
  {
    *(data_ptr[m]) = rhs[idx[m]];
  }
  
} // }}}

// get a slice of an array (for FD)
void SmartField::slice( const SmartField & rhs, const std::vector<UInt> idx )
{ // {{{

  if (idx.size() > FDSize )
  {
    std::cout << "***(Error in SmartField.slice)***";
    std::cout << "idx size too big for SmartField";
    exit(1);
  }

  if (&rhs == this)
  {
    std::cout << "***(Error in SmartField.slice)***";
    std::cout << "Cannot do in place operation.";
    exit(1);
  }
  
  for (UInt m=0; m<idx.size(); m++)
  {
    *(data_ptr[idx[m]]) = rhs[idx[m]];
  }
  
} // }}}

// (mapping to Field class member functions)

// initialization
void SmartField::fill(enum_fieldinitmethod filltype, Vec1dReal parameters, bool screenoutputtype)
{ //{{{
  for( UInt m=0; m<FDSize; m++)
  {
    data_ptr[m]->fill(filltype, parameters, screenoutputtype); 
  }
} // }}}

void SmartField::fillnoise(UInt noisetype)
{ //{{{
  for( UInt m=0; m<FDSize; m++)
  {
    data_ptr[m]->fillnoise(noisetype);
  }
} // }}}

// fft operations
void SmartField::fft_rtok( bool applyscale /*=true*/ )
{ //{{{
  for( UInt m=0; m<FDSize; m++)
  {
    data_ptr[m]->fft_rtok(applyscale);
  }
} // }}}

void SmartField::fft_rtonegk( bool applyscale /*=true*/ )
{ //{{{
  for( UInt m=0; m<FDSize; m++)
  {
    data_ptr[m]->fft_rtonegk(applyscale);
  }
} // }}}

void SmartField::fft_ktor()
{ //{{{
  for( UInt m=0; m<FDSize; m++)
  {
    data_ptr[m]->fft_ktor();
  }
} // }}}

// field manipulation/arithmetic setters
void SmartField::setaverage(const FieldType &value)
{ //{{{
  for( UInt m=0; m<FDSize; m++)
  {
    data_ptr[m]->setaverage(value);
  }
} // }}}

void SmartField::conjg()
{ //{{{
  for( UInt m=0; m<FDSize; m++)
  {
    data_ptr[m]->conjg();
  }
} // }}}

void SmartField::zero()
{ //{{{
  for( UInt m=0; m<FDSize; m++)
  {
    data_ptr[m]->zero();
  }
} // }}}

void SmartField::zeroreal()
{ //{{{
  for( UInt m=0; m<FDSize; m++)
  {
    data_ptr[m]->zeroreal();
  }
} // }}}

void SmartField::zeroimag()
{ //{{{
  for( UInt m=0; m<FDSize; m++)
  {
    data_ptr[m]->zeroimag();
  }
} // }}}

void SmartField::zeropos()
{ //{{{
  for( UInt m=0; m<FDSize; m++)
  {
    data_ptr[m]->zeropos();
  }
} // }}}

void SmartField::zeroneg()
{ //{{{
  for( UInt m=0; m<FDSize; m++)
  {
    data_ptr[m]->zeroneg();
  }
} // }}}

void SmartField::zero_k0()
{ // {{{
  for( UInt m=0; m<FDSize; m++)
  {
    // set the k=0 mode to zero
    //data_ptr[m]->setelement(FieldType(0.), 0);
    data_ptr[m]->setaverage(FieldType(0.));
  }
} // }}}

// data manipulation setters
void SmartField::setName(std::string szName_in)
{ //{{{
  for( UInt m=0; m<FDSize; m++)
  {
    data_ptr[m]->setName(szName_in);
  }
} // }}}

void SmartField::setflag_inrealspace(bool bInRealSpace)
{ //{{{
  for( UInt m=0; m<FDSize; m++)
  {
    data_ptr[m]->setflag_inrealspace(bInRealSpace);
  }
} // }}}

void SmartField::axpby_inplace(const SmartField &y, RealType a, RealType b)
{ //{{{
  for( UInt m=0; m<FDSize; m++)
  {
    data_ptr[m]->axpby_inplace(y.GetField(m), a, b);
  }
} // }}}

void SmartField::axpby_inplace(const SmartField &y, Vec1dReal a, Vec1dReal b)
{ //{{{
  for( UInt m=0; m<FDSize; m++)
  {
    data_ptr[m]->axpby_inplace(y.GetField(m), a[m], b[m]);
  }
} // }}}

void SmartField::xpby_inplace(const SmartField &y, RealType b)
{ //{{{
  for( UInt m=0; m<FDSize; m++)
  {
    data_ptr[m]->xpby_inplace(y.GetField(m),b);
  }
} // }}}

void SmartField::xpby_inplace(const SmartField &y, Vec1dReal b)
{ //{{{
  for( UInt m=0; m<FDSize; m++)
  {
    data_ptr[m]->xpby_inplace(y.GetField(m),b[m]);
  }
} // }}}

void SmartField::axpy_inplace(const SmartField &y, RealType a)
{ //{{{
  for( UInt m=0; m<FDSize; m++)
  {
    data_ptr[m]->axpy_inplace(y.GetField(m),a);
  }
} // }}}

void SmartField::axpy_inplace(const SmartField &y, Vec1dReal a)
{ //{{{
  for( UInt m=0; m<FDSize; m++)
  {
    data_ptr[m]->axpy_inplace(y.GetField(m),a[m]);
  }
} // }}}

void SmartField::axpby_inplace(const SmartField &y, FieldType a, FieldType b)
{ //{{{
  for( UInt m=0; m<FDSize; m++)
  {
    data_ptr[m]->axpby_inplace(y.GetField(m), a, b);
  }
} // }}}

void SmartField::axpby_inplace(const SmartField &y, Vec1dFieldType a, Vec1dFieldType b)
{ //{{{
  for( UInt m=0; m<FDSize; m++)
  {
    data_ptr[m]->axpby_inplace(y.GetField(m), a[m], b[m]);
  }
} // }}}

void SmartField::xpby_inplace(const SmartField &y, FieldType b)
{ //{{{
  for( UInt m=0; m<FDSize; m++)
  {
    data_ptr[m]->xpby_inplace(y.GetField(m),b);
  }
} // }}}

void SmartField::xpby_inplace(const SmartField &y, Vec1dFieldType b)
{ //{{{
  for( UInt m=0; m<FDSize; m++)
  {
    data_ptr[m]->xpby_inplace(y.GetField(m),b[m]);
  }
} // }}}

void SmartField::axpy_inplace(const SmartField &y, FieldType a)
{ //{{{
  for( UInt m=0; m<FDSize; m++)
  {
    data_ptr[m]->axpy_inplace(y.GetField(m),a);
  }
} // }}}

void SmartField::axpy_inplace(const SmartField &y, Vec1dFieldType a)
{ //{{{
  for( UInt m=0; m<FDSize; m++)
  {
    data_ptr[m]->axpy_inplace(y.GetField(m),a[m]);
  }
} // }}}


void SmartField::accumulateproduct_inplace(SmartField const &in1, SmartField const &in2, double scale)
{ //{{{
  for( UInt m=0; m<FDSize; m++)
  {
    data_ptr[m]->accumulateproduct_inplace(in1.GetField(m), in2.GetField(m), scale);
  }
} // }}}

void SmartField::accumulateproduct_inplace(SmartField const &in1, SmartField const &in2)
{ //{{{
  for( UInt m=0; m<FDSize; m++)
  {
    data_ptr[m]->accumulateproduct_inplace(in1.GetField(m), in2.GetField(m));
  }
} // }}}

void SmartField::exponentiate(const SmartField &in, RealType scale)
{ //{{{
  for( UInt m=0; m<FDSize; m++)
  {
    data_ptr[m]->exponentiate(in.GetField(m), scale);
  }
} // }}}

void SmartField::exponentiate(const SmartField &in)
{ //{{{
  for( UInt m=0; m<FDSize; m++)
  {
    data_ptr[m]->exponentiate(in.GetField(m));
  }
} // }}}

void SmartField::sqrt(const SmartField &in)
{ //{{{
  for( UInt m=0; m<FDSize; m++)
  {
    data_ptr[m]->sqrt(in.GetField(m));
  }
} // }}}

void SmartField::log(const SmartField &in, RealType scale)
{ //{{{
  for( UInt m=0; m<FDSize; m++)
  {
    data_ptr[m]->log(in.GetField(m), scale);
  }
} // }}}

void SmartField::log(const SmartField &in)
{ //{{{
  for( UInt m=0; m<FDSize; m++)
  {
    data_ptr[m]->log(in.GetField(m));
  }
} // }}}

void SmartField::logreal(const SmartField &in, RealType scale)
{ //{{{
  for( UInt m=0; m<FDSize; m++)
  {
    data_ptr[m]->logreal(in.GetField(m), scale);
  }
} // }}}

void SmartField::logreal(const SmartField &in)
{ //{{{
  for( UInt m=0; m<FDSize; m++)
  {
    data_ptr[m]->logreal(in.GetField(m));
  }
} // }}}

void SmartField::setelement(const FieldType &value, UInt indx)
{ // {{{
  
  // use a global index valid for either PS or FD field
  UInt i, j;
  i = std::floor(indx/NPW);
  j = indx - i*NPW;
  data_ptr[i]->setelement(value, j);

} // }}}

// data manipulation getters
bool SmartField::getflag_inrealspace() const
{ //{{{
  bool realspaceflag, flagref;
    
  realspaceflag = false;
  for( UInt m=0; m<FDSize; m++)
  {
    flagref = realspaceflag;
    realspaceflag = data_ptr[m]->getflag_inrealspace();
    if (m > 1 and flagref != realspaceflag)
    {
      std::cout << "Error: flag_inrealspace is different for different field elements.\n";
    }
  }
  return realspaceflag;
} // }}}

const FFTlayout& SmartField::getFFTlayout() const
{ //{{{
    return data_ptr[0]->getFFTlayout();
} // }}}

FieldType SmartField::getelement(UInt indx) const
{ //{{{

  // use a global index valid for either PS or FD field
  UInt i, j;
  i = std::floor(indx/NPW);
  j = indx - i*NPW;
  return data_ptr[i]->getelement(j);

} // }}}

void SmartField::copytoarray(std::vector< std::vector<FieldType> > &out) const
{ //{{{
  for( UInt m=0; m<FDSize; m++)
  {
    data_ptr[m]->copytoarray(out[m]);
  }
} // }}}

void SmartField::globalcopytoarray(std::vector< std::vector<FieldType> > &out) const
{ //{{{
  for( UInt m=0; m<FDSize; m++)
  {
    data_ptr[m]->globalcopytoarray(out[m]);
  }
} // }}}

RealType SmartField::l1norm() const
{ //{{{

  RealType l1sum = 0.;
  for( UInt m=0; m<FDSize; m++)
  {
    l1sum += std::abs( data_ptr[m]->l1norm() );
  }
  
  return l1sum;
  
} // }}}

RealType SmartField::l2norm() const
{ // {{{

  RealType l2sum = 0.;
  for( UInt m=0; m<FDSize; m++)
  {
    l2sum += std::pow( data_ptr[m]->l2norm(), 2);
  }
  l2sum = std::sqrt(l2sum);
  
  return l2sum;
} // }}}

FieldType SmartField::max() const
{ // {{{

  // returns raw element of maximum abs value
  RealType fieldmax = 0.;
  RealType localmax = 0.;
  UInt maxelem = 0;

  for( UInt m=0; m<FDSize; m++)
  {
    localmax = std::abs( data_ptr[m]->max() );
    if (localmax > fieldmax)
    {
      fieldmax = localmax;
      maxelem = m;
    }
  }
  return data_ptr[maxelem]->max();

} // }}}

FieldType SmartField::maxsigned() const
{ // {{{

  // returns element with maximum real value (respecting sign)
  FieldType fieldmax = data_ptr[0]->maxsigned();
  for( UInt m=1; m<FDSize; m++)
  {
    FieldType localmax = data_ptr[m]->maxsigned();
    if ( localmax.real() > fieldmax.real() )
      fieldmax = localmax;
  }
  return fieldmax;

} // }}}

FieldType SmartField::minsigned() const
{ //{{{

  // returns element with minimum real value (respecting sign)
  FieldType fieldmin = data_ptr[0]->minsigned();
  for( UInt m=1; m<FDSize; m++)
  {
    FieldType localmin = data_ptr[m]->minsigned();
    if ( localmin.real() < fieldmin.real() )
      fieldmin = localmin;
  }
  return fieldmin;

} // }}}

UInt SmartField::maxidx() const
{ //{{{

  // returns index of the element with maximum real value (respecting sign)
  UInt maxslice = 0;
  if(FDSize > 1)
  {
    // Determine which slice contains the maxval
    FieldType fieldmax = data_ptr[0]->maxsigned();
    for( UInt m=1; m<FDSize; m++)
    {
      FieldType localmax = data_ptr[m]->maxsigned();
      if ( localmax.real() > fieldmax.real() )
      {
        fieldmax = localmax;
        maxslice = m;
      }
    }
  }

  // return a global index valid for either PS or FD
  return maxslice*NPW + data_ptr[maxslice]->maxidx();

} // }}}

FieldType SmartField::getaverage() const
{ //{{{

  // (Thus is the same as integrate, but getaverage 
  //  is for k-space and integrate is for real-space)

  FieldType fieldavg = 0.;
  for( UInt m=0; m<FDSize; m++)
  {
    fieldavg += data_ptr[m]->getaverage();
  }
  fieldavg /= FDSize;

  return fieldavg;

} // }}}

FieldType SmartField::integrate() const
{ //{{{

  // Similar to Kris' definition: 
  // Integration => V/M * sum_i f_i; V is NOT INCLUDED
  // (Thus is the same as getaverage, but getaverage 
  //  is for k-space, this is for realspace.)

  FieldType fieldint = 0.;
  for( UInt m=0; m<FDSize; m++)
  {
    fieldint += data_ptr[m]->integrate();
  }
  fieldint /= FDSize;

  return fieldint;

} // }}}

FieldType SmartField::integrate_square() const
{ //{{{

  FieldType fieldint = 0.;
  for( UInt m=0; m<FDSize; m++)
  {
    fieldint += data_ptr[m]->integrate_square();
  }
  fieldint /= FDSize;

  return fieldint;

} // }}}

std::string SmartField::getName() const
{ //{{{
  // only find the name of the first element
  return data_ptr[0]->getName();
} // }}}

FieldType SmartField::sumelem() const
{ //{{{

  FieldType fieldsum = 0.;
  for( UInt m=0; m<FDSize; m++)
  {
    fieldsum += data_ptr[m]->sumelem();
  }
  
  return fieldsum;

} // }}}

FieldType SmartField::sumelem_weighted(const SmartField &weights) const
{ //{{{
  
  FieldType fieldsum = 0.;
  for( UInt m=0; m<FDSize; m++)
  {
    fieldsum += data_ptr[m]->sumelem_weighted(weights.GetField(m));
  }
  
  return fieldsum;

} // }}}

FieldType SmartField::sumelem_weighted(const Field<FieldType> &weights) const
{ //{{{

  FieldType fieldsum = 0.;
  for( UInt m=0; m<FDSize; m++)
  {
    data_ptr[m]->sumelem_weighted(weights);
  }
  
  return fieldsum;

} // }}}

// IO
void SmartField::writefield(std::string filename, bool omitimaginary)
{ //{{{

  bool iostatus=false;

  std::list<std::string> filenames;
  filenames.push_back(filename);

  if (BinaryOutFlag)
  {
    WriteFieldsUnformatted( filenames, data_ptr, iostatus);
  }
  else
  {
    WriteFields( filenames, data_ptr, omitimaginary, iostatus);
  }

} // }}}

void SmartField::readfield(std::string filename)
{ //{{{

  bool iostatus=false;
  ReadFields( filename, data_ptr, iostatus);

} // }}}

void SmartField::debug_print_nelements(UInt num) const
{ //{{{
  for( UInt m=0; m<FDSize; m++)
  {
    data_ptr[m]->debug_print_nelements(num);
  }
} // }}}

