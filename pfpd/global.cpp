/**********************************************************************************
 *
 * PolyFTS Project
 *
 * File created by Kris Delaney on 2011-10-05.
 * Copyright (c) 2011 University of California, Santa Barbara. All rights reserved.
 *
 **********************************************************************************/
#include "global.h"
#include "matrix_utils.h"

// IN THE GLOBAL SCOPE:
// Useful global functions
bool CompareEqual(double A, double B, double RelativeTol/*=-1.0*/)
{
  // Take the difference
  double diff = fabs(A-B);
  // Set up tolerance
  double T;
  if(RelativeTol < 0.)
    T = 100. * std::numeric_limits<double>::epsilon(); // Allow up to 100 epsilon difference by default
  else
    T = RelativeTol;

  // Compare for absolute difference close to zero
  if(diff<=T)
    return true;

  // Test based on RELATIVE error
  // Find the largest abs
  A = fabs(A);
  B = fabs(B);
  double largest = (B > A) ? B : A;
  if(diff <= largest * T)
    return true;

  return false;
}
bool CompareEqual(std::complex<double> A, double B, double tol/*=-1.0*/)
{
  return(CompareEqual(A.real(), B, tol) && CompareEqual(A.imag(),0.0, tol));
}
bool CompareEqual(std::complex<double> A, std::complex<double> B, double tol/*=-1.0*/)
{
  return(CompareEqual(A.real(), B.real(), tol) && CompareEqual(A.imag(), B.imag(), tol));
}
bool CompareEqual(double A, std::complex<double> B, double tol/*=-1.0*/)
{
  return(CompareEqual(A, B.real(), tol) && CompareEqual(0.0, B.imag(), tol));
}

bool CompareEqual(float A, float B, float RelativeTol/*=-1.0f*/)
{
  // Take the difference
  float diff = fabsf(A-B);
  // Set up tolerance
  float T;
  if(RelativeTol < 0.)
    T = 100.f * std::numeric_limits<float>::epsilon(); // Allow up to 100 epsilon difference by default
  else
    T = RelativeTol;

  // Compare for absolute difference close to zero
  if(diff<=T)
    return true;

  // Test based on RELATIVE error
  // Find the largest abs
  A = fabsf(A);
  B = fabsf(B);
  float largest = (B > A) ? B : A;
  if(diff <= largest * T)
    return true;

  return false;
}
bool CompareEqual(std::complex<float> A, float B, float tol/*=-1.0f*/)
{
  return(CompareEqual(A.real(), B, tol) && CompareEqual(A.imag(),0.0f, tol));
}
bool CompareEqual(std::complex<float> A, std::complex<float> B, float tol/*=-1.0f*/)
{
  return(CompareEqual(A.real(), B.real(), tol) && CompareEqual(A.imag(), B.imag(), tol));
}
bool CompareEqual(float A, std::complex<float> B, float tol/*=-1.0f*/)
{
  return(CompareEqual(A, B.real(), tol) && CompareEqual(0.0f, B.imag(), tol));
}

bool CompareEqual(UInt A, UInt B, float tol/*=-1.0f*/)
{
  if (A == B) 
    return true;
  else
    return false;
}

double gettime()
{
  struct timeval tv;
  gettimeofday(&tv,0);
  return(static_cast<double>(tv.tv_sec) + static_cast<double>(tv.tv_usec)*1e-6);
}

// --- DRT edits ---

Vec2dFieldType matmul( Vec2dFieldType A, Vec2dFieldType B)
{ // {{{

  Vec2dFieldType C(A.size(), Vec1dFieldType(B[0].size(), FieldType(0.)));
  
  for (UInt i=0; i<A.size(); ++i) // over rows in A
  {
    for (UInt j=0; j<B[0].size(); ++j) // columns in B
    {
      C[i][j] = FieldType(0.);
      for (UInt k=0; k<A[0].size(); ++k) // columns in A (or rows in B)
      {
        C[i][j] += A[i][k]*B[k][j];
      }
    }
  }
  return C;
} // }}}

Vec2dFieldType eye( UInt Nrow, UInt Ncol, FieldType scalar)
{ // {{{

  Vec2dFieldType out(Nrow, Vec1dFieldType(Ncol, FieldType(0.)));
  
  for (UInt i=0; i<Nrow; ++i)
  {
    for (UInt j=0; j<Ncol; ++j)
    {
      if ( i == j )
      {
        out[i][j] = scalar;
      }
      else
      {
        out[i][j] = FieldType(0.);
      }
    }
  }

  return out;

} // }}}

Vec2dFieldType kroneckerdelta( Vec1dReal A )
{ // {{{

  Vec2dFieldType B(A.size(), Vec1dFieldType(A.size(), FieldType(0.)));

  for (UInt i=0; i<A.size(); ++i) // over rows in A
  {
    B[i][i] = A[i];
  }
  return B;

} // }}}

/*
void make_posdef( Vec2dFieldType & A )
{ // {{{

  // Use LAPACK to do an eigenvalue decomposition of the matrix
  // check for negative eigenvalues, and then put the matrix
  // back together with only non-negative eigenvalues

  // Recall A = X Lamda X^{-1}, X is matrix with columns=eigenvectors

  lapack_int N = A.size();
  lapack_int info;
  char jobvl = 'N'; // no left eigenvectors
  char jobvr = 'V'; // yes right eigenvectors

  lapack_complex_double *AA = new lapack_complex_double[N*N]; // A
  lapack_complex_double *w = new lapack_complex_double[N]; // eigenvalues
  lapack_complex_double *vl = new lapack_complex_double[N*N]; 
  lapack_complex_double *vr = new lapack_complex_double[N*N]; // (right) eigenvectors
  lapack_complex_double *vrinv = new lapack_complex_double[N*N]; // vr^{-1}
  lapack_complex_double *lam= new lapack_complex_double[N*N]; // diagonal matrix of evals
  lapack_complex_double *lam_vrinv= new lapack_complex_double[N*N]; // lam.vr^{-1}
  lapack_int * ipiv = new lapack_int[N]; // pivots for LU

  // call lapack function for finding eigenvalues and eigenvectors
  // http://www.netlib.org/lapack/double/dgeev.f
  // https://software.intel.com/en-us/node/521147

  for (int i=0; i<N; i++)
  {
    for (int j=0; j<N; j++)
    {
      AA[i*N+j].real( A[i][j].real() );
      AA[i*N+j].imag( A[i][j].imag() );
    }
  }

  info = LAPACKE_zgeev(LAPACK_ROW_MAJOR, jobvl, jobvr, N, AA, N, w, vl, N, vr, N);

  if (info != 0)
  {
    std::cout << "*** Error in global.cpp:make_posdef ***\n";
    std::cout << "There was some problem with the LAPACK eigenvalue function.\n";
    std::cout << "info = " << info << std::endl;

    std::cout << "MATRIX: " << std::endl;
    for (int i=0; i<N; i++)
    {
      for (int j=0; j<N; j++)
        std::cout << A[i][j] << "\t";
      std::cout << std::endl;
    }

    exit(1);
  }

  // Find LU decompoistion of A:
  // http://www.netlib.org/lapack/double/dgetrf.f
  // https://software.intel.com/en-us/node/520877
  // and then get the matrix inverse:
  // http://www.netlib.org/lapack/double/dgetri.f
  // https://software.intel.com/en-us/node/520946

  for (int k=0; k<N*N; k++)
  {
    vrinv[k] = vr[k];
  }
  
  info = LAPACKE_zgetrf(LAPACK_ROW_MAJOR, N, N, vrinv, N, ipiv);

  if (info != 0)
  {
    std::cout << "*** Error in global.cpp:make_posdef ***\n";
    std::cout << "There was a problem with the LAPACK LU decomposition.\n";
    std::cout << "info = " << info << std::endl;
    exit(1);
  }

  info = LAPACKE_zgetri(LAPACK_ROW_MAJOR, N, vrinv, N, ipiv);

  if (info != 0)
  {
    std::cout << "*** Error in global.cpp:make_posdef ***\n";
    std::cout << "There was a problem with the LAPACK LU decomposition.\n";
    std::cout << "info = " << info << std::endl;
    exit(1);
  }

  // check that all eigenvalues are positive
  for (int i=0; i<N; i++)
  {
    w[i] = std::max(w[i].real(), 0.);
  }

  // rebuild matrix to ensure that it is positive semi-definite
  // A = vr . lam . vr_inv

  // 1) build lam
  for (int i=0; i<N; i++)
  {
    for (int j=0; j<N; j++)
    {
      if (i==j)
      {
        lam[i*N+j] = w[i];
      }
      else
      {
        lam[i*N+j] = 0.;
      }
    }
  } 

  // 2) w.vr_inv
  for (int i=0; i<N; i++)
  {
    for (int j=0; j<N; j++)
    {
      lam_vrinv[i*N+j]= 0.;
      for (int k=0; k<N; k++)
      {
        lam_vrinv[i*N+j] += lam[i*N+k]*vrinv[k*N+j];
      }
    }
  }

  RealType tmp;
  // 3) vr.(w.vr_inv)
  for (int i=0; i<N; i++)
  {
    for (int j=0; j<N; j++)
    {
      A[i][j] = 0;
      for (int k=0; k<N; k++)
      {
        tmp = A[i][j].real();
        tmp += (vr[i*N+k] * lam_vrinv[k*N+j]).real();
        A[i][j].real(tmp);
 
        tmp = A[i][j].imag();
        tmp += (vr[i*N+k] * lam_vrinv[k*N+j]).imag();
        A[i][j].imag(tmp);
      }
    }
  }

  delete[] AA;
  delete[] w;
  delete[] vl;
  delete[] vr;
  delete[] vrinv;
  delete[] lam;
  delete[] lam_vrinv;
  delete[] ipiv;

} // }}}
*/

void make_posdef_ReSymm( Vec2dFieldType & A )
{ // {{{
  // This is specifically for real symmetric matrices, and does not use LAPACK
  // For a real symmetric matrix, the matrix of eigenvectors is orthogonal, and no inversion is required.

  // Recall A = O d O^{T}, O is matrix with columns=eigenvectors of A

  int N = A.size();

  std::vector<std::vector<double> > ARe(N, std::vector<double>(N,0));
  std::vector<double> d(N,0);

  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++)
      ARe[i][j] = A[i][j].real();
  
  // Diagonalize the matrix. ARe will contain O after completion
  bool fail = eigensystem(ARe, d);
  if (fail)
  {
    std::cout << "*** Error in global.cpp:make_posdef_ReSymm ***\n";
    std::cout << "There was some problem with the eigenvalue solver.\n";
    std::cout << "Input matrix:" << std::endl;
    for (int i=0; i<N; i++)
      for (int j=0; j<N; j++)
        ARe[i][j] = A[i][j].real();
    MatrixScreenOutput(ARe);

    exit(1);
  }

  // Eliminate the negative eigenvalues
  for (int i=0; i<N; i++)
  {
#ifdef DEBUG
    if(d[i] < 0.) std::cout << "ELIMINATING A MODE" << std::endl;
#endif
    d[i] = std::max(d[i], 0.);
  }

  // Reconstruct the matrix with those eigenvalues eliminated
  for(int i=0; i<N; i++)
  {
    for(int j=0; j<N; j++)
    {
      A[i][j] = FieldType(0.);

      for(int k=0; k<N; k++)
      {
        A[i][j] += ARe[i][k] * ARe[j][k] * d[k];
      }
    }
  }

} // }}}


RealType linterp( RealType x_star, const Vec1dReal & x_data, const Vec1dReal & y_data )
{ // {{{

  // linear interpolation, 
  // assumptions: 
  // * equally spaced x-data
  // * x_data is monotonically increasing

  RealType dx = x_data[1] - x_data[0];
  int idx = floor( (x_star - x_data[0])/dx );
  RealType y_star;

  if (idx == x_data.size()-1)
  {
    y_star = y_data[idx];
  }
  else
  {
    y_star = (y_data[idx+1] - y_data[idx])*(x_star - x_data[idx])/dx + y_data[idx];
  }

  return y_star;

} // }}}

// --- End DRT edits ---

// --- JUG edits ---
bool jsonFile(const char* fname)
{//function returns true if fname is a json file
 //criteria: search for #, signature  of legacy input files
 //if not found, default to json 

  std::ifstream f1 (fname, std::ifstream::in);
  if ( f1.is_open() )
  {
    std::string tmp_str;
    std::getline(f1, tmp_str, '#');
    if(!f1.eof())
    {
      f1.close();
      std::cout << "\n";
      std::cout << " * Legacy format (07/2020) detected for " << fname << std::endl;
      std::cout << " * Format discontinued; use JSON input files for future runs." << std::endl;
      std::cout << "\n";
      return false;
    }
    else
    {
      f1.close();
      return true;
    }
  }
  else
  {
    std::cout << "(Error): Cannot open file, " << fname << std::endl;
    exit(1);
  }
  
}
// --- end: JUG edits ---

// Handle memory error: print location and exit code
void memerr(const char* file, int line)
{
  std::cerr << std::endl << std::endl << "  ***  Error allocating memory at line " << line << " in " << file << std::endl;
  exit(1);
}

// Abort due to code error or inconsistency
void codeerror_abort(const char *message, const char* file, int line)
{
  std::cerr << std::endl;
  std::cerr << " ERROR: code error or inconsistency" << std::endl;
  std::cerr << "  - Message = " << message << std::endl;
  std::cerr << "  - Received from line " << line << " in source file " << file << std::endl;
  abort();
}

// Exit due to unreasonable user behavior
void usererror_exit(const char *message, const char* file, int line)
{
  std::cerr << std::endl;
  std::cerr << " USAGE ERROR:" << std::endl;
  std::cerr << "  - Message = " << message << std::endl;
  std::cerr << "  - Received from line " << line << " in source file " << file << std::endl;
  exit(2); // Throw a non-zero return code
}

// Utility function for fast determination of power of 2 numbers
// used, for example, for fixing CUDA block sizes
unsigned int nextpow2(unsigned int v)
{
  // This algorithm works by propagating right-most set bit to all lower
  // bits, then adding one which carries the all bits to zero except
  // rightmost+1

  // ASSUMPTION: 'int' is 32 bits wide
  assert(sizeof(int)==4); // This check is disabled for release mode to improve speed

  v--;
  v |= v >> 1;
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
  v++;
  return(v);
}

// Utility function for fast determination of power of 2 numbers
unsigned int prevpow2(unsigned int v)
{
  // This algorithm works by propagating right-most set bit to all lower
  // bits, then adding one which carries the all bits to zero except
  // rightmost+1

  // ASSUMPTION: 'int' is 32 bits wide
  assert(sizeof(int)==4); // This check is disabled for release mode to improve speed
  // If already a power of 2, return original
  if((v != 0) && !(v & (v - 1)))
    return(v);
  else
    // Use nextpow2 to round to the next higher power of 2, then logic shift by 1
    return(nextpow2(v)>>1);
}

// Utility function for fast determination of multiples of 2
unsigned int nextmult2(unsigned int v)
{
  // Using integer rounding, we want v = (v+1)/2 * 2
  // However, i/2*2 can be replaced with (i & ~0x1)
  // Similarly, i/4*4 could be replaced with (i & ~0x3)
  unsigned int mask(0x1);
  mask = ~mask;
  return((v + 1) & mask);
}


#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define M 7
#define NSTACK 50
// Numerical Recipes in C, Sec. 8.4 - Quicksort Routine
void indexx(UInt n, const RealType arr[], UInt indx[])
// Indexes an array arr[1..n], i.e., outputs the array indx[1..n] such that arr[indx[j]] is
// in ascending order for j =1, 2,. ..,N. The input quantities n and arr are not changed.
{
  UInt i, indxt, ir=n, itemp, j, k, l=1;
  int jstack=0;
  RealType a;
  // Set up stack
  std::vector<int> istack(NSTACK);
  // Initialize indx to input order
  for(j=1;j<=n;j++)
    indx[j]=j-1;

  for(;;)
  {
    if(ir-l < M)
    {
      for(j=l+1;j<=ir;j++)
      {
        indxt=indx[j];
        a=arr[indxt];
        for(i=j-1;i>=l;i--)
        {
          if(arr[indx[i]] <= a)
            break;
          indx[i+1]=indx[i];
        }
        indx[i+1]=indxt;
      }
      if(jstack == 0)
        break;
      ir=istack[jstack--];
      l=istack[jstack--];
    }
    else
    {
      k=(l+ir) >> 1;
      SWAP(indx[k],indx[l+1]);
      if(arr[indx[l]] > arr[indx[ir]])
      {
        SWAP(indx[l],indx[ir])
      }
      if(arr[indx[l+1]] > arr[indx[ir]])
      {
        SWAP(indx[l+1],indx[ir])
      }
      if(arr[indx[l]] > arr[indx[l+1]])
      {
        SWAP(indx[l],indx[l+1])
      }
      i=l+1;
      j=ir;
      indxt=indx[l+1];
      a=arr[indxt];
      for (;;)
      {
        do i++; while(arr[indx[i]] < a);
        do j--; while(arr[indx[j]] > a);
        if (j < i)
          break;
        SWAP(indx[i],indx[j])
      }
      indx[l+1]=indx[j];
      indx[j]=indxt;
      jstack += 2;
      if(jstack > NSTACK)
        codeerror_abort("NSTACK too small in indexx.",__FILE__,__LINE__);
      if(ir-i+1 >= j-l)
      {
        istack[jstack]=ir;
        istack[jstack-1]=i;
        ir=j-1;
      }
      else
      {
        istack[jstack]=j-1;
        istack[jstack-1]=l;
        l=i;
      }
    }
  }
}

void writeStatus(enum enum_runstatus RunStatus)
{
  FILE *statfile=fopen("STATUS","w");
  fprintf(statfile,"%d",static_cast<int>(RunStatus));
  fclose(statfile);
}

bool checkStopStatus()
{
  FILE *status=fopen("STOP","r");
  if(status==NULL)
    return false;
  else
  {
    fclose(status);
    return true;
  }
}

void clearStopStatus()
{
  remove("STOP");
}

namespace iobase{
  // Create an instance of the nullstream class that will reside in this namespace for the process duration
#ifdef DEBUG
  dbgstreambuf dbgstream(std::cout);
  std::ostream cdbg(&dbgstream); // If debugging is on, send cdbg to cout with a stamp at the start of each line.
#else
  nullstreambuf nullstream;
  std::ostream cdbg(&nullstream); // If debugging is not on, send cdbg to oblivion
#endif
}
