/**********************************************************************************
 *
 * File adaped from PolyFTS Project
 * Originally created by Kris Delaney on 2011-04-16.
 * Edited by Doug Tree on 2016-05-23.
 * Copyright (c) 2011,2016 University of California, Santa Barbara. All rights reserved.
 *
 **********************************************************************************/
#ifndef _GLOBAL_H_
#define _GLOBAL_H_

// A global include file - will be included everywhere in the code.

// Some extra debug output and checks are enabled if DEBUG is defined
// If it is NOT defined, we defined NDEBUG (no debug) which additionally
// optimizes away all assert() checks in the code - for higher performance
#ifndef DEBUG
#define NDEBUG
#endif

/*---------------------------------------------------------*
 *  Define the default datatypes used throughout the code. *
 *  Allows simple precision changing without rewrite       *
 *---------------------------------------------------------*/
#include <complex>
#include <vector>
#include <list>
#include <map>
#include "rapidjson/document.h"
#include "rapidjson/istreamwrapper.h"

// Also fetch std integer types
#include <inttypes.h>
typedef uint64_t UInt; // unsigned integer used throughout for iterating and sizing objects
#define IOFLOATDIGITS 10

// Typedefs that very with global precision settings
// Global precision is the prec. used for input/output and collecting final
// results. Actual intensive computation has precision set at runtime.
typedef double RealType;
// General typedefs - not dependent on precision
typedef std::complex<RealType> Complex;
typedef std::vector<RealType> Vec1dReal;
typedef std::vector<std::vector<RealType> > Vec2dReal;
typedef std::vector<Complex> Vec1dCmplx;
typedef std::vector<std::vector<Complex> > Vec2dCmplx;
// --- DRT edits ---
typedef std::complex<RealType> FieldType; // can be RealType or Complex
typedef std::vector<FieldType> Vec1dFieldType; // can be RealType or Complex
typedef std::vector< std::vector<FieldType> > Vec2dFieldType; // can be RealType or Complex
typedef std::vector<int> Vec1dInt; 
typedef std::vector< std::vector<int> > Vec2dInt;
typedef std::vector<UInt> Vec1dUInt;
typedef std::vector< std::vector<UInt> > Vec2dUInt;
typedef std::vector<Vec2dReal> Vec3dReal;
//#define lapack_complex_float std::complex<float>
//#define lapack_complex_double std::complex<double>
//#ifdef __MKL__
//#include <mkl.h>
//#else
//#include <lapacke.h>
//#include <cblas.h>
//#endif
// --- End DRT edits ---
// Utility includes for stream handling, floating comparisons and timers
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <stdlib.h>
#include <cerrno>

#include <math.h>   // Math operations
#include <cmath>
#include <assert.h> // Debugging using assertions
#include <limits>
#include <algorithm>

// Timers & other system-specific
#include <ctime>
#include "time.h"      // For ctime()
#include <unistd.h>    // For gethostname
#include <sys/param.h> // For gethostname
#include <sys/types.h> // For gethostname
#include <sys/time.h>
#include <sys/times.h>


// FFTW
#ifndef __GPU__
#include <fftw3.h>
#ifdef __MPI__
#include <fftw3-mpi.h>
#endif
#endif

// MPI parallel
#ifdef __MPI__
#include <mpi.h>
#endif

// OMP Threading
#ifdef __OMP__
#include <omp.h>
#endif

// BLAS routines
#ifdef BLAS
#include "cblas.h"
#endif

// --- DRT edits ---
//#include "lapackpp.h"
// --- end DRT edits ---

// Literal constants
#define PI 3.14159265358979323846
#define TPI 6.28318530717958647693

enum enum_runstatus {RUNNING=0, FAILED=1, DONE=2, TIMEOUT=3};

// Forward declarations of global data and functions.
bool CompareEqual(double A, double B, double tol=-1.0); // Compare two floating point numbers for equality within machine precision
bool CompareEqual(std::complex<double> A, double B, double tol=-1.0); // Compare two floating point numbers for equality within machine precision
bool CompareEqual(double A, std::complex<double> B, double tol=-1.0); // Compare two floating point numbers for equality within machine precision
bool CompareEqual(std::complex<double> A, std::complex<double> B, double tol=-1.0); // Compare two floating point numbers for equality within machine precision
bool CompareEqual(float A, float B, float tol=-1.0f); // Compare two floating point numbers for equality within machine precision
bool CompareEqual(std::complex<float> A, float B, float tol=-1.0f); // Compare two floating point numbers for equality within machine precision
bool CompareEqual(float A, std::complex<float> B, float tol=-1.0f); // Compare two floating point numbers for equality within machine precision
bool CompareEqual(std::complex<float> A, std::complex<float> B, float tol=-1.0f); // Compare two floating point numbers for equality within machine precision
bool CompareEqual(UInt A, UInt B, float tol=-1.0f); // Compare two floating point numbers for equality within machine precision

// --- DRT edits ---
Vec2dFieldType matmul( Vec2dFieldType A, Vec2dFieldType B);
Vec2dFieldType eye( UInt Nrow, UInt Ncol, FieldType scalar );
Vec2dFieldType kroneckerdelta( Vec1dReal );
//void make_posdef( Vec2dFieldType &A );
void make_posdef_ReSymm( Vec2dFieldType &A );
RealType linterp( RealType x_star, const Vec1dReal & x_data, const Vec1dReal & y_data );

// --- End DRT edits ---
bool jsonFile(const char* fname);

#ifdef __cplusplus
extern "C" {
#endif
  void memerr(const char* file, int line);      // Report location of a memory allocation error
  void codeerror_abort(const char* message, const char* file, int line);      // abort due to some code inconsistency
  void usererror_exit(const char* message, const char* file, int line);      // exit due to user input
  unsigned int nextpow2(unsigned int in); // Round up to next integer power of 2
  unsigned int prevpow2(unsigned int in); // Round down to next integer power of 2
  unsigned int nextmult2(unsigned int in); // Round up to next multiple of 2
  double gettime(); // Get the current time in seconds
  void indexx(UInt n, const RealType arr[], UInt indx[]);
  void writeStatus(enum enum_runstatus RunStatus);
  bool checkStopStatus();
  void clearStopStatus();
#ifdef __cplusplus
}
#endif

namespace iobase {
  // Create a nullstream class that discards all content sent to it
  class nullstreambuf : public std::streambuf
  {
    protected:
      virtual int overflow(int c) {return 0;};
  };
  // Create a debug stream that writes messages with "!DEBUG " at
  // the start of each line and sends the remaining part of the message to the
  // ostream object that was passed in the ctor arg.
  class dbgstreambuf : public std::streambuf
  {
    public:
      dbgstreambuf(std::ostream &out) : _out(out), _pSink(0), _newline(true)
      {
        _pSink = _out.rdbuf(); // Take a copy of the output's streambuf
      };
    protected:
      virtual int_type overflow(int_type c = traits_type::eof())
      {
        if(traits_type::eq_int_type(c, traits_type::eof()))
          return _pSink->pubsync() == -1 ? c : traits_type::not_eof(c); // Convention: sync buffer if EOF
        if(_newline)
        {
          std::ostream str(_pSink);
          if(!(str<<"!DEBUG "))
            return traits_type::eof();
        }
        _newline = traits_type::to_char_type(c) == '\n';
        return _pSink->sputc(c); // Send char to output
      };
    private:
      std::ostream &_out;
      std::streambuf *_pSink;
      bool _newline;
  };
  // Specify that the cdbg object exists in this namespace scope
  // to any translation unit that includes global.h
  extern std::ostream cdbg;
}

#endif
