/**********************************************************************************
 *
 * PolyFTS Project
 *
 * File created by Kris Delaney on 2012-05-29.
 * Copyright (c) 2012 University of California, Santa Barbara. All rights reserved.
 *
 **********************************************************************************/
#ifndef _MATRIX_UTILS_H_
#define _MATRIX_UTILS_H_
#include "global.h"
#include <math.h>

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

// Interface routines
template<typename T> bool eigensystem(std::vector<std::vector<T> > &A, std::vector<T> &d);
template<typename T> bool invert(std::vector<std::vector<T> > &A);
template<typename T> bool solvelinsystem(std::vector<std::vector<T> > &A, std::vector<T> &b);

// Numerical Recipes backend utilities
template<typename T> void tred2(std::vector<std::vector<T> > &a, int n, std::vector<T> &d, std::vector<T> &e);
template<typename T> void tqli(std::vector<T> &d, std::vector<T> &e, int n, std::vector<std::vector<T> > &z);
template<typename T> T pythag(T a, T b);
template<typename T> void ludcmp(std::vector<std::vector<T> > &a, UInt n, std::vector<UInt> &indx, T &d);
template<typename T> void lubksb(std::vector<std::vector<T> > &a, UInt n, std::vector<UInt> &indx, std::vector<T> &b);
template<typename T> void MatrixScreenOutput(std::vector<std::vector<T> > const &M);

// Additional functions (JUG)
template<typename T> void VectorScreenOutput(std::vector<T>  const &v);
template<typename T> void rootFind2 (T const &b, T const &c, T &r0, T &r1);
template<typename T> void eigSys2 (std::vector<std::vector<T> > &A, std::vector<T> &d);
template<typename T> bool eigTrue (std::vector<std::vector<T> > const &A, std::vector<T> const &v, T const &l);

#endif
