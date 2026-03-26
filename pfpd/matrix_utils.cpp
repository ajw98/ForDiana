/**********************************************************************************
 *
 * PolyFTS Project
 *
 * File created by Kris Delaney on 2012-05-29.
 * Copyright (c) 2012 University of California, Santa Barbara. All rights reserved.
 *
 **********************************************************************************/
#include "matrix_utils.h"
#include <iostream>

#define TINY 1.0e-20 // A small number

//#define DEBUGMATUTILS


template<typename T>
bool eigensystem(std::vector<std::vector<T> > &A, std::vector<T> &d)
{//Solves the eigenvalue problem A*l = l*v.
 //Changes A to contain its eigenvectors (v) as column vectors. 
 //Changes d to contain the eigenvalues(l). 
 //Integrity of input values are NOT preserved.
 //For a real symmetric matrix A: A=P*DP; note: P* = P-transpose 
 //D is a diagonal matrix with the eigenvalues (l) of A.
 //P columns are the eigenvectors (v) of A.
 //Returns true if any eigenvalue is zero or is not finite (INF, NAN).

  UInt n=A.size();
  
  if (n == 2)//direct calculation
  {
    eigSys2 (A, d);
  }
  else //general solver from numerical recipes
  {
    std::vector<T> e(n);
    tred2(A, n, d, e);
    tqli(d, e, n, A);
  }
  //singularity and validity check
  for(UInt i=0; i<n; i++)
    if(!std::isnormal(d[i]))
      return true;
  return false;
}

template<typename T>
bool invert(std::vector<std::vector<T> > &A)
{
  UInt n=A.size();
  if(n<1)
    return false;
  UInt np=A[0].size();
  if(n != np)
  {
    std::cout << " Inversion failed: non-square matrix" << std::endl;
    return false;
  }

#ifdef DEBUGMATUTILS
  std::cout << "Matrix inversion using LU decomposition and backsubstitution" << std::endl;
  std::cout << " - Input matrix: " << std::endl;
  for(UInt i=0; i<n; i++)
  {
    std::cout << "(";
    for(UInt j=0; j<n; j++)
      std::cout <<std::setw(9) << std::right << A[i][j] << " ";
    std::cout << ")"<<std::endl;
  }
  std::cout << std::endl;
  std::vector<std::vector<T> > Acopy=A;
#endif

  // LU decomposition
  T d; // Parity from row switching
  std::vector<UInt> indx(n,0);
  ludcmp(A,n,indx,d);
  if(CompareEqual(d,T(0.0))) // By convention, ludcmp has been modified to return d=0 for singular matrices
  {
    std::cout << " Inversion failed: singular matrix" << std::endl;
    return false;
  }

  // Inversion by backsubstitution
  std::vector<T> col(n);
  std::vector<std::vector<T> > Ainv=A;
  for(UInt j=0; j<n; j++)
  {
    for(UInt i=0; i<n; i++)
      col[i] = 0.0;
    col[j] = 1.0;
    lubksb(A,n,indx,col);
    for(UInt i=0; i<n; i++)
      Ainv[i][j] = col[i];
  }
  // A is destroyed anyway in LUdcmp; replace with inverse
  A=Ainv;

#ifdef DEBUGMATUTILS
  std::cout << " - Inverse matrix: " << std::endl;
  for(UInt i=0; i<n; i++)
  {
    std::cout << "(";
    for(UInt j=0; j<n; j++)
      std::cout <<std::setw(9) << std::right << A[i][j] << " ";
    std::cout << ")"<<std::endl;
  }
  std::cout << std::endl;
  // Test that the product of the original and inverse matrices is the identity; dump result into Ainv for this debug test
  for(UInt i=0; i<n; i++)
    for(UInt j=0; j<n; j++)
    {
      Ainv[i][j] = 0.0;
      for(UInt k=0; k<n; k++)
        Ainv[i][j] += A[i][k] * Acopy[k][j];
//        Ainv[i][j] += A[k][j] * Acopy[i][k]; // Can test either mult order
    }
  // Output result
  std::cout << " - Product of original and inverse: " << std::endl;
  for(UInt i=0; i<n; i++)
  {
    std::cout << "(";
    for(UInt j=0; j<n; j++)
      std::cout <<std::setw(9) << std::right << Ainv[i][j] << " ";
    std::cout << ")"<<std::endl;
  }
  std::cout << std::endl;

#endif

  return true; // success
}

// Solve for x in A.x=b
// 'b' is replaced with 'x', and "A" is replaced with its packed LU decomposition (routine destroys integrity of both inputs)
template<typename T>
bool solvelinsystem(std::vector<std::vector<T> > &A, std::vector<T> &b)
{
  UInt n=A.size();
  if(n<1)
    return false;
  UInt np=A[0].size();
  if(n != np)
  {
    std::cout << " Linear solver failed: non-square matrix" << std::endl;
    return false;
  }
  if(b.size() != n)
  {
    std::cout << " Linear solver failed: input 'b' is wrong size" << std::endl;
    return false;
  }

#ifdef DEBUGMATUTILS
  std::cout << "Solve linear system A.x=b for x" << std::endl;
  std::cout << " - Input matrix: " << std::endl;
  for(UInt i=0; i<n; i++)
  {
    std::cout << "    (";
    for(UInt j=0; j<n; j++)
      std::cout <<std::setw(9) << std::right << A[i][j] << " ";
    std::cout << ")"<<std::endl;
  }
  std::cout << std::endl;
  std::cout << " - Input vector: " << std::endl;
  std::cout << "(";
  for(UInt j=0; j<n; j++)
    std::cout <<std::setw(9) << std::right << b[j] << " ";
  std::cout << ")"<<std::endl;
  std::vector<std::vector<T> > Acopy=A;
  std::vector<T> bcopy=b;
#endif

  std::vector<UInt> indx(n);
  T d; // Parity
  ludcmp(A,n,indx,d);
  if(CompareEqual(d,T(0.0)))
  {
    std::cout << "Linear system solver: singular matrix" << std::endl;
    return false;
  }
  lubksb(A,n,indx,b);

#ifdef DEBUGMATUTILS
  std::cout << " Test solution: " << std::endl;
  std::vector<T> Ax(n,0.0);
  for(UInt i=0; i<n; i++)
    for(UInt k=0; k<n; k++)
      Ax[i] += Acopy[i][k] * b[k];
  std::cout << " - A.x: " << std::endl;
  std::cout << "    (";
  for(UInt j=0; j<n; j++)
    std::cout <<std::setw(9) << std::right << Ax[j] << " ";
  std::cout << ")"<<std::endl;
  std::cout << " - b (orginal): " << std::endl;
  std::cout << "    (";
  for(UInt j=0; j<n; j++)
    std::cout <<std::setw(9) << std::right << bcopy[j] << " ";
  std::cout << ")"<<std::endl;
#endif


  return true;
}


// Householder reduction of a real, symmetric matrix to tridiagonal form.
// From Numerical Recipes in C. The routine expects array indices to run [1..n]
// whereas ours run from [0..n-1]. Hence we shift indices when accessing arrays
// without otherwise altering the code.
//
// This routine is destructive. Upon completion, the matrix a is replaced by
// the orthogonal matrix that brings the original a into tridiagonal form.
// d is the vector of diagonal entries of the reduced form, and e is the off diagaonals (with e[0] = 0)
template<typename T>
void tred2(std::vector<std::vector<T> > &a, int n, std::vector<T> &d, std::vector<T> &e)
{
  int l, k, j, i;
  T scale, hh, h, g, f;

  for( i=n ; i>=2 ; i-- )
  {
    l = i-1;
    h=scale=0.0;

    if(l > 1)
    {
      for(k=1 ; k<=l ; k++)
      {
        scale += fabs(a[i-1][k-1]);
      }
      if(scale == 0.0)
      {
        e[i-1] = a[i-1][l-1];
      }
      else // scale != 0
      {
        for(k=1 ; k<=l ; k++)
        {
          a[i-1][k-1] /= scale;
          h += a[i-1][k-1]*a[i-1][k-1];
        }
        f=a[i-1][l-1];
        g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
        e[i-1]=scale*g;
        h -= f*g;
        a[i-1][l-1]=f-g;
        f=0.0;
        for(j=1;j<=l;j++)
        {
          a[j-1][i-1]=a[i-1][j-1]/h;
          g=0.0;
          for(k=1;k<=j;k++)
            g += a[j-1][k-1]*a[i-1][k-1];
          for(k=j+1;k<=l;k++)
            g += a[k-1][j-1]*a[i-1][k-1];
          e[j-1]=g/h;
          f += e[j-1]*a[i-1][j-1];
        }
        hh=f/(h+h);
        for(j=1;j<=l;j++)
        {
          f=a[i-1][j-1];
          e[j-1]=g=e[j-1]-hh*f;
          for(k=1;k<=j;k++)
            a[j-1][k-1] -= (f*e[k-1]+g*a[i-1][k-1]);
        }
      }
    }
    else // if l == 0
    {
      e[i-1]=a[i-1][l-1];
    }
    d[i-1]=h;
  } // i

  d[0]=0.0;
  e[0]=0.0;

  for(i=1;i<=n;i++)
  {
    l = i-1;
    if(d[i-1] != 0.0)
    {
      for(j=1;j<=l;j++)
      {
        g=0.0;
        for(k=1;k<=l;k++)
        {
          g += a[i-1][k-1]*a[k-1][j-1];
        }
        for(k=1;k<=l;k++)
        {
          a[k-1][j-1] -= g*a[k-1][i-1];
        }
      }
    }
    d[i-1]=a[i-1][i-1];
    a[i-1][i-1]=1.0;
    for(j=1;j<=l;j++)
    {
      a[j-1][i-1] = a[i-1][j-1] = 0.0;
    }
  }

}
/* (C) Copr. 1986-92 Numerical Recipes Software )k45D,&4. */
// As for the routine above, we take the Num. Rec. in C version of
// this code directly, which assumes arrays are indexed 1..n. Rather
// than modify loop counters and flow of the algorithm to account for our
// zero-based indexing, appropriately, we simply shift all array accesses by -1
template<typename T>
void tqli(std::vector<T> &d, std::vector<T> &e, int n, std::vector<std::vector<T> > &z)
{
  int m,l,iter,i,k;
  T s,r,p,g,f,dd,c,b;

  for(i=2;i<=n;i++)
    e[i-2]=e[i-1];
  e[n-1]=0.0;

  for(l=1;l<=n;l++)
  {
    iter=0;
    do
    {
      for(m=l;m<=n-1;m++)
      {
        dd=fabs(d[m-1])+fabs(d[m]);
        if ((double)(fabs(e[m-1])+dd) == dd)
          break;
      }
      if(m != l)
      {
        if (iter++ == 30)
	 {
          std::cerr << "Too many iterations in tqli" << std::endl;
          if (!std::isnormal(d[m]))
	  {
            std::cout << " *** Error: Eigenvalue calculation failed for a ";
	    std::cout << n << "x" << n << " system." << std::endl;
	    std::cout << " *** Calculated eigenvalue = " << d[m] << std::endl;
	    exit(1);
	  } 
	 }
	g=(d[l]-d[l-1])/(2.0*e[l-1]);
        r=pythag(g,T(1.0));
        g=d[m-1]-d[l-1]+e[l-1]/(g+SIGN(r,g));
        s=c=1.0;
        p=0.0;
        for(i=m-1 ; i >= l ; i--)
        {
          f=s*e[i-1];
          b=c*e[i-1];
          e[i]=(r=pythag(f,g));
          if(r == 0.0) {
            d[i] -= p;
            e[m-1]=0.0;
            break;
          }
          s=f/r;
          c=g/r;
          g=d[i]-p;
          r=(d[i-1]-g)*s+2.0*c*b;
          d[i]=g+(p=s*r);
          g=c*r-b;
          for(k=1;k<=n;k++)
          {
            f=z[k-1][i];
            z[k-1][i]=s*z[k-1][i-1]+c*f;
            z[k-1][i-1]=c*z[k-1][i-1]-s*f;
          }
        }
        if(r == 0.0 && i >= l)
          continue;
        d[l-1] -= p;
        e[l-1]=g;
        e[m-1]=0.0;
      } // if(m!=l)
    } while (m != l);
  }
}
/* (C) Copr. 1986-92 Numerical Recipes Software )k45D,&4. */

template<typename T>
T pythag(T a, T b)
{//Return c where c*c = a*a + b*b
  
  T absa,absb;
  absa=fabs(a);
  absb=fabs(b);

  if (absa > absb) return absa*sqrt(1.0+(absb*absb)/(absa*absa));
  else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+(absa*absa)/(absb*absb)));

}
/* (C) Copr. 1986-92 Numerical Recipes Software )k45D,&4. */


// LU decomposition from Numerical Recipes in C, 2nd Edition
// Translated to use STL containers for matrices and vectors
// Input: square matrix, a, of size, n
// Output: a is replaces by its LU decomposition (arranged as in Num. Rec. in C 2nd Ed. Eq. 2.3.14)
//         d - row permutation parity
//         indx - vector of row permutations imposed by pivoting procedure
template<typename T>
void ludcmp(std::vector<std::vector<T> > &a, UInt n, std::vector<UInt> &indx, T &d)
{
  UInt i,imax,j,k;
  T big,dum,sum,temp;
  std::vector<T> vv(n); // vv stores the implicit scaling of each row
  d=1.0; // No row interchanges as yet => parity +1
  for (i=0;i<n;i++)  // Loop over rows to get the implicit scaling information
  {
    big=0.0;
    for (j=0;j<n;j++)
      if ((temp=fabs(a[i][j])) > big)
        big=temp;

    if (CompareEqual(big,T(0.0)))
    {
      std::cout << "Singular matrix in routine ludcmp" << std::endl; // No nonzero largest element.
      // MODIFY...
      d=0.0; // Signature to caller that matrix is singular.
      return;
    }

    vv[i]=1.0/big;  // Save the scaling
  }
  for (j=0;j<n;j++) // This is the loop over columns of Crout’s method.
  {
    for (i=0;i<j;i++) // This is equation (2.3.12) except for i = j.
    {
      sum=a[i][j];
      // Save the scaling.
      for (k=0;k<i;k++)
        sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    big=0.0;
    for (i=j;i<n;i++)
    {
      sum=a[i][j];
      // Initialize for the search for largest pivot element.
      // This is i = j of equation (2.3.12) and i = j+1...N
      // of equation (2.3.13).
      for (k=0;k<j;k++)
        sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
      if ((dum=vv[i]*fabs(sum)) >= big)
      {
        // Is the figure of merit for the pivot better than the best so far?
        big=dum;
        imax=i;
      }
    }
    // Do we need to interchange rows?
    if (j != imax)
    {
      // Yes, do so...
      for (k=0;k<n;k++)
      {
        dum=a[imax][k];
        a[imax][k]=a[j][k];
        a[j][k]=dum;
      }
      // ...and change the parity of d. Also interchange the scale factor.
      d = -d;
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    // If the pivot element is zero the matrix is singular (at least to the precision of the algorithm).
    // For some applications on singular matrices, it is desirable to substitute TINY for zero.
    if (CompareEqual(a[j][j],T(0.0)))
    {
      a[j][j]=TINY;
      // MODIFY...
      d=0.0; // Inform the calling routine that the matrix is singular
    }
    if (j != n-1)
    {
      // Now, finally, divide by the pivot element.
      dum=1.0/(a[j][j]);
      for (i=j+1;i<n;i++)
        a[i][j] *= dum;
    }
  }  // Go back for the next column in the reduction.
}

// Num. Recipes in C, 2nd Ed.
// Forward/backward substitution, as in Eq. 2.3.6 and 2.3.7
// - Solves the set of n linear equations A·X = B. Here a[0..n-1][0..n-1] is input,
//   not as the matrix A but rather as its LU decomposition, determined by the routine ludcmp.
//   indx[0..n-1] is input as the permutation vector returned by ludcmp.
//   b[0..n-1] is input as the right-hand side vector B, and returns with the solution vector X.
//   a, n, and indx are not modified by this routine and can be left in place for successive
//   calls with different right-hand sides b. This routine takes into account the possibility
//   that b will begin with many zero elements, so it is efficient for use in matrix inversion.
template<typename T>
void lubksb(std::vector<std::vector<T> > &a, UInt n, std::vector<UInt> &indx, std::vector<T> &b)
{
  int ii(-1);
  T sum;
  for (UInt i=0;i<n;i++)
  {
    // When ii is set to a nonnegative value, it will become the index of
    // the first nonvanishing element of b. We now do the forward substitution,
    // equation (2.3.6). The only new wrinkle is to unscramble the permutation as we go.
    UInt ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    if (ii>=0)
      for (UInt j=static_cast<UInt>(ii);j<=i-1;j++)
        sum -= a[i][j]*b[j];
    else if (!CompareEqual(sum,T(0.0))) // A nonzero element was encountered, so from now on we will have to do the sums in the loop above.
      ii=static_cast<int>(i);
    b[i]=sum;
  }

  UInt i=n;
  for (;i-- > 0;) // i=[n-1,0]
  {
    // Now we do the backsubstitution, equation (2.3.7).
    sum=b[i];
    for (UInt j=i+1;j<n;j++)
      sum -= a[i][j]*b[j];
    b[i]=sum/a[i][i]; // Store a component of the solution vector X.
  }
}


template<typename T>
void MatrixScreenOutput(std::vector<std::vector<T> > const &M)
{
  std::cout << std::showpos;
  for(UInt i=0; i<M.size(); i++)
  {
    std::cout << "               ( ";
    for(UInt j=0; j<M[0].size(); j++)
      std::cout << std::setw(8) << std::right << M[i][j] << " ";
    std::cout << ")" << std::endl;
  }
  std::cout << std::noshowpos;
  std::cout << std::endl;
}

template<typename T>
void VectorScreenOutput(std::vector<T>  const &v)
{
  std::cout << std::showpos;
  for(UInt i=0; i<v.size(); i++)
  {
    std::cout << "               ( ";
    std::cout << std::setw(8) << std::right << v[i] << " ";
    std::cout << ")" << std::endl;
  }
  std::cout << std::noshowpos;
  std::cout << std::endl;
}

template<typename T>
void rootFind2 (T const &b, T const &c, T &r0, T &r1)
{//returns the 2 roots, r0 and r1, of 'x**2 + bx + c = 0'
 //function works for Real and Imaginary roots as long as complex.h is included.
 
  //sqrt(b*b-4c)
  T det = b;
  T tmp = c;
  tmp *= 4.;
  det *= det;
  det -= tmp;
  det = sqrt(det); 

  //r0 = (-b - det)/2 
  r0 = -b;
  r0 -= det;
  r0 /= 2.;

  //r1 = (-b + det)/2 
  r1 = -b;
  r1 += det;
  r1 /= 2.;
}

template<typename T>
void eigSys2 (std::vector<std::vector<T> > &A, std::vector<T> &d)
{//Solves the eigenvalue problem A*l = l*v for a Real 2x2 matrix 
 //Changes A to contain its eigenvectors (v) as column vectors. 
 //Changes d to contain the eigenvalues(l). 
 //Does not require symmetry to work, but symmetric matrices have nice properties.
 //For a real symmetric matrix A: A=P*DP; note: P* = P-transpose 
 //D is a diagonal matrix with the eigenvalues (l) of A.
 //P columns are the eigenvectors (v) of A.
 
  if (A.size() != 2 || A[0].size() != 2)
  {
    std::cout << "***Error: Called eigSys2 for a non-2x2 matrix***" << std::endl;
    exit(1); 
  }
  
  //For testing--paired with Validate solution
  //std::vector<std::vector<T> > Aorig = A; 
 
  //Set up the characteristic equation
  T trc = A[0][0]; trc += A[1][1]; 
  T det = A[0][0]; det *= A[1][1];  
  T tmp = A[0][1]; tmp *= A[1][0];
  det -= tmp;
   
  //Solve the eigenvalues
  rootFind2 (-trc, det, d[0], d[1]);
  
  //Calculate eigenvectors using Cayley-Hamilton Theorem
  //https://en.wikipedia.org/wiki/Eigenvalue_algorithm
  
  //eigenvector for d[0]
  A[0][0] -= d[1];  
  //A[1][0] = A[1][0];//unchanged; left for completeness.  
  
  //eigenvector for d[1]
  //A[0][1] = A[0][1]; //unchanged; left for completeness.  
  A[1][1] -= d[0];  
  
  //alternate eigenvector for d[0]
  //A[0][0] = A[0][1];  
  //A[1][0] = A[1][1]-d[1];  
  //alternate eigenvector for d[1]
  //A[0][1] = A[0][0] - d[0];  
  //A[1][1] = A10;  

  //Normalize the column vectors (eigenvectors now saved in A) 
  for(UInt i=0; i<2; i++)
  {
    T norm = pythag (A[0][i], A[1][i]);
    A[0][i] /= norm;
    A[1][i] /= norm;
  }
  
  //Validate solution---commented out for performance 
  //for(UInt i=0; i<2; i++)
  //{
  // std::vector<T> v;
  // v.reserve(2);
  // v.emplace_back(A[0][i]);
  // v.emplace_back(A[1][i]);
  // if(!eigTrue(Aorig, v, d[i]))  
  //   exit(1);  
  //}
}

template<typename T>
bool eigTrue (std::vector<std::vector<T> > const &A, std::vector<T> const &v, T const &l)
{//Returns true if A*l = l*v
 //Useful for tests and debugging.	

 if (A[0].size() !=  v.size())
 {
   std::cout << "***Error: eigTrue inputs are not compatible.***" << std::endl;
   exit(1); 
 }
 
 for(UInt i=0; i < v.size(); i++)
 {
   T lhs = 0; 
   T rhs = 0; 
   for(UInt j=0; j < v.size(); j++)
   {  
     lhs += A[i][j]*v[j];
   }
   rhs = l*v[i];
   if ( fabs(lhs-rhs) > 1.e-10 )
   {
     std::cout << "Eigensystem test failed with lhs-rhs = " << fabs(lhs-rhs) << std::endl;
     std::cout << "A =" << std::endl;
     MatrixScreenOutput (A);
     std::cout << "v =" << std::endl;
     VectorScreenOutput (v); 
     std::cout << "l = " << l << std::endl; 
     std::cout << "Av = " << lhs << std::endl;
     std::cout << "lv = " << rhs << std::endl;
     return false;
   }
 }
 return true;
}


//template specialization stubs
template bool eigensystem(std::vector<std::vector<double> > &, std::vector<double> &);
template bool eigensystem(std::vector<std::vector<float> > &, std::vector<float> &);
template bool invert(std::vector<std::vector<double> > &);
template bool invert(std::vector<std::vector<float> > &);
template bool solvelinsystem(std::vector<std::vector<double> > &, std::vector<double> &);
template bool solvelinsystem(std::vector<std::vector<float> > &, std::vector<float> &);
template void tred2(std::vector<std::vector<double> > &, int, std::vector<double> &, std::vector<double> &);
template void tred2(std::vector<std::vector<float> > &, int, std::vector<float> &, std::vector<float> &);
template void tqli(std::vector<double> &, std::vector<double> &, int, std::vector<std::vector<double> > &);
template void tqli(std::vector<float> &, std::vector<float> &, int, std::vector<std::vector<float> > &);
template double pythag(double, double);
template float pythag(float, float);
template void ludcmp(std::vector<std::vector<double> > &, UInt, std::vector<UInt> &, double &);
template void ludcmp(std::vector<std::vector<float> > &, UInt, std::vector<UInt> &, float &);
template void lubksb(std::vector<std::vector<double> > &, UInt, std::vector<UInt> &, std::vector<double> &);
template void lubksb(std::vector<std::vector<float> > &, UInt, std::vector<UInt> &, std::vector<float> &);
template void MatrixScreenOutput(std::vector<std::vector<double> > const &);
template void MatrixScreenOutput(std::vector<std::vector<float> > const &);
template void VectorScreenOutput(std::vector<RealType>  const &);
template void VectorScreenOutput(std::vector<FieldType>  const &);
template void rootFind2 (RealType const &, RealType const &, RealType &, RealType &);
template void rootFind2 (FieldType const &, FieldType const &, FieldType &, FieldType &);
template void eigSys2 (std::vector<std::vector<RealType> > &, std::vector<RealType> &);
template bool eigTrue (std::vector<std::vector<RealType> > const &, std::vector<RealType> const &, RealType const &);
