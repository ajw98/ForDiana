/***************************************************************
*
* Class for Interleaved Banded Field Matrices
* 
* DRT -- Fri, 20 Mar 2017
*
****************************************************************/

#include "MatInterFieldMat.h"

// ----------------- Constructor/Destructor -----------------

MatInterFieldMat::MatInterFieldMat( int nrow, int ncol, int nfmrow, int nfmcol ) :
_Nrow(nrow),
_Ncol(ncol),
_NFMrow(nfmrow),
_NFMcol(nfmcol),
_LU( false ),
_Data( nrow, Vec1dIFM( ncol, InterFieldMat( nfmrow, nfmcol) ))
{ // {{{
  _BW = (_Nrow-1)/2;
  _NFD = _Ncol;
  _NIcomp = _NFMrow;
} // }}}

MatInterFieldMat::MatInterFieldMat( const MatInterFieldMat &rhs ):
_Nrow( rhs._Nrow ),
_Ncol( rhs._Ncol ),
_NFMrow( rhs._NFMrow ),
_NFMcol( rhs._NFMcol ),
_BW( rhs._BW ),
_NFD( rhs._NFD ),
_NIcomp( rhs._NIcomp ),
_LU( rhs._LU),
_Data( rhs._Data )
{ // {{{
} // }}}

MatInterFieldMat::~MatInterFieldMat()
{ // {{{
} // }}}

// ----------------- Operators  -----------------

MatInterFieldMat & MatInterFieldMat::operator=( const MatInterFieldMat & rhs )
{ // {{{

  this->_Nrow = rhs.Nrow();
  this->_Ncol = rhs.Ncol();
  this->_NFMrow = rhs.NFMrow();
  this->_NFMcol = rhs.NFMcol();
  this->_BW = rhs.BW();
  this->_NFD = rhs.NFD();
  this->_NIcomp = rhs.NIcomp();
  this->_LU = rhs.LU();
  this->_Data = rhs._Data;
  
  return *this;

} // }}}

// ----------------- Member Functions -----------------

void MatInterFieldMat::zero()
{ // {{{

  for (int i = 0; i < _Nrow; i++)
  {
    for (int j = 0; j < _Ncol; j++)
    {
      _Data[i][j] = FieldType(0.);
    }
  }
  _LU = false;

} // }}}

void MatInterFieldMat::setflag_inrealspace( bool realspaceflag )
{ // {{{

  for (int i = 0; i < _Nrow; i++)
  {
    for (int j = 0; j < _Ncol; j++)
    {
      _Data[i][j].setflag_inrealspace(realspaceflag);
    }
  }

} // }}}

bool MatInterFieldMat::getflag_inrealspace() const
{ // {{{

  bool rspace_flag = true;
  bool tmp_flag = true;

  for (int i = 0; i < _Nrow; i++)
  {
    for (int j = 0; j < _Ncol; j++)
    {
      tmp_flag = _Data[i][j].getflag_inrealspace();
      if (i > 0 and j > 0 and rspace_flag != tmp_flag)
      {
        std::cout << "Warning: SmartFieldVec is mixed real space/fourier space. Use component functions.";
      }
      rspace_flag = tmp_flag;
    }
  }

  return rspace_flag;

} // }}}

void MatInterFieldMat::smart_to_inter( SmartFieldOpMat & T )
{ // {{{

  if ( _BW != T.BW() or _Ncol != T.FDSize() or _NFMrow != T.Nrow() )
  {
    std::cout << "***(Error in MatInterFieldMat::interleave)***\n";
    std::cout << "Incompatible Matrix Sizes\n";
    exit(1);
  }

  SmartFieldMat tmp( T.Nrow(), T.Ncol() );
  tmp.zero();
  tmp.setflag_inrealspace( T.getflag_inrealspace() );

  this->zero();
  this->setflag_inrealspace( T.getflag_inrealspace() );

  for( int i = 0; i < _Nrow; i++ )
  {
    T.GetBand( i-_BW, tmp );

    for ( int j = 0; j < _NFMrow; j++ )
    {
      for ( int k = 0; k < _NFMcol; k++ )
      {
        for ( int l = 0; l < _Ncol; l++ )
        {
          _Data[i][l][j][k] = tmp(j,k)[l];
        }
      }
    }
  
  }

} // }}}

void MatInterFieldMat::inter_to_smart( SmartFieldOpMat & T )
{ // {{{

  T.zero();
  T.setflag_inrealspace( this->getflag_inrealspace() );

  SmartFieldMat tmp( _NFMrow, _NFMcol );
  tmp.zero();
  tmp.setflag_inrealspace( this->getflag_inrealspace() );

  for( int i = 0; i < _Nrow; i++ )
  {
    for ( int j = 0; j < _NFMrow; j++ )
    {
      for ( int k = 0; k < _NFMcol; k++ )
      {
        for ( int l = 0; l < _Ncol; l++ )
        {
          tmp(j,k)[l] = _Data[i][l][j][k];
        }
      }
    }
  
    T.AddBand( i-_BW, tmp );

  }

} // }}}

void MatInterFieldMat::band_LU()
{ // {{{

  // ************************************************************************
  //   B and X are allowed to be matrices so that one can do multiple
  //   solutions simultaneously for different rhs vectors. This is useful
  //   for the solution of periodic systems, where we need to solve multiple
  //   banded systems with the same matrix, but different rhs vectors.
  // ************************************************************************

  InterFieldMat tmpmat(_NIcomp, _NIcomp);
  InterFieldMat tmpmat2(_NIcomp, _NIcomp);
  InterFieldMat M(_NIcomp, _NIcomp); // elimination multiplier

  tmpmat.setflag_inrealspace(false);
  tmpmat2.setflag_inrealspace(false);
  M.setflag_inrealspace(false);
/*
  if (_LU == true)
  { 
    std::cout << "** WARNING: matrix has already been decomposed into L & U" << std::endl;
  }
*/
  // LU decomposition for a banded system
  // See pp. 45--48 in Hoffman, Numerical Methods for Engineers and Scientists

  int max_i = 0;
  int max_k = 0;

  for (int j = 0; j<_NFD-1; j++) // columns
  { 

    max_i = std::min(j+_BW+1, _NFD);

    for (int i = j+1; i<max_i; i++) // rows
    {

      // M = T[i, j] . T[j, j]^{-1}
      tmpmat = _Data[_BW][j];
      tmpmat2.invert(tmpmat);
      M.dot(_Data[j-i+_BW][j], tmpmat2);

      // T[i, j] = M
      _Data[j-i+_BW][j] = M;

      max_k = std::min( j+_BW+1, _NFD );

      for (int k = j+1; k<max_k; k++) // columns
      {

        // T[i, k] = T[i, k] - M.T[j, k]
        tmpmat.dot(M, _Data[k-j+_BW][k]);
        _Data[k-i+_BW][k] -= tmpmat;

      } 

    }
  }

  _LU = true;

} // }}}

// when solving for a RHS vector (non-periodic)
// T = this
void MatInterFieldMat::band_solve( VecInterFieldVec & b, VecInterFieldVec & x )
{ // {{{

  InterFieldVec tmpvec(_NIcomp);
  InterFieldVec tmpsum(_NIcomp);

  tmpvec.setflag_inrealspace(false);
  tmpsum.setflag_inrealspace(false);

  if (_LU == false)
  {
    std::cout << "*** WARNING: Cannot solve Tx = b"; 
    std::cout << " unless LU decomposition is performed first.***" << std::endl;
  }
  else // _LU == true
  {

    int min_j = 0;
    int max_j = 0;

    // --- forward substitution (Ly = b) ---
    //  * do this in place so that b <- y
    // i==0 is unchanged

    for (int i=1; i<_NFD; i++)
    {

      tmpsum.zero();
      min_j = std::max(0, i-_BW);

      for (int j=min_j; j<i; j++)
      {
        tmpvec.dot( _Data[j-i+_BW][j], b[j] );
        tmpsum += tmpvec;
      }
      b[i] -= tmpsum;
    }

    // --- back substitution (Ux = y) ---

    for (int i=_NFD-1; i>=0; i--)
    {

      tmpsum = b[i];
      max_j = std::min(i+_BW+1, _NFD);

      for (int j=i+1; j<max_j; j++)
      {
        // tmpsum = tmpsum - T[i, j] . x[j]
        tmpvec.dot( _Data[j-i+_BW][j], x[j] );
        tmpsum -= tmpvec;
      }

      // x[m] = T[i, i]^{-1} . tmpsum
      x[i].setflag_inrealspace(false);
      x[i].linsolve(_Data[_BW][i], tmpsum);

    }

  }

} // }}}

// when solving for a RHS matrix (periodic)
// T = this
void MatInterFieldMat::band_solve( MatInterFieldMat & B, MatInterFieldMat & X )
{ // {{{

  InterFieldMat tmpmat(_NIcomp, _NIcomp);
  InterFieldMat tmpsum(_NIcomp, _NIcomp);

  tmpmat.setflag_inrealspace(false);
  tmpsum.setflag_inrealspace(false);

  // *** Recall:  T[i'][j'] = T'[j',j'+i'-2] ***

  if (_LU == false)
  {
    std::cout << "*** WARNING: Cannot solve Tx = b unless LU decomposition is performed first.***" << std::endl;
  }
  else // _LU == true
  {

    // --- loop over all columns of B ---
    for( int n = 0; n < B.Ncol(); n++ )
    {

      int min_j = 0;
      int max_j = 0;

      // --- forward substitution (LY = B) ---
      //  * do this in place so that B <- Y
      // i==0 is unchanged

      for (int i=1; i<_NFD; i++)
      {

        tmpsum.zero();
        min_j = std::max(0, i-_BW);

        for (int j=min_j; j<i; j++)
        {
          tmpmat.dot( _Data[j-i+_BW][j], B[j][n] );
          tmpsum += tmpmat;
        }
        B[i][n] -= tmpsum;
      }

      // --- back substitution (UX = Y) ---

      for (int i=_NFD-1; i>=0; i--)
      {

        tmpsum = B[i][n];
        max_j = std::min(i+_BW+1, _NFD);

        for (int j=i+1; j<max_j; j++)
        {
          // tmpsum = tmpsum - T[i, j] . x[j]
          tmpmat.dot( _Data[j-i+_BW][j], X[j][n] );
          tmpsum -= tmpmat;
        }

        // X[i][n] = T[i, i]^{-1} . tmpsum
        X[i][n].setflag_inrealspace(false);
        tmpmat.invert( _Data[_BW][i] );
        X[i][n].dot(tmpmat, tmpsum);
      }

    }

  }

} // }}}

void MatInterFieldMat::matmul( MatInterFieldMat & A, MatInterFieldMat & B )
{ // {{{

  InterFieldMat tmp(A.NFMrow(), B.NFMcol());

  if ( A.getflag_inrealspace() != B.getflag_inrealspace() or
       A.getflag_inrealspace() != this->getflag_inrealspace() )
  {
    std::cout << "***(Error in MatInterFieldMat::matmul)***\n";
    std::cout << "Incompatible Fourier Representation\n";
    exit(1);
  }

  if ( A.Ncol() != B.Nrow() or A.Nrow() != _Nrow or B.Ncol() != _Ncol )
  {
    std::cout << "***(Error in MatInterFieldMat::matmul)***\n";
    std::cout << "Incompatible Matrix Sizes\n";
    exit(1);
  }

  bool rspaceflag = A.getflag_inrealspace();

  for (UInt i=0; i<A.Nrow(); ++i) // over rows in A
  {
    for (UInt j=0; j<B.Ncol(); ++j) // columns in B
    {
      tmp.zero();
      tmp.setflag_inrealspace(rspaceflag);
      for (UInt k=0; k<B.Nrow(); ++k) // columns in A (or rows in B)
      {
        tmp.dot(A[i][k], B[k][j]);
        _Data[i][j] += tmp;
      }
    }
  }

} // }}}

void MatInterFieldMat::transpose()
{ // {{{

  InterFieldMat tmp(_NFMrow, _NFMcol);

  for( UInt i=0; i<_Nrow; i++)
  {
    for( UInt j=i+1; j<_Ncol; j++)
    {
      tmp = _Data[i][j];
      _Data[i][j] = _Data[j][i];
      _Data[j][i] = tmp;
    }
  }

} // }}}

FieldType MatInterFieldMat::max()
{ // {{{

  // returns element with maximum abs value
  RealType fieldmax = 0.;
  RealType localmax = 0.;
  UInt max_i = 0;
  UInt max_j = 0;

  for( UInt i=0; i<_Nrow; i++)
  {
    for( UInt j=0; j<_Ncol; j++)
    {
      localmax = std::abs( _Data[i][j].max() );
      if (localmax > fieldmax)
      {
        fieldmax = localmax;
        max_i = i;
        max_j = j;
      }
    }
  }
  return _Data[max_i][max_j].max();

} // }}}

void MatInterFieldMat::ReportState()
{ // {{{

  FILE* file1;
  file1 = fopen("T.dat", "w");

  std::cout << "Reporting: " << "T.dat" << "\n";
  //std::cout << "Nrow: " << _Nrow << "\n";
  //std::cout << "Ncol: " << _Ncol << "\n";
  //std::cout << "NFMrow: " << _NFMrow << "\n";
  //std::cout << "NFMcol: " << _NFMcol << "\n";
  //std::cout << "LU: " << _LU << "\n\n";

  fprintf(file1, "# Nrow: %d\n",  _Nrow);
  fprintf(file1, "# Ncol: %d\n",  _Ncol);
  fprintf(file1, "# NFMrow: %d\n", _NFMrow);
  fprintf(file1, "# NFMcol: %d\n", _NFMcol);
  fprintf(file1,"# LU: %d\n", _LU);

  for (int i = 0; i<_Nrow; i++)
  {
    //std::cout << "Band: " << i << "\n";
    for (int j = 0; j<_Ncol; j++)
    {
      //std::cout << "x_j = " << j << ": ";
      for (int k = 0; k<_NFMrow; k++)
      {
        //if (k > 0) std::cout << "        ";
        for (int l = 0; l<_NFMrow; l++)
        {
          //std::cout << "(" << i << ", " << j << ", " << k << ", " << l << "): "; 
          //std::cout << _Data[i][j][k][l].getelement(0) <<"\n";

          fprintf(file1, "%25.16e", _Data[i][j][k][l].getelement(0).real());
          //std::cout << _Data[i][j][k][l].getelement(0).real();
          //if (l < _NFMrow) std::cout << ", ";
        }
        //std::cout << "\n";
      }
    }
    fprintf(file1, "\n");
  }

  fclose(file1);

} // }}}

