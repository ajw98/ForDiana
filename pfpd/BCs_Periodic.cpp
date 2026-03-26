/***************************************************************
*
* Periodic boundary condition solver
* 
* DRT -- Sat 15 Apr 2017
*
****************************************************************/

#include "BCs_Periodic.h"

// ----------------- Constructor/Destructor -----------------

BCs_Periodic :: BCs_Periodic( std::string filename, int nicomp,
                              Grid * gridarg ) : 
BCs_Base( nicomp, gridarg )
{ // {{{

  // Read in Boundary conditions
  std::string tmp_str;
  
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
      std::cout << "cannot open " << filename << std::endl;
      exit(1);
    }
  }
  else
  {//Legacy file format
    std::ifstream f2(filename.c_str());
    if ( f2.is_open() )
    {

      std::getline(f2, tmp_str, '#');
      _BCFlag = atoi(tmp_str.c_str());
      std::getline(f2, tmp_str, '\n');

    }
    else
    {
      std::cout << "cannot open " << filename << std::endl;
      exit(1);
    }
  }

  std::cout << "\n" << std::scientific;
  std::cout << " * Boundary Conditions Initiated\n" << std::scientific;
  std::cout << "   - BC Parameters:\n" << std::scientific;
  std::cout << "      BCFlag             = " << _BCFlag << std::endl;

} // }}}

BCs_Periodic :: ~BCs_Periodic()
{ // {{{
} // }}}

void BCs_Periodic::solve_BCs( SmartFieldOpMat & T, 
                              SmartFieldVec & rhs, 
                              SmartFieldVec & x,
                              RealType t)
{ // {{{

  // *********************
  //    solve T.x = rhs
  // *********************
  // where T is a cyclic, block-banded matrix of bandwidth, BW 
  // 
  // T_int: incoming intercalated block-banded matrix 
  //        (Nbands, Nx, NIcomp, NIcomp, Ny/Nz)
  //
  // [ B C        ...   A ]
  // [ A B C      ...     ]
  // [   A B C    ...     ]
  // [         ...        ]
  // [     ...    A B C   ]
  // [     ...      A B C ]
  // [ C   ...        A B ]
  //
  // where:
  //
  // Nbands = 2*BW + 1
  //
  // A (corner) = upper triangular (BW x BW) block
  //  ** Note: each diagonal D_i is a NIcomp x NIcomp block
  // 
  // [   Dn  Dn-1  Dn-2  ...       D1 ]
  // [    0    Dn  Dn-1  ...       D2 ]
  // [    0     0    Dn  ...       D3 ]
  // [                   ...          ]
  // [    0     0     0  ...       Dn ]  
  // 
  // B (diagonal) = square (BW x BW) block
  // 
  // [   D0   D1   D2   ...      Dn-2  Dn-1 ]
  // [   D1   D0   D1   D2   ...       Dn-2 ]
  // [   D2   D1   D0   D1   D2   ...       ]
  // [                  ...                 ]
  // [ Dn-1 Dn-2        ... Dn-2 Dn-1    D0 ]
  // 
  // C (corner) = lower triangular (BW x BW) block
  // 
  // [   Dn     0    0   ...         0 ]
  // [ Dn-1    Dn    0   ...         0 ]
  // [ Dn-2  Dn-1   Dn   ...         0 ]
  // [                   ...           ]
  // [   D1    D2   D3   ...        Dn ]
  //
  // For SMW, we will construct the following
  //
  // T = T' + UV^{T}
  //
  // U = [ s I ]
  //     [  0  ]
  //     [ ... ]
  //     [  0  ]
  //     [  C  ]
  //
  // V^{T} = [ I 0 ... 0 A/s ]
  //
  // T' = T - UV^{T}
  // 
  // [ B-sI  C          ...           0 ]
  // [    A  B  C       ...           0 ]
  // [       A  B  C    ...           0 ]
  // [             ...                  ]
  // [    0     ...    C  B  A          ]
  // [    0     ...       C  B  A       ]
  // [    0     ...          C  B-C.A/s ]
  //
  // and then solve:
  // T . y = rhs
  // T . U = Z
  // x = y - Z.(I + V^{T}.Z)^{-1}.V^{T}.y

  UInt NIcomp = T.Nrow();
  UInt BW = T.BW();
  UInt FDSize = T.FDSize();

  MatInterFieldMat T_int( 2*BW+1, FDSize, NIcomp, NIcomp );
  VecInterFieldVec rhs_int( FDSize, NIcomp );
  VecInterFieldVec x_int( FDSize, NIcomp );

  T_int.smart_to_inter( T );
  rhs_int.smart_to_inter( rhs );

  // In this case the system is diagonal. Much easier!
  if (BW == 0)
  {
  
    T_int.band_LU();
    T_int.band_solve( rhs_int, x_int ); // T.x = rhs
    x_int.inter_to_smart( x );

  }
  else
  {

    // --- get A, C, C.A ---
    MatInterFieldMat A( BW, BW, NIcomp, NIcomp);
    MatInterFieldMat C( BW, BW, NIcomp, NIcomp);
    MatInterFieldMat CdotA( BW, BW, NIcomp, NIcomp);

    A.setflag_inrealspace(false);
    C.setflag_inrealspace(false);
    CdotA.setflag_inrealspace(false);

    // Get A (below diagonal)
    int bandnum = 0;
    for (int i=0; i<BW; i++)
    {
      for(int j=0; j<BW; j++)
      {

        bandnum = j-i;
        if (bandnum >= 0 )
        {
          A[i][j] = T_int[bandnum][j];
        }
        else
        {
          A[i][j].zero();
        }

      }
    }

    // find max(A), and scale A by max value
    RealType s = std::abs( A.max() );
    for (int i=0; i<BW; i++)
    {
      for(int j=0; j<BW; j++)
      {
        A[i][j] *= 1./s;
      }
    }

    // Get C (above diagonal)
    for (int i=0; i<BW; i++)
    {
      for(int j=0; j<BW; j++)
      {

        bandnum = j-i+BW + BW;
        if (bandnum < (2*BW + 1) )
        {
          C[i][j] = T_int[bandnum][j];
        }
        else
        {
          C[i][j].zero();
        }

      }
    }

    CdotA.matmul( C, A );

    // --- modify T->A, get U, get V ---
    
    // Modify T

    // top half: -s I
    InterFieldMat tmp( NIcomp, NIcomp );
    tmp.setflag_inrealspace( false );
    tmp.eye( s );

    for(int j=0; j<BW; j++)
    {
      T_int[BW][j] -= tmp;
    }

    //  bottom half: - CT.C/s
    for(int i=0; i<BW; i++)
    {
      for(int j=0; j<BW; j++)
      {
        T_int[j-i+BW][FDSize-BW+j] -= CdotA[i][j];
      }
    }

    // set U
    MatInterFieldMat U( FDSize, BW, NIcomp, NIcomp);
    U.setflag_inrealspace(false);
    U.zero();

    for(int j=0; j<BW; j++)
    {
      U[j][j].eye(s);
    }

    for(int i=0; i<BW; i++)
    {
      for(int j=0; j<BW; j++)
      {
        U[FDSize-BW+i][j] = C[i][j];
      }
    }

    // set V
    MatInterFieldMat VT( BW, FDSize, NIcomp, NIcomp);
    VT.setflag_inrealspace(false);
    VT.zero();

    for(int j=0; j<BW; j++)
    {
      VT[j][j].eye(FieldType(1.));
    }

    for(int i=0; i<BW; i++)
    {
      for(int j=0; j<BW; j++)
      {
        VT[i][FDSize-BW+j] = A[i][j];
      }
    }

    // --- solve the linear sytem using a banded solver ---

    MatInterFieldMat Z( FDSize, BW, NIcomp, NIcomp );
    VecInterFieldVec y( FDSize, NIcomp );

    Z.setflag_inrealspace(false); 
    y.setflag_inrealspace(false); 

    T_int.band_LU();
    T_int.band_solve( rhs_int, y ); // T.y = rhs
    T_int.band_solve( U, Z ); // T.U = Z

    // --- find x using the Sherman-Morrison-Woodbury Formula ---
    // x = y - Z.(I + V^{T}.Z)^{-1}.V^{T}.y
    solve_SMW( x_int, y, Z, VT );
    x_int.inter_to_smart( x );
  }

} // }}}

void BCs_Periodic::solve_SMW( VecInterFieldVec & x, 
                              VecInterFieldVec & y, 
                              MatInterFieldMat & Z, 
                              MatInterFieldMat & VT )
{ // {{{

  // In:    y, Z, V^T
  // Todo:  Use Sherman-Morrison-Woodbury formula to calculate
  //            x = y - Z.H^{-1}.V^{T}.y
  //        where,
  //            V^T : (BW x NFD) 
  //            V^T = [ I 0 ... 0 C/s ]
  //        and,
  //            H : (BW x BW) 
  //            H = I + V^T.Z
  //
  // Out:   x

  UInt NFD = Z.Nrow();
  UInt BW = Z.Ncol();
  UInt NIcomp = Z.NFMrow();

  // BW x BW matrix
  Vec2dIFM H( BW, Vec1dIFM( BW, InterFieldMat( NIcomp, NIcomp )));
  Vec2dIFM Hinv( BW, Vec1dIFM( BW, InterFieldMat( NIcomp, NIcomp )));
  Vec2dIFM BigEye( BW, Vec1dIFM( BW, InterFieldMat( NIcomp, NIcomp )));

  // --- find H ---
  // H = I + V^T.Z (BW x BW matrix)
  {
    InterFieldMat tmpmat( NIcomp, NIcomp);
    tmpmat.setflag_inrealspace(false);

    for ( UInt i = 0; i < BW; i++ )
    {
      for ( UInt j = 0; j < BW; j++ )
      {
        H[i][j].setflag_inrealspace(false);
        if (i == j)
        {
          H[i][j].eye( FieldType(1.) );
        }
        else
        {
          H[i][j].zero();
        }
        for ( UInt m = 0; m < NFD; m++ )
        {
          tmpmat.dot(VT[i][m], Z[m][j]);
          H[i][j] += tmpmat;
        }
      }
    }
  }

  // --- Invert H ---
  // Solve H Hinv = Eye

  // Initialize BigEye
  for ( UInt i = 0; i< BW; ++i )
  {
    for ( UInt j = 0; j< BW; ++j )
    {
      BigEye[i][j].setflag_inrealspace(false);
      if (i == j)
      {
        BigEye[i][j].eye(FieldType(1.));
      }
      else
      {
        BigEye[i][j].zero();
      }
    }
  }

  // Do LU decomposition
  // L (U Hinv) = Eye
  // let: U Hinv = Y
  {
    InterFieldMat tmpmat( NIcomp, NIcomp);
    InterFieldMat tmpmat2( NIcomp, NIcomp);
    tmpmat.setflag_inrealspace(false);
    tmpmat2.setflag_inrealspace(false);

    InterFieldMat M( NIcomp, NIcomp ); // multipliers
    M.setflag_inrealspace(false);

    for( UInt k = 0; k < BW-1; k++ )
    {
      for (UInt i = k+1; i < BW; i++ )
      {

        tmpmat2 = H[k][k]; 
        tmpmat.invert( tmpmat2 );
        M.dot(H[i][k], tmpmat);
        H[i][k] = M;

        for (UInt j = k+1; j < BW; j++ )
        {
          tmpmat.dot(M, H[k][j]);
          H[i][j] -= tmpmat;
        }
      }
    }
  }

  // Forward substitution
  // L Y = Eye
  // (Do in place Eye <- Y)
  // Y_ij = Eye_ij - Sum_{k=0}^{i-1} L_ik . Eye_kj
  {
    InterFieldMat tmpmat( NIcomp, NIcomp);
    InterFieldMat tmpsum( NIcomp, NIcomp);
    tmpmat.setflag_inrealspace(false);
    tmpsum.setflag_inrealspace(false);

    for( UInt j = 0; j < BW; j++ ) // loop over columns of Eye
    {
      // loop over rows of L
      // start at i = 1, since the i=0 row is unchanged
      for( UInt i = 1; i < BW; i++ ) 
      {
        tmpsum = BigEye[i][j];
        // loop over columns L (from 0 to before diagonal)
        for( UInt k = 0; k < i; k++ )
        {
          tmpmat.dot( H[i][k], BigEye[k][j] );
          tmpsum -= tmpmat;
        }
        BigEye[i][j] = tmpsum;
      }
    }
  }

  // Back substitution
  // U Hinv = Y
  // Hinv_ij = Y_ij - Sum_{k=i+1}^{N-1} U_ik . Y_kj
  {
    InterFieldMat tmpmat( NIcomp, NIcomp);
    InterFieldMat tmpsum( NIcomp, NIcomp);
    tmpmat.setflag_inrealspace(false);
    tmpsum.setflag_inrealspace(false);

    for( UInt j = 0; j < BW; j++ ) // loop over columns of Eye
    {
      // loop (backwords) over rows of L
      // get the bottom corner, because we don't need a sum.
      Hinv[BW-1][j].setflag_inrealspace(false);
      tmpmat.invert(H[BW-1][BW-1]);
      Hinv[BW-1][j].dot(tmpmat, BigEye[BW-1][j]);

      // loop (up) over rest of the rows
      // (need int or enters infinte loop)
      for( int i = BW-2; i >= 0; i-- ) 
      {
        tmpsum = BigEye[i][j];
        // loop over columns U (from diag+1 to N-1)
        for( UInt k = i+1; k < BW; k++ )
        {
          tmpmat.dot(H[i][k],Hinv[k][j]);
          tmpsum -= tmpmat;
        }
        Hinv[i][j].setflag_inrealspace(false);
        tmpmat.invert( H[i][i] );
        Hinv[i][j].dot(tmpmat, tmpsum);
      }
    }
  }

  // --- find Z.H^{-1} ---
  //  (NFD x BW).(BW x Bw) = (NFD x Bw)
  Vec2dIFM ZHinv(NFD, Vec1dIFM(BW, InterFieldMat( NIcomp, NIcomp)) );
  {

    InterFieldMat tmpmat( NIcomp, NIcomp);
    InterFieldMat tmpsum( NIcomp, NIcomp);
    tmpmat.setflag_inrealspace(false);
    tmpsum.setflag_inrealspace(false);

    for ( UInt i = 0; i < NFD; i++ )
    {
      for ( UInt j = 0; j < BW; j++ )
      {
        tmpsum = FieldType(0.);
        for ( UInt m = 0; m < BW; m++ )
        {
          tmpmat.dot(Z[i][m], Hinv[m][j]);
          tmpsum += tmpmat;
        }
        ZHinv[i][j].setflag_inrealspace(false);
        ZHinv[i][j] = tmpsum;
      }
    }
  }
  
  // --- find V^{T}.y ---
  //  (BW x NFD ).(NFD x 1) = (_Bw x 1)
  Vec1dIFV VTy ( BW, InterFieldVec(NIcomp)); // BW x 1
  {
    InterFieldVec tmpvec(NIcomp);
    InterFieldVec vecsum(NIcomp);

    tmpvec.setflag_inrealspace(false);
    vecsum.setflag_inrealspace(false);

    for ( UInt i = 0; i < BW; i++ )
    {
      vecsum = FieldType(0.);
      for ( UInt m = 0; m < NFD; m++ )
      {
        tmpvec.dot(VT[i][m], y[m]);
        vecsum += tmpvec;
      }
      VTy[i].setflag_inrealspace(false);
      VTy[i] = vecsum;
    }

    // --- find x ---
    // x = y - Z.H^{-1}.V^{T}.y
    // (NFD x 1) = (NFD x 1) - (NFD x BW) . (BW x 1)
    for ( UInt i = 0; i < NFD; i++ )
    {
      x[i] = y[i];

      vecsum = FieldType(0.);
      for ( UInt m = 0; m < BW; m++ )
      {
        tmpvec.dot(ZHinv[i][m], VTy[m]);
        vecsum += tmpvec;
      }

      x[i] -= vecsum;
    }

  }

} // }}}

