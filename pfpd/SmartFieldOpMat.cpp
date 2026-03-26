/***************************************************************
*
* Class for the time integration scheme
* 
* DRT -- Fri, 13 Mar 2015
*
****************************************************************/

#include "SmartFieldOpMat.h"

// ----------------- Constructor/Destructor -----------------

SmartFieldOpMat::SmartFieldOpMat( UInt nrow, UInt ncol ) :
_Nrow( nrow ),
_Ncol( ncol ),
_Data( nrow, Vec1dSFO( ncol, SmartFieldOp() ) )
{ // {{{
} // }}}

SmartFieldOpMat::SmartFieldOpMat( const SmartFieldOpMat &rhs ):
_Ncol( rhs._Ncol ),
_Nrow( rhs._Nrow ),
_Data( rhs._Data )
{ // {{{
} // }}}

SmartFieldOpMat::~SmartFieldOpMat()
{ // {{{
} // }}}

// ----------------- Operators  -----------------

SmartFieldOpMat & SmartFieldOpMat::operator=( const SmartFieldOpMat & rhs )
{ // {{{

  this->_Ncol = rhs._Ncol;
  this->_Nrow = rhs._Nrow;
  this->_Data = rhs._Data;

  return *this;

} // }}}

SmartFieldOpMat & SmartFieldOpMat::operator+=( const SmartFieldOpMat & rhs)
{ // {{{

  if (this->_Nrow != rhs._Nrow )
  {
    std::cout << "*** Error in SmartFieldOpMat::operator += ***\n";
    std::cout << "Sizes are not compatible\n";
    exit(1);
  }

  if (this->getflag_inrealspace() != rhs.getflag_inrealspace() )
  {
    std::cout << "*** Error in SmartFieldOpMat::operator += ***\n";
    std::cout << "Fourier representations are not compatible\n";
    exit(1);
  }

  for (UInt i=0; i<_Nrow; ++i)
  {
    for (UInt j=0; j<_Ncol; ++j)
    {
      this->_Data[i][j] += rhs._Data[i][j];
    }
  }

  return *this;

} // }}}

SmartFieldOpMat & SmartFieldOpMat::operator*=( const FieldType scalar )
{ // {{{

  for (UInt i=0; i<_Nrow; ++i)
  {
    for (UInt j=0; j<_Ncol; ++j)
    {
      this->_Data[i][j] *= scalar;
    }
  }

  return *this;

} // }}}

// ----------------- Member Functions -----------------

void SmartFieldOpMat::AddBand( int bandnum, const SmartFieldMat & Band )
{ // {{{

  if ( _Nrow != Band.get_Nrow() )
  {
    std::cout << "*** Error in SmartFieldOpMat::AddBand()***\n";
    std::cout << "Operator and SmartFieldMat sizes are not compatible\n";
    exit(1);
  }

  for (int i=0; i<_Nrow; i++)
  {
    for (int j=0; j<_Ncol; j++)
    {
      _Data[i][j].AddBand(bandnum, Band(i,j));
    }
  }

} // }}}

void SmartFieldOpMat::AddBand( int bandnum, const Vec2dFieldType & Band )
{ // {{{

  if ( _Nrow != Band.size() )
  {
    std::cout << "*** Error in SmartFieldOpMat::AddBand()***\n";
    std::cout << "Operator and Vec2dFieldType sizes are not compatible\n";
    exit(1);
  }

  for (int i=0; i<_Nrow; i++)
  {
    for (int j=0; j<_Ncol; j++)
    {
      _Data[i][j].AddBand(bandnum, Band[i][j]);
    }
  }

} // }}}

void SmartFieldOpMat::GetBand( int bandnum, SmartFieldMat & Band ) const
{ // {{{

  if ( _Nrow != Band.get_Nrow() )
  {
    std::cout << "*** Error in SmartFieldOpMat::GetBand()***\n";
    std::cout << "Operator and SmartFieldMat sizes are not compatible\n";
    exit(1);
  }

  for (int i=0; i<_Nrow; i++)
  {
    for (int j=0; j<_Ncol; j++)
    {
      _Data[i][j].GetBand(bandnum, Band(i,j));
    }
  }

} // }}}

void SmartFieldOpMat::zero()
{ // {{{

  for (int i=0; i<_Nrow; i++)
  {
    for (int j=0; j<_Ncol; j++)
    {
      _Data[i][j].zero();
    }
  }

} // }}}

void SmartFieldOpMat::setflag_inrealspace( bool realspaceflag )
{ // {{{

  for (int i=0; i<_Nrow; i++)
  {
    for (int j=0; j<_Ncol; j++)
    {
      _Data[i][j].setflag_inrealspace(realspaceflag);
    }
  }

} // }}}

bool SmartFieldOpMat::getflag_inrealspace() const
{ // {{{

  bool rspace_flag = true;
  bool tmp_flag = true;

  for(UInt i=0; i<_Nrow; ++i)
  {
    for(UInt j=0; j<_Ncol; ++j)
    {
      tmp_flag = this->_Data[i][j].getflag_inrealspace();
      if (i > 0 and j > 0 and rspace_flag != tmp_flag)
      {
        std::cout << "Warning: SmartFieldVec is mixed real space/fourier space. Use component functions.";
      }
      rspace_flag = tmp_flag;
    }
  }

  return rspace_flag;

} // }}}

void SmartFieldOpMat::OpDot( const SmartFieldVec & in, SmartFieldVec & out )
{ // {{{

  if (this->_Ncol != in.get_Nelem() or this->_Ncol != out.get_Nelem())
  {
    std::cout << "***(Error) in SmartFieldOpMat::OpDot***\n";
    std::cout << "Incompatible operator and fieldvec size.\n";
    exit(1);
  }

  bool rspaceflag = this->getflag_inrealspace();

  if ( rspaceflag != in[0].getflag_inrealspace())
  {
    std::cout << "***(Error) in SmartFieldOpMat::OpDot***\n";
    std::cout << "Both operator and field must be in the same Fourier representation.\n";
    exit(1);
  }

  SmartField tmp;
  tmp.setflag_inrealspace( rspaceflag );

  out.setflag_inrealspace( rspaceflag );
  out = FieldType(0.);

  for (int i=0; i<this->_Nrow; i++)
  {
    for (int j=0; j<this->_Ncol; j++)
    {
      tmp = FieldType(0.);
      (*this)[i][j].OpDot( in[j], tmp );
      out[i] += tmp;
    }
  }

} // }}}

void SmartFieldOpMat::ReportState()
{ // {{{

  for (int i=0; i<_Nrow; i++)
  {
    for (int j=0; j<_Ncol; j++)
    {
      std::cout << "row: " << i << ", col: " << j << "\n";
      _Data[i][j].ReportState();
    }
  }
  
} // }}}

