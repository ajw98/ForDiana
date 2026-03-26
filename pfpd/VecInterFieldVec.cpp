/***************************************************************
*
* Class for Interleaved Banded Field Matrices
* 
* DRT -- Fri, 20 Mar 2017
*
****************************************************************/

#include "VecInterFieldVec.h"

// ----------------- Constructor/Destructor -----------------

VecInterFieldVec::VecInterFieldVec( int nfd, int ncomp ) :
_NFD(nfd),
_NIcomp(ncomp),
_Data( nfd, InterFieldVec( ncomp) )
{ // {{{
} // }}}

VecInterFieldVec::VecInterFieldVec( const VecInterFieldVec &rhs ):
_NFD( rhs._NFD ),
_NIcomp( rhs._NIcomp ),
_Data( rhs._Data )
{ // {{{
} // }}}

VecInterFieldVec::~VecInterFieldVec()
{ // {{{
} // }}}

// ----------------- Operators  -----------------

VecInterFieldVec & VecInterFieldVec::operator=( const VecInterFieldVec & rhs )
{ // {{{

  this->_NFD = rhs.NFD();
  this->_NIcomp = rhs.NIcomp();
  this->_Data = rhs._Data;

  return *this;

} // }}}

// ----------------- Member Functions -----------------

void VecInterFieldVec::zero()
{ // {{{

  for (int j = 0; j < _NFD; j++)
  {
    _Data[j] = FieldType(0.);
  }

} // }}}

void VecInterFieldVec::setflag_inrealspace( bool realspaceflag )
{ // {{{

  for (int j = 0; j < _NFD; j++)
  {
    _Data[j].setflag_inrealspace(realspaceflag);
  }

} // }}}

bool VecInterFieldVec::getflag_inrealspace() const
{ // {{{

  bool rspace_flag = true;
  bool tmp_flag = true;

  for (int j = 0; j < _NFD; j++)
  {
    tmp_flag = _Data[j].getflag_inrealspace();
    if (j > 0 and rspace_flag != tmp_flag)
    {
      std::cout << "Warning: SmartFieldVec is mixed real space/fourier space. Use component functions.";
    }
    rspace_flag = tmp_flag;
  }

  return rspace_flag;

} // }}}

void VecInterFieldVec::smart_to_inter( SmartFieldVec & v )
{ // {{{

  if ( _NFD != v[0].GetSize() and _NIcomp != v.get_Nelem() )
  {
    std::cout << "***(Error in VecInterFieldVec::smart_to_inter)***\n";
    std::cout << "Incompatible Matrix Sizes\n";
    exit(1);
  }

  this->zero();
  this->setflag_inrealspace( v.getflag_inrealspace() );

  for ( int i = 0; i < _NIcomp; i++ )
  {
    for ( int j = 0; j < _NFD; j++ )
    {
      _Data[j][i] = v[i][j];
    }
  }

} // }}}

void VecInterFieldVec::inter_to_smart( SmartFieldVec & v )
{ // {{{

  v.zero();
  v.setflag_inrealspace( this->getflag_inrealspace() );

  for ( int i = 0; i < _NIcomp; i++ )
  {
    for ( int j = 0; j < _NFD; j++ )
    {
      v[i][j] = _Data[j][i];
    }
  }

} // }}}

void VecInterFieldVec::ReportState(std::string filename )
{ // {{{

  FILE* file1;
  file1 = fopen(filename.c_str(), "w");

  fprintf(file1, "# NIcomp: %d\n",  _NIcomp);
  fprintf(file1, "# NFD: %d\n", _NFD);

  for (int i = 0; i<_NFD; i++)
  {
    for (int j = 0; j<_NIcomp; j++)
    {
      fprintf(file1, "%25.16e", _Data[i][j].getelement(0).real());
    }
  }

  fclose(file1);

} // }}}

