/***************************************************************
*
* Class for the time integration scheme
* 
* DRT -- Fri, 13 Mar 2015
*
****************************************************************/

#include "SmartFieldOp.h"

// ----------------- Constructor/Destructor -----------------

SmartFieldOp::SmartFieldOp() :
_BW(0),
_BandStencil( 1, 0 ),
_Data( 1, SmartField() )
{ // {{{
} // }}}

SmartFieldOp::SmartFieldOp( const SmartFieldOp &rhs ):
_BW( rhs._BW ),
_BandStencil( rhs._BandStencil ),
_Data( rhs._Data )
{ // {{{
} // }}}

SmartFieldOp::~SmartFieldOp()
{ // {{{
} // }}}

// ----------------- Operators  -----------------

SmartFieldOp & SmartFieldOp::operator=( const SmartFieldOp & rhs )
{ // {{{

  this->_BW = rhs.BW();
  this->_BandStencil = rhs._BandStencil;
  this->_Data = rhs._Data;

  return *this;

} // }}}

SmartFieldOp & SmartFieldOp::operator+=( const SmartFieldOp & rhs)
{ // {{{

  if (this->_Data[0].getflag_inrealspace() != rhs._Data[0].getflag_inrealspace() )
  {
    std::cout << "*** Error in SmartFieldOp::operator += ***\n";
    std::cout << "Fourier representations are not compatible\n";
    exit(1);
  }

  for (UInt j=0; j<rhs._Data.size(); ++j)
  {
    this->AddBand( rhs._BandStencil[j], rhs._Data[j]);
  }

  return *this;

} // }}}

SmartFieldOp & SmartFieldOp::operator*=( const FieldType scalar )
{ // {{{

  for (size_t j=0; j<_Data.size(); j++)
  {
    _Data[j] *= scalar;
  }

  return *this;

} // }}}

// ----------------- Member Functions -----------------

int SmartFieldOp::GetStencilPos( int bandnum ) const
{ // {{{

  for (int i=0; i<_BandStencil.size(); i++)
  {
    if (bandnum == _BandStencil[i])
    {
      return i;
    }
  }
  return -1;

} // }}}

void SmartFieldOp::SetBand( int bandnum, SmartField & Band )
{ // {{{

  int stencilpos = GetStencilPos(bandnum);

  // if band is new, make space and add it
  if ( stencilpos < 0 )
  {
    _BandStencil.push_back( bandnum );
    _Data.push_back( Band );
  }
  else // if not, add it to existing band
  {
    _Data[ stencilpos ] = Band;
  }

  // get new bandwidth
  for (int i=0; i<_BandStencil.size(); i++)
  {
    _BW = std::max(_BW, std::abs(_BandStencil[i]));
  }

} // }}}

void SmartFieldOp::AddBand( int bandnum, const SmartField & Band )
{ // {{{

  int stencilpos = GetStencilPos(bandnum);
  
  // if band is new, make space and add it
  if ( stencilpos < 0 )
  {
    _BandStencil.push_back( bandnum );
    _Data.push_back( Band );
  }
  else // if not, add it to existing band
  {
    _Data[ stencilpos ] += Band;
  }

  // get new bandwidth
  for (int i=0; i<_BandStencil.size(); i++)
  {
    _BW = std::max(_BW, std::abs(_BandStencil[i]));
  }

} // }}}

void SmartFieldOp::AddBand( int bandnum, FieldType Band )
{ // {{{

  int stencilpos = GetStencilPos(bandnum);
  
  // if band is new, make space and add it
  if ( stencilpos < 0 )
  {
    _BandStencil.push_back( bandnum );
    _Data.push_back( Band );
  }
  else // if not, add it to existing band
  {
    _Data[ stencilpos ] += Band;
  }

  // get new bandwidth
  for (int i=0; i<_BandStencil.size(); i++)
  {
    _BW = std::max(_BW, std::abs(_BandStencil[i]));
  }

} // }}}

void SmartFieldOp::GetBand( int bandnum, SmartField & Band ) const
{ // {{{

  Band = _Data[ GetStencilPos(bandnum) ];

} // }}}

void SmartFieldOp::OpDot( const SmartField & in, SmartField & out )
{ // {{{

  if (this->_Data[0].getflag_inrealspace() != in[0].getflag_inrealspace())
  {
    std::cout << "***(Error) in SmartFieldOp::OpDot***\n";
    std::cout << "Both operator and field must be in the same Fourier representation.\n";
    exit(1);
  }

  out.setflag_inrealspace( _Data[0].getflag_inrealspace() );
  out = FieldType(0.);

  SmartField tmp;
  tmp.setflag_inrealspace( _Data[0].getflag_inrealspace() );
  std::vector<UInt> stencil(in.GetSize(), 0);
  int idx = 0;

  for (int i=0; i<this->_BandStencil.size(); i++)
  {
    // get a stencil for evaluating "in"
    for (int j=0; j<in.GetSize(); j++)
    {
      idx = (_BandStencil[i] + j + in.GetSize())%in.GetSize();
      stencil[j] = UInt(idx);
    }

    tmp = FieldType(0.);
    tmp.subscript( in, stencil );
    tmp *= _Data[i];
    out += tmp;
  }

} // }}}

void SmartFieldOp::zero()
{ // {{{

  bool realspace_flag = _Data[0].getflag_inrealspace(); // keep realspace_flag of the operator we are zeroing

  while( ! _Data.empty() )
  {
    _Data.pop_back();
    _BandStencil.pop_back();
  }

  _Data.push_back( SmartField() );
  _Data[0].zero();
  _Data[0].setflag_inrealspace(realspace_flag);

  _BandStencil.push_back( 0 );
  _BW = 0;

} // }}}

void SmartFieldOp::setflag_inrealspace( bool realspaceflag )
{ // {{{
  for (UInt i=0; i<_Data.size(); i++)
  {
    _Data[i].setflag_inrealspace(realspaceflag);
  }
} // }}}

bool SmartFieldOp::getflag_inrealspace() const
{ // {{{

  bool rspace_flag = true;
  bool tmp_flag = true;

  for(UInt i=0; i<this->_Data.size(); ++i)
  {
    tmp_flag = this->_Data[i].getflag_inrealspace();
    if (i > 0 and rspace_flag != tmp_flag)
    {
      std::cout << "Warning: SmartFieldVec is mixed real space/fourier space. Use component functions.";
    }
    rspace_flag = tmp_flag;
  }

  return rspace_flag;

} // }}}

// ----------------- Reporters -----------------

void SmartFieldOp::ReportState()
{ // {{{

  std::cout << "Bandwidth = " << _BW << std::endl;
  std::cout << "Data Size = " << _Data.size() << std::endl;

  std::cout << "BandStencil" << std::endl;
  for (int i=0; i < _BandStencil.size(); i++) 
    std::cout << i << "\t";
  std::cout << "\n";
  for (int i=0; i < _BandStencil.size(); i++) 
    std::cout << _BandStencil[i] << "\t";
  std::cout << "\n\n";

  std::cout << "Data" << std::endl;
  std::cout << "\t";
  for (int j=0; j < _Data[0].GetSize(); j++) 
    std::cout << j << "\t";
  std::cout << "\n";
  for (int i=0; i < _Data.size(); i++)
  {
    std::cout << i << "\t";
    for (int j=0; j<_Data[0].GetSize(); j++ )
    {
      std::cout << "\t" << _Data[i].getelement(j);
    }
    std::cout << "\n";
  }

} // }}}

void SmartFieldOp::writefield( std::string filename, bool omitimaginary )
{ // {{{

  char buffer[80];
  std::string newname;

  for( int i=0; i<Nband(); i++ )
  {
    sprintf( buffer, "%i", _BandStencil[i]);
    newname = filename + "_" + buffer + ".dat";
    _Data[i].writefield( newname, omitimaginary );
  }
 
} // }}}

