/***************************************************************
*
* Class to read grid parameters and initialize the grid
* 
* DRT -- Mon 16 May 2016
*
****************************************************************/

#include "Grid.h"

// ----------------- Constructor/Destructor -----------------

Grid :: Grid( pllhandler & pll, bool is_periodic ) :
_Pll( &pll ),
_IsPeriodic( is_periodic )
{ // {{{

  // Read in grid parameters

  std::string tmp_str;
  std::ifstream fid;

  if (jsonFile("params_Grid.in"))
  {
    fid.open("params_Grid.in");
    if ( fid.is_open() )
    {
      rapidjson::IStreamWrapper isw(fid);
      rapidjson::Document d;
      d.ParseStream(isw); 
      _GridFlag = d["GridFlag"].GetInt();

      _Dim = d["dim"].GetInt();

      _Lx = d["Lx"].GetDouble();

      _Ly = d["Ly"].GetDouble();

      _Lz = d["Lz"].GetDouble();

      _Nx = d["Nx"].GetInt();

      _Ny = d["Ny"].GetInt();

      _Nz = d["Nz"].GetInt();

      std::cout << std::endl;
      std::cout << " * Grid initiated\n";
      std::cout << "   - Grid Parameters:\n" << std::scientific;
      std::cout << "      GridFlag        = " << _GridFlag << "\n";
      std::cout << "      dim             = " << _Dim << "\n";
      std::cout << "      Lx              = " << _Lx << "\n";
      std::cout << "      Ly              = " << _Ly << "\n";
      std::cout << "      Lz              = " << _Lz << "\n";
      std::cout << "      Nx              = " << _Nx << "\n";
      std::cout << "      Ny              = " << _Ny << "\n";
      std::cout << "      Nz              = " << _Nz << "\n";

    }
    else
    {
      std::cout << "(Error): No Grid Parameters\n";
      exit(1);
    }
    fid.close();
  }
  
  else
  {//legacy format
    fid.open("params_Grid.in");
    if ( fid.is_open() )
    {
      std::getline(fid, tmp_str, '#');
      _GridFlag = atoi(tmp_str.c_str());
      std::getline(fid, tmp_str, '\n');

      std::getline(fid, tmp_str, '#');
      _Dim = atoi(tmp_str.c_str());
      std::getline(fid, tmp_str, '\n');

      std::getline(fid, tmp_str, '#');
      _Lx = atof(tmp_str.c_str());
      std::getline(fid, tmp_str, '\n');

      std::getline(fid, tmp_str, '#');
      _Ly = atof(tmp_str.c_str());
      std::getline(fid, tmp_str, '\n');

      std::getline(fid, tmp_str, '#');
      _Lz = atof(tmp_str.c_str());
      std::getline(fid, tmp_str, '\n');

      std::getline(fid, tmp_str, '#');
      _Nx = atoi(tmp_str.c_str());
      std::getline(fid, tmp_str, '\n');

      std::getline(fid, tmp_str, '#');
      _Ny = atoi(tmp_str.c_str());
      std::getline(fid, tmp_str, '\n');

      std::getline(fid, tmp_str, '#');
      _Nz = atoi(tmp_str.c_str());
      std::getline(fid, tmp_str, '\n');

      std::cout << std::endl;
      std::cout << " * Grid initiated\n";
      std::cout << "   - Grid Parameters:\n" << std::scientific;
      std::cout << "      GridFlag        = " << _GridFlag << "\n";
      std::cout << "      dim             = " << _Dim << "\n";
      std::cout << "      Lx              = " << _Lx << "\n";
      std::cout << "      Ly              = " << _Ly << "\n";
      std::cout << "      Lz              = " << _Lz << "\n";
      std::cout << "      Nx              = " << _Nx << "\n";
      std::cout << "      Ny              = " << _Ny << "\n";
      std::cout << "      Nz              = " << _Nz << "\n";
    }
    else
    {
      std::cout << "(Error): No Grid Parameters\n";
      exit(1);
    }
    fid.close();
  }


  if (_GridFlag == 0) // PS
  {

    _PSDim = _Dim;

    Vec2dReal cell_tensor(_Dim, Vec1dReal(_Dim, 0.)); // cell tensor
    cell_tensor[0][0] = _Lx;
    if (_Dim > 1) cell_tensor[1][1] = _Ly;
    if (_Dim > 2) cell_tensor[2][2] = _Lz;
    _Cell = new Cell(_Dim, cell_tensor, true);

    std::vector<UInt> Ngrid(_Dim, 0);
    Ngrid[0] = _Nx;
    if (_Dim > 1) Ngrid[1] = _Ny;
    if (_Dim > 2) Ngrid[2] = _Nz;
    _CurrLayout = new FFTlayout( *_Cell, Ngrid, *_Pll, false );

    _NPW = _CurrLayout->getNPWglobal();
    _FDSize = 1;

  }
  else // GridFlag == 1, FD
  {

    _PSDim = _Dim - 1;

    if (_PSDim < 1)
    {
      std::cout << std::endl << std::endl;
      std::cout << "Error: 1D simulations of finite difference are not supported at this time";
      std::cout << std::endl << std::endl;
      exit(1);
    }

    // set up a simulation cell in ND-1 dimensions
    Vec2dReal cell_tensor(_PSDim, Vec1dReal(_PSDim, 0.)); // cell tensor
    cell_tensor[0][0] = _Ly;
    if (_PSDim > 1) cell_tensor[1][1] = _Lz;
    _Cell = new Cell(_PSDim, cell_tensor, true);

    // set up FFT layout in ND-1 dimensions
    std::vector<UInt> Ngrid(_PSDim, 0);
    Ngrid[0] = _Ny;
    if (_PSDim > 1) Ngrid[1] = _Nz;
    _CurrLayout = new FFTlayout( *_Cell, Ngrid, *_Pll, false );

    _NPW = _CurrLayout->getNPWglobal();
    _FDSize = _Nx;

  }

  if (_IsPeriodic) // periodic
  { 
    _Dx = _Lx/_Nx;
  }
  else // not periodic
  { 
    _Dx = _Lx/(_Nx-1);
  }

  if (_GridFlag == 0) // PS
  {
    write_PS_grid();
  }
  else // FD
  {
    write_FD_grid();
  }

  _DxDim.push_back(_Dx);
  if (_Dim>0) _DxDim.push_back(_Ly/_Ny);
  if (_Dim>1) _DxDim.push_back(_Lz/_Nz);

} // }}}

Grid :: ~Grid()
{ // {{{

  delete _CurrLayout;
  delete _Cell;

} // }}}

// ----------------- Other useful functions -----------------

void Grid::write_PS_grid()
{ // {{{

  // Write grid.dat
  Vec1dReal xcart(_Dim);
  std::ofstream gridfile;
  gridfile.open("grid.dat");
  gridfile << "# Dimensions = " << _Dim << "\n" << std::scientific;
  gridfile << "# Nx = " << _Nx << ", Ny = " << _Ny << ", Nz = " << _Nz << "\n";
  gridfile << "# Lx = " << _Lx << ", Ly = " << _Ly << ", Lz = " << _Lz << "\n";
  switch (_Dim)
  {
    case 1: gridfile << "# Columns: n X kx\n";
            for ( int m = 0; m < _NPW; m++ )
            {
              xcart = _CurrLayout->rvec_of_indx_0corner(m, false);
              gridfile << m << "\t";
              gridfile << xcart[0] << "\t";
              gridfile << _CurrLayout->getKVecMap()[0][m] << "\n";
            }
            break;

    case 2: gridfile << "# Columns: n X Y kx ky\n";
            for ( int m = 0; m < _NPW; m++ )
            {
              xcart = _CurrLayout->rvec_of_indx_0corner(m, false);
              gridfile << m << "\t";
              gridfile << xcart[0] << "\t" << xcart[1] << "\t";
              gridfile << _CurrLayout->getKVecMap()[0][m] << "\t";
              gridfile << _CurrLayout->getKVecMap()[1][m] << "\n";
            }
            break;

    case 3: gridfile << "# Columns: n X Y Z kx ky kz\n";
            for ( int m = 0; m < _NPW; m++ )
            {
              xcart = _CurrLayout->rvec_of_indx_0corner(m, false);
              gridfile <<  m << "\t";
              gridfile << xcart[0] << "\t" << xcart[1] << "\t" << xcart[2] << "\t";
              gridfile << _CurrLayout->getKVecMap()[0][m] << "\t";
              gridfile << _CurrLayout->getKVecMap()[1][m] << "\t";
              gridfile << _CurrLayout->getKVecMap()[2][m] << "\n";
            }
            break;
  }
  gridfile.close();

} // }}}

void Grid::write_FD_grid()
{ // {{{

  // Write grid.dat
  Vec1dReal xcart(_PSDim);
  std::ofstream gridfile;
  gridfile.open("grid.dat");
  gridfile << "# Dimensions = " << _Dim << std::scientific;
  gridfile << ", (FD = 1, PS = " << _PSDim << ")\n";
  gridfile << "# Nx = " << _Nx << ", Ny = " << _Ny << ", Nz = " << _Nz << "\n";
  gridfile << "# Lx = " << _Lx << ", Ly = " << _Ly << ", Lz = " << _Lz << "\n";
  switch (_Dim)
  {
    case 1: gridfile << "# Columns: n X(FD)\n";
            for (int m = 0; m < _Nx; ++m)
            {
              gridfile << m << "\t";
              gridfile << m*_Dx << "\t";
            }
            break;

    case 2: gridfile << "# Columns: n X(FD) Y(PS) ky\n";
            for ( int m = 0; m < _Nx; m++ )
            {
              for ( int n = 0; n < _NPW; n++ )
              {
                xcart = _CurrLayout->rvec_of_indx_0corner(n, false);
                gridfile << m*_NPW+n << "\t";
                gridfile << m*_Dx << "\t" << xcart[0] << "\t";
                gridfile << _CurrLayout->getKVecMap()[0][n] << "\n";
              }
            }
              break;

    case 3: gridfile << "# Columns: n X(FD) Y(PS) Z(PS) ky kz\n";
            
            for ( int m = 0; m < _Nx; m++ )
            {
              for ( int n = 0; n < _NPW; n++ )
              {
                xcart = _CurrLayout->rvec_of_indx_0corner(n, false);
                gridfile << m*_NPW+n << "\t";
                gridfile << m*_Dx << "\t" << xcart[0] << "\t" << xcart[1] << "\t";
                gridfile << _CurrLayout->getKVecMap()[0][n] << "\t";
                gridfile << _CurrLayout->getKVecMap()[1][n] << "\n";
              }
            }
            break;
  }
  gridfile.close();

} // }}}

