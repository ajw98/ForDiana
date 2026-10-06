#!/usr/bin/env python
"""
Module for loading a field from file.
DRT, 19 Aug 2016

This module is broken out and specialized to handle the increased need 
for multiple file formats across all of the plotting and animation functions
that have been developed.
"""

import numpy as np
import sys
import struct
import os
import json

def jsonFile(filename):
  f = open(filename,"r");
  f1 = f.readlines();
  if (f1[0].find("#") == -1): return True;
  else: return False;

if (jsonFile("params_ICs.in")): 
  # JSON file format
  # Determine which file format
  params = open("params_ICs.in");
  jsonp = json.load(params)
  BinaryInFlag = int(jsonp["BinaryIn_flag"])
  BinaryOutFlag = int(jsonp["BinaryOut_flag"])
  params.close()
else:
  # Legacy file format
  # Determine which file format
  params = np.loadtxt("params_ICs.in");
  BinaryInFlag = int(params[0])
  BinaryOutFlag = int(params[1])

if (jsonFile("params_Grid.in")):
  # JSON file format
  # read in grid dimensions
  params = open('params_Grid.in');
  jsonp = json.load(params)
  GridFlag = int(jsonp["GridFlag"])
  Dim = int(jsonp["dim"])
  Lx = float(jsonp["Lx"])
  Ly = float(jsonp["Ly"])
  Lz = float(jsonp["Lz"])
  Nx = int(jsonp["Nx"])
  Ny = int(jsonp["Ny"])
  Nz = int(jsonp["Nz"])
else: 
  # Legacy file format
  # read in grid dimensions
  params = np.loadtxt('params_Grid.in');
  GridFlag = int(params[0])
  Dim = int(params[1])
  Lx = float(params[2])
  Ly = float(params[3])
  Lz = float(params[4])
  Nx = int(params[5])
  Ny = int(params[6])
  Nz = int(params[7])

# load the grid
if os.path.isfile("grid.dat"):
  grid = np.loadtxt("grid.dat")
  if (Dim == 1):
    X = grid[:,1];
    NPW = Nx;
    dx = float(Lx)/float(Nx)
  if (Dim == 2):
    X = grid[:, 1].reshape(Nx, Ny);
    Y = grid[:, 2].reshape(Nx, Ny);
    NPW = Nx*Ny;
    dx = float(Lx)/float(Nx)
    dy = float(Ly)/float(Ny)
  if (Dim == 3):
    X = grid[:, 1].reshape(Nx, Ny, Nz);
    Y = grid[:, 2].reshape(Nx, Ny, Nz);
    Z = grid[:, 3].reshape(Nx, Ny, Nz);
    NPW = Nx*Ny*Nz;
    dx = float(Lx)/float(Nx)
    dy = float(Ly)/float(Ny)
    dz = float(Lz)/float(Nz)
#elif GridFlag==0 # PS has grid info
else:
  print ('*** Error No grid information ***')

# load the time info
if (os.path.isfile("time.dat")):
  time = np.loadtxt("time.dat")
  if (len(time.shape)>1):
    t = time[:, 1]
  else:
    t = time;

def read(filename):
  
  name, ext = os.path.splitext(filename)
  if ( (BinaryOutFlag and (ext == '.out' or ext == '.dat'))
      or (BinaryInFlag and ext == '.in') ):

    infilehndl = open(filename, 'rb')
    infilehndl.seek(0)
    contents = infilehndl.read()

    # Check header and file version number
    header = struct.unpack_from("@8s",contents)
    if header[0] != b'FieldBin':
      sys.stderr.write("\nError: Not an unformatted Field file\n")
      sys.exit(1)

    pos = 9
    version = struct.unpack_from("@I",contents,offset=pos)
    pos = pos + 4
    #print ("Found unformatted field file with version {}".format(version[0]))
    if version[0] != 51:
      sys.stderr.write("\nError: Only version 51 is currently supported\n")
      sys.exit(1)

    nfields = struct.unpack_from("@i",contents,offset=pos)[0]
    pos = pos + 4
    #print ("# fields in file = {}".format(nfields))

    NDim = struct.unpack_from("@I",contents,offset=pos)[0]
    pos = pos + 4
    #print ("Spatial dimensionality = {}".format(NDim))

    griddim = struct.unpack_from("@{}L".format(NDim),contents,offset=pos)
    pos = pos + 8*NDim
    M = 1
    for i in range(NDim):
      M = M * griddim[i]
    #print ("PW grid : ",griddim,"\tTotal # PWs = ",M)

    (kspacedata,complexdata) = struct.unpack_from("@2?",contents,offset=pos)
    pos = pos + 2
    #print ("k space? {}\tComplex container? {}".format(kspacedata,complexdata))

    harray = struct.unpack_from("@{}d".format(NDim*NDim),contents,offset=pos)
    pos = pos + 8*NDim*NDim
    h = np.reshape(harray,(NDim,NDim))
    #print ("Cell tensor:\n",h)

    elsize = struct.unpack_from("@L",contents,offset=pos)[0]
    pos = pos + 8
    #print ("# bytes per element = {}".format(elsize))

    if elsize == 4 and not complexdata:
      #print (" * Single precision")
      fielddata = struct.unpack_from("@{}f".format(M*nfields),contents,offset=pos)
    elif elsize == 8 and complexdata:
      #print (" * Single precision")
      fielddata = struct.unpack_from("@{}f".format(2*M*nfields),contents,offset=pos)
    elif elsize == 8 and not complexdata:
      #print (" * Double precision")
      fielddata = struct.unpack_from("@{}d".format(M*nfields),contents,offset=pos)
    elif elsize == 16 and complexdata:
      #print (" * Double precision")
      fielddata = struct.unpack_from("@{}d".format(2*M*nfields),contents,offset=pos)
    else:
      sys.stderr.write("\nError: Unknown element size")
      sys.exit(1)

    # Because of how the data is read it, it appears it doesn't matter if
    # it is laid out like PS data or FD data.

    if (complexdata):
      # ignore complex part for now, only care about real fields
      data = np.array(fielddata).reshape(NPW,2)
      if Dim == 1:
        #return data[:,0] + 1j*data[:,1]
        return data[:,0]
      elif Dim == 2:
        #return data[:,0].reshape(Nx, Ny) + 1j*data[:,1].reshape(Nx, Ny)
        return data[:,0].reshape(Nx, Ny)
      elif Dim == 3:
        #return data[:,0].reshape(Nx, Ny, Nz) + 1j*data[:,1].reshape(Nx, Ny, Nz)
        return data[:,0].reshape(Nx, Ny, Nz)
    else:
      data = np.array(fielddata).reshape(NPW)
      if Dim == 1:
        return data
      elif Dim == 2:
        return data.reshape(Nx, Ny)
      elif Dim == 3:
        return data.reshape(Nx, Ny, Nz)

  else: # not BinaryOutFlag

    if (GridFlag==0):

      if (Dim == 1):
        data = np.loadtxt(filename)[:, 1];
        return data;
      elif (Dim == 2):
        data = np.loadtxt(filename)[:, 2];
        return data.reshape(Nx, Ny);
      elif (Dim == 3):
        data = np.loadtxt(filename)[:, 3];
        return data.reshape(Nx, Ny, Nz);
      else:
        print ('(Error) -- Reading PS data')

    else: # GridFlag == 1

      if (Dim == 1):
        data = np.zeros(X.shape);
        print ('(Error) Reading data for FD in 1D is not yet supported')
        return data;
      elif (Dim == 2):
        data = np.loadtxt(filename)[:, 1:].T;
        return data;
      elif (Dim == 3):
        data = np.loadtxt(filename)[:, 2:].T;
        return data.reshape(Nx, Ny, Nz);
      else:
        print ('(Error) -- Reading FD data')

