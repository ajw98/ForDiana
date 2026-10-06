#!/usr/bin/env python3
"""
  Python script for plotting combined density and velocity fields

  Options:
    -h, --help: show help
    -b, --batch: use matplotlib "Agg" which works on braid, etc.
    -f, --half: plot half of the domain
    -n, --nfields: which fields to plot
    -s, --nosave: do not save the figure
    -v, --verbose: show the plot on-screen (doesn't work with -b)
    -c, --colormap: use different component volume fraction colormap
    -u, --vcolormap: use different velocity profile colormap
"""

# --- import modules and load command line options ---

import numpy as np
import sys
import os
import ReadField as field
import argparse
import json
# {{{
parser = argparse.ArgumentParser(description='Plot combined density and velocity fields')

parser.add_argument("filename", metavar="filename", type=str, nargs='+', 
  help="list of numbers of the density fields to plot (or 'in', 'out' or 'err')")

parser.add_argument("-b", "--batch", action="store_true",
  help='use for creating plots on clusters (negates -v)')

parser.add_argument("-f", "--half", action="store_true",
  help="plot only half of the domain (0, Ly/2)")

parser.add_argument("-l", "--slice", nargs=1, 
  help="plot a slice of the data with the form: 'x/y/z = VAL' (with spaces)")

parser.add_argument("-a", "--average", action="store_true", 
  help="used in conjuction with --slice. Makes the slice an average.")

DefArg = '12'
if os.path.isfile('params_TimeInt.in'):
  f1 = open('params_TimeInt.in',"r").readlines()
  if (f1[0].find("#") == -1):
    TimeIntFlag = json.load(open('params_TimeInt.in'))["TimeIntFlag"];
  else:
    TimeIntFlag = np.loadtxt('params_TimeInt.in')[0];#legacy file
  if (TimeIntFlag == 1):
    DefArg='12v';

parser.add_argument("-n", "--nfields", nargs=1, type=str,
  default=[DefArg],
  help="which frames frames to plot (e.g. '13' plots phi1/3, '124v' plots phi1/2/4 + velocity)")

parser.add_argument("-r", "--resolution", nargs=1, type=int,
  help="decrease the resolution of the vector field to every nth value")

parser.add_argument("-s", "--nosave", action="store_true",
  help="do not save the figure to file")

parser.add_argument("-t", "--streamlines", action="store_true",
  help="use streamlines instead of a quiver plot for the velocity")

parser.add_argument("-v", "--verbose", action="store_true",
  help="show the figures as created (doesn't work with -b)")

parser.add_argument("-c", "--colormap", nargs=1, type=str,
  default=['viridis'],
  help="use a different colormap for component volume fraction according to https://matplotlib.org/examples/color/colormaps_reference.html. Perceptually uniform sequential colormaps are recommended (viridis, plasma, inferno, magma)")

parser.add_argument("-u", "--vcolormap", nargs=1, type=str,
  default=['bone'],
  help="use a different colormap for velocity profile according to https://matplotlib.org/examples/color/colormaps_reference.html. Perceptually uniform sequential colormaps are recommended (viridis, plasma, inferno, magma)")

args = parser.parse_args();

if (args.batch):
  import matplotlib
  matplotlib.use('Agg')

import matplotlib.pyplot as plt

if (args.resolution != None):
  res = args.resolution[0]
else:
  res = 1

if (res < 0 or 
   (field.Dim == 2 and res > min(field.Nx, field.Ny)) or
   (field.Dim == 3 and res > min(field.Nx, field.Ny, field.Nz)) ):
  print ('(Error) -- resolution must be a positive integer ', \
        'less than min(Nx, Ny, Nz)')
  sys.exit();

if ( not args.nosave ):
  os.system('mkdir -p vel/')

if (field.Dim == 1):
  grid_dict = {'x' : field.X}
elif (field.Dim == 2):
  grid_dict = {'x' : field.X, 'y' : field.Y}
else: # field.Dim == 3
  grid_dict = {'x' : field.X, 'y' : field.Y, 'z' : field.Z}

if (args.slice != None):

  slice_dir = args.slice[0].split()[0]
  slice_val = float(args.slice[0].split()[-1])

  if (field.Dim == 1):
    print ('(Error) -- Cannot make a slice in 1D')
    sys.exit()

  if (field.Dim ==2 and not np.any(slice_dir == np.array(['x', 'y']))):
    print ("(Error) -- Slice direction must be either 'x' or 'y'")
    sys.exit()

  if (field.Dim == 3 and not np.any(slice_dir == np.array(['x', 'y', 'z']))):
    print ("(Error) -- Slice direction must be either 'x', 'y', or 'z'")
    sys.exit()

  slice_idx = np.where(grid_dict[slice_dir] == slice_val)

  if (len(grid_dict[slice_dir][slice_idx]) == 0):
    print ('(Error) -- No values stored at ', args.slice[0])
    sys.exit()

if (args.slice == None and args.average == True):
  print ('(Error) -- Can only use AVERAGE with SLICE')
  sys.exit()

if os.path.isfile('params_Keys.in'):
  f1 = open('params_TimeInt.in',"r").readlines()
  if (f1[0].find("#") == -1):
    NComp = json.load(open('params_Keys.in'))["Ncomp"];
  else:#legacy file format
    NComp = np.loadtxt('params_Keys.in')[0];
# }}}

# --- read nfields option ---

phi1_arg = 0
phi2_arg = 0
phi3_arg = 0
phi4_arg = 0
vel_arg = 0

for i in range(0, len(args.nfields[0])):
  if (args.nfields[0][i] == '1'):
    phi1_arg = 1
  elif (args.nfields[0][i] == '2'):
    phi2_arg = 1
  elif (args.nfields[0][i] == '3'):
    if (NComp >= 3):
      phi3_arg = 1
    else:
      print ('(Error) -- No 3rd component for 2-component system')
      sys.exit()
  elif (args.nfields[0][i] == '4'):
    if (NComp >= 4):
      phi4_arg = 1
    else:
      print ('(Error) -- No 4th component for 2- or 3-component system')
      sys.exit()
  elif (args.nfields[0][i] == 'v'):
    vel_arg = 1
  else:
    print ('(Error) -- Unable to read nfields option. Acceptable characters include 1,2,3,4,v')
    sys.exit()

# --- make plots ---

for i in range(0, len(args.filename)):

  # {{{
  if (args.filename[i] == "in"):
    if (NComp == 2):
      phi1_name = "phi1.in"
    elif (NComp == 3):
      phi1_name = "phi1.in"
      phi2_name = "phi2.in"
    elif (NComp == 4):
      phi1_name = "phi1.in"
      phi2_name = "phi2.in"
      phi3_name = "phi3.in"
    vx_name = "vx.in"
    vy_name = "vy.in"
    vz_name = "vz.in"
    plot_fname = "nfields_in"
  elif (args.filename[i] == "out"):
    if (NComp == 2):
      phi1_name = "phi1.out"
    elif (NComp == 3):
      phi1_name = "phi1.out"
      phi2_name = "phi2.out"
    elif (NComp == 4):
      phi1_name = "phi1.out"
      phi2_name = "phi2.out"
      phi3_name = "phi3.out"
    vx_name = "vx.out"
    vy_name = "vy.out"
    vz_name = "vz.out"
    plot_fname = "nfields_out"
  elif (args.filename[i] == "err"):
    if (NComp == 2):
      phi1_name = "phi1_ERR.dat"
    elif (NComp == 3):
      phi1_name = "phi1_ERR.dat"
      phi2_name = "phi2_ERR.dat"
    elif (NComp == 4):
      phi1_name = "phi1_ERR.dat"
      phi2_name = "phi2_ERR.dat"
      phi3_name = "phi3_ERR.dat"
    vx_name = "vx_ERR.dat"
    vy_name = "vy_ERR.dat"
    vz_name = "vz_ERR.dat"
    plot_fname = "nfields_err"
  else:
    if (NComp == 2):
      phi1_name = "phi1_%05d"%int(args.filename[i])+".dat"
    if (NComp == 3):
      phi1_name = "phi1_%05d"%int(args.filename[i])+".dat"
      phi2_name = "phi2_%05d"%int(args.filename[i])+".dat"
    if (NComp == 4):
      phi1_name = "phi1_%05d"%int(args.filename[i])+".dat"
      phi2_name = "phi2_%05d"%int(args.filename[i])+".dat"
      phi3_name = "phi3_%05d"%int(args.filename[i])+".dat"
    vx_name = "vx_%05d"%int(args.filename[i])+".dat"
    vy_name = "vy_%05d"%int(args.filename[i])+".dat"
    vz_name = "vz_%05d"%int(args.filename[i])+".dat"
    plot_fname = "nvel_%05d"%int(args.filename[i])

  print ("plotting: ", plot_fname)

  if (NComp == 2):
    if os.path.isfile(phi1_name):
      phi1 = field.read(phi1_name)
    else:
      print ('(Error) -- cannot open:', phi1_name)
      sys.exit()
    phi2 = 1 - phi1
  elif (NComp == 3):
    if os.path.isfile(phi1_name):
      phi1 = field.read(phi1_name)
    else:
      print ('(Error) -- cannot open:', phi1_name)
      sys.exit()
    if os.path.isfile(phi2_name):
      phi2 = field.read(phi2_name)
    else:
      print ('(Error) -- cannot open:', phi2_name)
      sys.exit()
    phi3 = 1 - phi1 - phi2
  elif (NComp == 4):
    if os.path.isfile(phi1_name):
      phi1 = field.read(phi1_name)
    else:
      print ('(Error) -- cannot open:', phi1_name)
      sys.exit()
    if os.path.isfile(phi2_name):
      phi2 = field.read(phi2_name)
    else:
      print ('(Error) -- cannot open:', phi2_name)
      sys.exit()
    if os.path.isfile(phi3_name):
      phi3 = field.read(phi3_name)
    else:
      print ('(Error) -- cannot open:', phi3_name)
      sys.exit()
    phi4 = 1 - phi1 - phi2 - phi3
  else:
    print ('(Error) -- systems with more than 4 components is not currently supported')
    sys.exit()

  # conditionally load the vx, vy files
  if (field.Dim > 1 and vel_arg == 1):
    if os.path.isfile(vx_name):
      vx = field.read(vx_name)
    else:
      print ('(Error) -- cannot open:', vx_name)
      sys.exit()

    if os.path.isfile(vy_name):
      vy = field.read(vy_name)
    else:
      print ('(Error) -- cannot open:', vy_name)
      sys.exit()

    if (field.Dim == 2):
      s = np.sqrt(vx**2 + vy**2)
    else: # field.Dim == 3
      if os.path.isfile(vz_name):
        vz = field.read(vz_name)
      else:
        print ('(Error) -- cannot open:', vz_name)
      s = np.sqrt(vx**2 + vy**2 + vz**2)
    max_s = max(np.amax(s),1e-10)

    if (field.Dim == 2):
      vel_dict = {'x' : vx, 'y' : vy}
    if (field.Dim == 3):
      vel_dict = {'x' : vx, 'y' : vy, 'z' : vz}
  # }}}

  # 1D plots
  # {{{
  if ((field.Dim == 1 and args.slice == None) or (field.Dim == 2 and args.slice != None)):

    # If doing a slice from 2D, get the index
    if (args.slice):

      if (slice_dir == 'x'):
        pltdir = 'y'
        x_len = field.Ly
      elif (slice_dir == 'y'):
        pltdir = 'x'
        x_len = field.Lx

      idx = slice_idx
      n_xtic = 5

    elif (args.half):

      pltdir = 'x'
      x_len = field.Lx/2
      idx = range(0, field.Nx/2)
      n_xtic = 5;

    else:

      pltdir = 'x'
      x_len = field.Lx
      idx = range(0, field.Nx)
      n_xtic = 5

    name_list = []
    phi_list = []
    fmt_list = []

    if (phi1_arg):
      name_list.append('polymer')
      phi_list.append(phi1)
      fmt_list.append('b-')
    if (phi2_arg):
      name_list.append('nonsolvent')
      phi_list.append(phi2)
      fmt_list.append('r--')
    if (phi3_arg):
      name_list.append('solvent')
      phi_list.append(phi3)
      fmt_list.append('g.-')
    if (phi4_arg):
      name_list.append('additive')
      phi_list.append(phi4)
      fmt_list.append('m^-')

    if (args.slice and args.average):
      for j in range(len(phi_list)):
        if (slice_dir == 'x'):
          phi_av = np.sum(phi_list[j], axis=0)/field.Nx
          phi_list[j] = np.tile(phi_av, (field.Nx, 1))
        elif (slice_dir == 'y'):
          phi_av = np.sum(phi_list[j], axis=1)/field.Ny
          phi_list[j] = np.tile(phi_av.reshape(field.Nx, 1), (1, field.Ny))

    plt.figure();
    for j in range(len(phi_list)):
      plt.plot(grid_dict[pltdir][idx], phi_list[j][idx], fmt_list[j], label=name_list[j]);

    plt.xlim([0, x_len])
    plt.ylim([0, 1])
    plt.xticks(np.linspace(0, x_len, n_xtic));
    plt.legend(loc=2)
    plt.xlabel("$x/R_{0}$", fontsize=10)
    plt.ylabel("volume fraction", fontsize=6)
    if (args.filename[i] == 'in' or args.filename[i] == 'out' or args.filename[i] == 'err'):
      plt.title(args.filename[i])
    else:
      plt.title("t = %f"%field.t[int(args.filename[i])])

  # }}}

  # 2D plots 
  if ( (field.Dim == 2 and args.slice == None)
       or (field.Dim == 3 and args.slice != None) ):
  # {{{

    if (args.slice):
      if (slice_dir == 'x'):
        x_dir = 'y'
        y_dir = 'z'
        x_lo = 0
        x_hi = field.Ly
        y_lo = 0
        y_hi = field.Lz
        nx = field.Ny
        ny = field.Nz
      elif (slice_dir == 'y'):
        x_dir = 'x'
        y_dir = 'z'
        x_lo = 0
        x_hi = field.Lx
        y_lo = 0
        y_hi = field.Lz
        nx = field.Nx
        ny = field.Nz
      elif (slice_dir == 'z'):
        x_dir = 'x'
        y_dir = 'y'
        x_lo = 0
        x_hi = field.Lx
        y_lo = 0
        y_hi = field.Ly
        nx = field.Nx
        ny = field.Ny

      idx = slice_idx
      X_plt = grid_dict[x_dir][idx].reshape(nx,ny)[::res, ::res].T
      Y_plt = grid_dict[y_dir][idx].reshape(nx,ny)[::res, ::res].T
      if (vel_arg):
        vx_plt = vel_dict[x_dir][idx].reshape(nx,ny)[::res, ::res].T
        vy_plt = vel_dict[y_dir][idx].reshape(nx,ny)[::res, ::res].T
        s_plt = s[idx].reshape(nx,ny)[::res, ::res].T

    elif (args.half):

      x_lo = 0.
      x_hi = field.Lx
      y_lo = field.Ly//2
      y_hi = field.Ly
      nx = field.Nx
      ny = field.Ny//2
      n_ytic = 5
      idx = np.meshgrid( np.arange(0, field.Nx), 
                         np.arange(field.Ny//2, field.Ny), 
                         indexing='ij' )

      X_plt = field.X[::res, field.Ny//2::res].T
      Y_plt = field.Y[::res, field.Ny//2::res].T
      if (vel_arg):
        vx_plt = vx[::res, field.Ny//2::res].T
        vy_plt = vy[::res, field.Ny//2::res].T
        s_plt = s[::res, field.Ny//2::res].T

    else:

      x_lo = 0.
      x_hi = field.Lx
      y_lo = 0.
      y_hi = field.Ly
      nx = field.Nx
      ny = field.Ny
      n_ytic = 5
      idx = np.meshgrid( np.arange(0, field.Nx), 
                         np.arange(0, field.Ny), 
                         indexing='ij' )

      X_plt = field.X[::res, ::res].T
      Y_plt = field.Y[::res, ::res].T
      if (vel_arg):
        vx_plt = vx[::res, ::res].T
        vy_plt = vy[::res, ::res].T
        s_plt = s[::res, ::res].T

    name_list = []
    phi_list = []

    if (phi1_arg):
      name_list.append('polymer')
      phi_list.append(phi1)
    if (phi2_arg):
      name_list.append('nonsolvent')
      phi_list.append(phi2)
    if (phi3_arg):
      name_list.append('solvent')
      phi_list.append(phi3)
    if (phi4_arg):
      name_list.append('additive')
      phi_list.append(phi4)
    if (vel_arg):
      name_list.append('velocity')
    
    if (len(args.nfields[0]) < 4):
      plt_shape = [1,len(args.nfields[0])]
    else:
      plt_shape = [2,len(args.nfields[0])-2]

    if (args.slice and args.average):
      print ('Sorry, averaged slicing is not available for 2D slices')
      sys.exit()
      
    plt.figure(figsize=(1.725*plt_shape[1], 1.725*plt_shape[0]));
    #print(plt_shape)
    n = 0
    for j in range(plt_shape[0]): 
      for k in range(plt_shape[1]):

        if ( n < len(phi_list) ):
          print('1')
          plt.subplot(plt_shape[0], plt_shape[1], n+1)
          plt.imshow( phi_list[n][tuple(idx)].reshape(nx,ny).T, 
            origin='lower', aspect='equal', \
            extent=[x_lo, x_hi, y_lo, y_hi], \
            interpolation='bicubic', cmap=args.colormap[0], vmin=0, vmax=1)
            # other interpolation keywords: bicubic, gaussian, nearest
        elif ( n > len(name_list)-1 ):
        	continue
        else:
          print('2')
          plt.subplot(plt_shape[0], plt_shape[1], n+1)
          if (args.streamlines):
            plt.streamplot(X_plt, Y_plt, vx_plt, vy_plt, \
                            linewidth=0.8, color='gray', density=[1, 4])
            plt.imshow(s_plt, origin='lower', aspect='equal', \
                      extent=[x_lo, x_hi, y_lo, y_hi], \
                      interpolation='bicubic', cmap=args.vcolormap[0], vmin=0, vmax=max_s)
          else:
            vmax = max(np.amax(np.sqrt(vx_plt**2 + vy_plt**2)), 1e-6)
            X_sub  = X_plt[::64, ::64]
            Y_sub  = Y_plt[::64, ::64]
            vx_sub = vx_plt[::64, ::64]
            vy_sub = vy_plt[::64, ::64]
            s_sub = s_plt[::64,::64]
            vmax = 0.015
            plt.quiver(X_sub, Y_sub, vx_sub, vy_sub, s_sub, 
                      pivot='mid', angles='xy', scale=2.*vmax, units='inches', 
                      width=0.02, headwidth=1, minlength=1e-8, clim=[0,vmax], \
                      cmap=args.vcolormap[0])

            plt.xlim([x_lo, x_hi])
            plt.ylim([y_lo, y_hi])
            plt.gca().set_aspect('equal')

        #plt.colorbar(fraction=0.046, pad=0.1)
        #plt.xlabel("$x/R_{0}$", fontsize=6,labelpad=1)
        #plt.ylabel("$y/R_{0}$", fontsize=6,labelpad = 0.5)
        plt.xticks(np.linspace(0,field.Lx,3),fontsize=6);
        plt.yticks([512,768,1024],[0,256,512],fontsize=6);
        #plt.title(" ")#name_list[n])
        
        n += 1

    if (args.filename[i] == 'in' or args.filename[i] == 'out' \
         or args.filename[i] == 'err'):
      plt.suptitle(args.filename[i], fontsize=6)
    else:
      plt.suptitle(" ")#"t = %f"%field.t[int(args.filename[i])], fontsize=6)
    fig = plt.gcf()
    #ax.set_aspect('equal')
    for ax in fig.axes:
        print('ax')
        ax.set_position([0.25, 0.2, 0.55, 0.55])
        for text in ax.findobj(match=plt.Text):
            text.set_fontsize(6)

    for text in fig.texts:
        text.set_fontsize(6)

    fig.canvas.draw()  # Make sure the layout is finalized

    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())

    width, height = bbox.width, bbox.height
    print(f"Graph width: {width:.3f} inches")
    print(f"Graph height: {height:.3f} inches")

    print(fig.get_size_inches())


    plt.tight_layout()
    # }}}

  if (field.Dim == 3 and args.slice == None):
    print ('Sorry, non-sliced 3D plots are not supported at this time.')

  if (not args.nosave):
    #fig, ax = plt.subplots(figsize=(1.625, 1.625), constrained_layout=True)

    plt.savefig("vel/"+plot_fname+".png", dpi=600)

  if (not args.verbose):
    plt.close()

if (args.verbose):
  plt.show();

