#!/usr/bin/env python3

import numpy as np
import sys
import os
import struct
import json

import os

if (len(sys.argv) == 2):
  runid = int(sys.argv[1]); # should be an integer
else:
  runid = 0 # if nothing supplied, run it for zero

print('Writing initial conditions for runid: ', runid )

# ----------------------------------------------
#  Initial Conditions
# ----------------------------------------------
# BinaryIn_flag (==0, ascii input files, ==1 binary input files)
# BinaryOut_flag (==0, ascii output files, ==1 binary output files)
# PS_flag (==1, files formatted for PS, ==0 files formatted for FD)
# WriteVel_flag (==0, don't write vel files, ==1 write vel files)

# Note: Binary files are significantly smaller and faster. 
# Use them unless you have a good reason otherwise.
BinaryIn_flag=1
BinaryOut_flag=1

f = open('params_Grid.in')
grid_params = json.load(f)

# start params_IC.in
# {{{
PS_flag=int(grid_params["GridFlag"])

WriteVel_flag = 0
BodyForceFlag = 0
if (os.path.isfile('params_TimeInt.in')):
  g = open('params_TimeInt.in')
  TimeInt_params = json.load(g)
  TimeIntFlag = TimeInt_params["TimeIntFlag"]
  if (TimeIntFlag == 1 or TimeIntFlag == 2 or TimeIntFlag == 4):
    WriteVel_flag = 1;
    BodyForceFlag = TimeInt_params["BodyForce_flag"]

print('  - Binary Input?', BinaryIn_flag==1)
print('  - Binary Output?', BinaryOut_flag==1)

# i/o
#f1 = open('params_ICs.in', 'w')
json_dict = {'BinaryIn_flag':BinaryIn_flag, 
                  'BinaryOut_flag':BinaryOut_flag}

# Number of components, other keys
keys = open('params_Keys.in')
jsonkeys = json.load(keys)
Ncomp = int(jsonkeys["Ncomp"])
NIcomp = int(jsonkeys["NIcomp"])
# }}}

# ------------------------
# set up grid
# ------------------------
# {{{
dim = int(grid_params["dim"])
Lx = grid_params["Lx"]
Ly = grid_params["Ly"]
Lz = grid_params["Lz"]
Nx = int(grid_params["Nx"])
Ny = int(grid_params["Ny"])
Nz = int(grid_params["Nz"])

if (dim == 1): Ngrid = Nx
if (dim == 2): Ngrid = Nx*Ny
if (dim == 3): Ngrid = Nx*Ny*Nz

if (dim >= 1):
  x = np.linspace(0, Lx, Nx+1)[0:-1]
  X = x
if (dim >= 2):
  y = np.linspace(0, Ly, Ny+1)[0:-1]
  (X, Y) = np.meshgrid(x, y)
  X = X.T
  Y = Y.T
if (dim == 3):
  z = np.linspace(0, Lz, Nz+1)[0:-1]
  (X, Y, Z) = np.meshgrid(x, y, z)
  X = np.transpose(X, axes=(1,0,2))
  Y = np.transpose(Y, axes=(1,0,2))
  Z = np.transpose(Z, axes=(1,0,2))
# }}}

# ----------------------------------------------
#  Pick type of concentration IC
# ----------------------------------------------
# (000 Series) Ternary NIPS
#   0: homogeneous conditions
#   1: uniform noise
#   2: guassian noise
#   3: tanh profile in y
#   4: sinusoids in x and y
#   5: concentration gradient in y
#   6: 1D tanh profile in x
#   7: tanh profile in y with roll cells 
#   8: read from file, add gradient
#   9: 1D sinusoid in x or y
#  10: gaussian pulse in 2D
#  11: 3D sinusoid
#  12: Circular droplet in 2D
#  13: Multiple spheres at random locations in 3D (**TODO: Make this periodic **)
#  14: Retracting Droplet
#  15: Interfacial Droplet
#  16: 2D tanh profile in x (Pseudo-1D VIPS Dirichlet BC)
#  17: Multiple spheres at random locations in 2D
# (100 Series) Quaternary NIPS
# 100: homogeneous conditions
# 101: uniform noise
# 102: guassian noise
# 103: tanh profile in y
# 104: sinusoids in x and y
# 105: concentration gradient in y
# 106: 1D tanh profile in x
# 107: tanh profile in y with roll cells
# 108: read from file, add gradient
# 109: 1D sinusoid in x or y
# 110: gaussian pulse in 2D
# 111: 3D sinusoid
# 112: Circular droplet in 2D
# (200 Series) Block Polymers
# 200: homogeneous conditions
# 201: uniform noise
# 202: guassian noise (**Not written**)
# 203: tanh profile in y (**Not written**)
# 204: sinusoids in x and y (**Not written**)
# 205: concentration gradient in y (**Not written**)
# 206: 1D tanh profile in x
# 207: tanh profile in y with roll cells (**Not written**)
# 208: read from file, add velocity perturbation
# 209: 1D sinusoid in x or y
# 210: gaussian pulse in 2D (**Not written**)
# 211: 3D sinusoid (**Not written**)
# 212: Circular droplet in 2D


IC_flag =3 

Y=X

print('  - Creating initial condition: ', IC_flag)

phi = np.zeros((NIcomp, *X.shape))

# --- (000 Series) Ternary NIPS ---
if (IC_flag == 0): # homogeneous initial condition
# {{{
  # params
  phi1_0 = 0.1
  phi2_0 = 0.1

  # phi
  phi[0] = np.zeros(X.shape) + phi1_0
  phi[1] = np.zeros(X.shape) + phi2_0

  # i/o
  ic_dict = {'IC_flag':IC_flag, \
             'phi1_0':phi1_0, \
             'phi2_0':phi2_0}
# }}}
elif (IC_flag == 1): # uniform noise about homogeneous initial conditions
# {{{
  # params
  phi1_0 = 0.1 #0.46
  phi2_0 = 0.1 #0.44 #0.45 + 0.05*runid
  sig = 1e-2 # std. deviation for random noise
  my_seed = 1234 #np.random.randint(0, 4294967295) #1234

  # phi
  np.random.seed(my_seed)
  phi1 = sig*np.random.rand(*X.shape)
  phi2 = sig*np.random.rand(*X.shape)
  phi1_k = np.fft.fftn(phi1)
  phi2_k = np.fft.fftn(phi2)
  phi1_k[np.unravel_index(0, X.shape)] = phi1_0*Ngrid # set 0 element to mean
  phi2_k[np.unravel_index(0, X.shape)] = phi2_0*Ngrid
  phi[0] = np.fft.ifftn(phi1_k).real
  phi[1] = np.fft.ifftn(phi2_k).real

  # i/o
  ic_dict = {'IC_flag':IC_flag, \
             'phi1_0':phi1_0, \
             'phi2_0':phi2_0, \
             'sig':sig, \
             'my_seed':my_seed}
# }}}
elif (IC_flag == 2): # guassian random initial conditions
# {{{
  # params
  phi1_0 = 0.2
  phi2_0 = 0.6
  sig = 1e-2 # std. deviation for random noise (IC_flag == 1, 2, 4 only)
  my_seed = 1234 #np.random.randint(0, sys.maxint) #1234

  # mask for high k-modes
  dx = Lx/Nx;
  dy = Ly/Ny;
  kx = 2*np.pi*np.fft.fftfreq(X.shape[0])/dx
  ky = 2*np.pi*np.fft.fftfreq(X.shape[1])/dy
  Nyquist_kx = np.pi/dx; # maximum wavevector in x
  Nyquist_ky = np.pi/dy; # maximum wavevector in y
  x_mask = np.less_equal(abs(kx), Nyquist_kx/2.);
  y_mask = np.less_equal(abs(ky), Nyquist_ky/2.);
  mask = np.outer(x_mask, y_mask);

  # phi
  np.random.seed(my_seed)
  phi1 = sig*np.random.randn(*X.shape)
  phi2 = sig*np.random.randn(*X.shape)
  phi1_k = np.fft.fftn(phi1)
  phi2_k = np.fft.fftn(phi2)
  phi1_k[np.unravel_index(0, X.shape)] = phi1_0*Ngrid
  phi2_k[np.unravel_index(0, X.shape)] = phi2_0*Ngrid
  phi1_k = phi1_k*mask;
  phi2_k = phi2_k*mask;
  phi[0] = np.fft.ifftn(phi1_k).real
  phi[1] = np.fft.ifftn(phi2_k).real

  # i/o
  ic_dict = {'IC_flag':IC_flag, \
             'phi1_0':phi1_0, \
             'phi2_0':phi2_0, \
             'sig':sig, \
             'my_seed':my_seed}
# }}}
elif (IC_flag == 3): # tanh profile in y
# {{{

  # ** Read from phase diagram **
  #eq_phi = np.loadtxt("phase_diagram_bin.in");
  #alpha = eq_phi[:, :3]
  #beta  = eq_phi[::-1, 3:]
  #phi1_in = alpha[60, 0]
  #phi2_in = alpha[60, 1]
  #phi1_out = beta[60, 0]
  #phi2_out = beta[60, 1]
#0.13562690857598853 0.36007543905275574 0.5042976523712558
  #2 0.15347631604319042 0.46633331821008717 0.3801903657467224 1.0 0.0037177086299875536 0.9800000000018911 0.0149999999987398 0.9987177086306185
  #params
  phi1_in = 0.3 #0.1534763160431904 #polymer volume fraction in inner region
  phi2_in = 0.68#0.46633331821008717 # non-solvent volume fraction in inner region
  phi1_out =  0.005#0.0037177086299875536 # polymer volume fraction in outer region
  phi2_out =  0.99#0.9800000000018911 # non-solvent volue fraction in outer region
  sig = 1e-3; # noise strength
  w = .9 # boundary thickness
  f = 1./2. # fraction of inner region area/total area
  my_seed = np.random.randint(0, 4294967295) #4294... is the maximum seed value

  # basic phi step parallel to the x-axis
  lb = 0.5 - 0.5*f; # lower bound
  ub = 0.5 + 0.5*f; # upper bound
  Y_lb = (Y-lb*Lx)/w#(Y-lb*Ly)/w
  Y_ub = (ub*Lx -Y)/w#(ub*Ly-Y)/w
  dphi1 = phi1_in - phi1_out
  dphi2 = phi2_in - phi2_out
  np.random.seed(my_seed)
  phi1_0 = 0.5*dphi1*np.tanh(Y_lb) + 0.5*dphi1*np.tanh(Y_ub) + phi1_out 
  phi2_0 = 0.5*dphi2*np.tanh(Y_lb) + 0.5*dphi2*np.tanh(Y_ub) + phi2_out

  # get random perturbation
  phi1_rand = sig*(np.random.rand(*Y.shape)-0.5)
  phi2_rand = sig*(np.random.rand(*Y.shape)-0.5)

  # mask for high k-modes
  dx = float(Lx)/Nx;
  dy = float(Ly)/Ny;
  kx = 2*np.pi*np.fft.fftfreq(X.shape[0])/dx
  ky = 2*np.pi*np.fft.fftfreq(X.shape[0])/dy
  Nyquist_kx = np.pi/dx; # maximum wavevector in x
  Nyquist_ky = np.pi/dy; # maximum wavevector in y
  x_mask = np.less_equal(abs(kx), Nyquist_kx/2.);
  y_mask = np.less_equal(abs(ky), Nyquist_ky/2.);
  mask = np.outer(x_mask, y_mask);

  # apply low-pass filter to noise
  phi1_rand_k = np.fft.fftn(phi1_rand)
  phi2_rand_k = np.fft.fftn(phi2_rand)
  phi1_rand_k = phi1_rand_k*mask;
  phi2_rand_k = phi2_rand_k*mask;
  phi1_rand = np.fft.ifftn(phi1_rand_k).real
  phi2_rand = np.fft.ifftn(phi2_rand_k).real
  # ensure phi1 and phi2 are symmetric in y after adding noise
  phi1_rand = 0.5*(phi1_rand + phi1_rand[:, ::-1])
  phi2_rand = 0.5*(phi2_rand + phi2_rand[:, ::-1])

  phi[0] = phi1_0 #+phi1_rand
  phi[1] = phi2_0 #+phi2_rand

  # i/o
  ic_dict = {'IC_flag':IC_flag, \
             'phi1_in':phi1_in, \
             'phi2_in':phi2_in, \
             'phi1_out':phi1_out, \
             'phi2_out':phi2_out, \
             'sig':sig, \
             'w':w, \
             'f':f, \
             'my_seed':my_seed}
# }}}
elif (IC_flag == 4): # sinusoids in x and y
# {{{
  #params
  phi1_0 = 0.25 # avg concentration of polymer
  phi2_0 = 0.60 # avg concentration of non-solvent
  Tx = 1 # periods for sinusoid in x-dir
  Ty = 3 # periods for sinusoid in y-dir
  A = 0.05 # sinusoid amplitude

  # phi
  phi[0] = A*np.sin(2*Tx*np.pi*X/Lx)*np.sin(2*Ty*np.pi*Y/Ly) + phi1_0;
  phi[1] = A*np.sin(2*Tx*np.pi*X/Lx)*np.sin(2*Ty*np.pi*Y/Ly) + phi2_0;

  # i/o
  ic_dict = {'IC_flag':IC_flag, \
             'phi1_0':phi1_0, \
             'phi2_0':phi2_0, \
             'Tx':Tx, \
             'Ty':Ty, \
             'A':A}
# }}}
elif (IC_flag == 5): # gradient in y-direction
# {{{
  #params
  phi1_N = 0.30 # concentration of polymer on north end
  phi1_S = 0.30 # concentration of polymer on south end
  phi2_N = 0.65 # concentration of non-solvent on north end
  phi2_S = 0.50 # concentration of non-solvent on south end
  sig = 5e-3 # std. deviation for random noise
  my_seed = np.random.randint(0, 4294967295) #1234

  # phi
  np.random.seed(my_seed)
  phi1 = sig*np.random.rand(*X.shape)
  phi2 = sig*np.random.rand(*X.shape)
  phi1_k = np.fft.fftn(phi1)
  phi2_k = np.fft.fftn(phi2)
  phi1_k[np.unravel_index(0, X.shape)] = 0 # set mean to 0
  phi2_k[np.unravel_index(0, X.shape)] = 0
  phi1 = np.fft.ifftn(phi1_k).real
  phi2 = np.fft.ifftn(phi2_k).real

  half_Ny = int(Ny/2)

  phi[0, :, :half_Ny] = phi1[:, :half_Ny] - 2*(phi1_N-phi1_S)*Y[:, :half_Ny]/Ly + phi1_N;
  phi[1, :, :half_Ny] = phi2[:, :half_Ny] - 2*(phi2_N-phi2_S)*Y[:, :half_Ny]/Ly + phi2_N;
  phi[0, :, half_Ny:] = phi1[:, half_Ny:] + 2*(phi1_N-phi1_S)*(Y[:, half_Ny:]/Ly-0.5) + phi1_S;
  phi[1, :, half_Ny:] = phi2[:, half_Ny:] + 2*(phi2_N-phi2_S)*(Y[:, half_Ny:]/Ly-0.5) + phi2_S;

  # i/o
  ic_dict = {'IC_flag':IC_flag, \
             'phi1_N':phi1_N, \
             'phi1_S':phi1_S, \
             'phi2_N':phi2_N, \
             'phi2_S':phi2_S, \
             'sig':sig, \
             'my_seed':my_seed}
# }}}
elif (IC_flag == 6): # 1D tanh profile in x
# {{{

  #phi_params = np.loadtxt("../../setup/run"+str(runid)+"/params_Equilibrium.in")
#params
  phi1_in = 0.1 #phi_params[0] # polymer volume fraction in inner region
  phi2_in = 0.01 #phi_params[1] # non-solvent volume fraction in inner region
  phi1_out = 0.01 #phi_params[2] # polymer volume fraction in outer region
  phi2_out = 0.9 #phi_params[3] # non-solvent volue fraction in outer region
  phi1_part = 0.1
  phi2_part = 0.01
  f = 0.5 #phi_params[4] #1./2. # fraction of inner region area/total area
  sig = 0.5e-3 #0.5e-3 # noise strength
  w = 0.5 # boundary thickness
  my_seed = np.random.randint(0, 4294967295) #4294... is the maximum seed value
  
  R  = 0.1*Lx
  Rx = R          # droplet x-axis
  Ry = R          # droplet y-axis
  Y_c = (X - Lx/2)
#  Z_c = (Z - Lz/2)
  
  # mask for high k-modes
  dx = float(Lx)/Nx;
  kx = 2*np.pi*np.fft.fftfreq(X.shape[0])/dx
  Nyquist_kx = np.pi/dx; # maximum wavevector in x
  x_mask = np.less_equal(abs(kx), Nyquist_kx/2.);
  
  # phi
  lb = 0.5 - 0.5*f; # lower bound
  ub = 0.5 + 0.5*f; # upper bound
  X_lb = (X-lb*Lx)
  X_ub = (ub*Lx-X)
  X_c = X_lb - 1/2*R
  dphi1 = phi1_in - phi1_out
  dphi2 = phi2_in - phi2_out
  dphi1p = phi1_part - phi1_out
  dphi2p = phi2_out - phi2_part
  np.random.seed(my_seed)
  phi1 = -0.5*dphi1p*(np.tanh((X_c**2/R**2 + Y_c**2/R**2  - 1)/w) + 1) + phi1_part
  phi2_part = -0.5*dphi2p*(np.tanh((X_c**2/R**2 + Y_c**2/R**2 - 1)/w) + 1) + phi2_out
  phi2 = 0.5*dphi2*(np.tanh(X_ub/w)+np.tanh(X_lb/w))*(0.5*dphi1p*(np.tanh((X_c**2/R**2 + Y_c**2/R**2 - 1)/w) + 1)) + phi2_out
  phi2 = phi2-phi2_part#phi1 = 0.5*dphi1*np.tanh(X_lb) + 0.5*dphi1*np.tanh(X_ub) + phi1_out + sig*(np.random.rand(*X.shape)-0.5)
  #phi2 = 0.5*dphi2*np.tanh(X_lb) + 0.5*dphi2*np.tanh(X_ub) + phi2_out + sig*(np.random.rand(*X.shape)-0.5)
  #phi1 += 0.5*dphi1p*(np.tanh((X_lb**2 + Y_c**2 - 1)/w) + 1) + phi1_part
  #phi2 += 0.5*dphi2p*(np.tanh((X_lb**2 + Y_c**2 - 1)/w) + 1) + phi2_part

  # apply low-pass filter
#  phi1_k = np.fft.fft(phi1)
#  phi2_k = np.fft.fft(phi2)
#  phi1_k = phi1_k*x_mask;
#  phi2_k = phi2_k*x_mask;
#  phi1 = np.fft.ifft(phi1_k).real
#  phi2 = np.fft.ifft(phi2_k).real

  # ensure phi1 and phi2 are symmetric in x after adding noise
  phi[0] = phi1#0.5*(phi1 + phi1[::-1])
  phi[1] = phi2#0.5*(phi2 + phi2[::-1])

  # i/o
  ic_dict = {'IC_flag':IC_flag, \
             'phi1_in':phi1_in, \
             'phi2_in':phi2_in, \
             'phi1_out':phi1_out, \
             'phi2_out':phi2_out, \
             'sig':sig, \
             'w':w, \
             'f':f, \
             'my_seed':my_seed}
# }}}
elif (IC_flag == 7): # tanh profile in y with roll cells
# {{{
  #params
  phi1_in =  0.12 #polymer volume fraction in inner region
  phi2_in = 0.875 # non-solvent volume fraction in inner region
  phi1_out = 0.005 # polymer volume fraction in outer region
  phi2_out =  0.99 # non-solvent in outer region
  sig = 9.9e-3 #0.99e-2 # noise strength
  w = 4 # boundary thickness
  f = 1./4. #5./12. # fraction of inner region area/total area
  my_seed = np.random.randint(0, 4294967295) #4294... is the maximum seed value

  # basic phi step parallel to the x-axis
  lb = 0.5 - 0.5*f; # lower bound
  ub = 0.5 + 0.5*f; # upper bound
  Y_lb = (Y-lb*Ly)/w
  Y_ub = (ub*Ly-Y)/w
  dphi1 = phi1_in - phi1_out
  dphi2 = phi2_in - phi2_out
  np.random.seed(my_seed)
  phi1_0 = 0.5*dphi1*np.tanh(Y_lb) + 0.5*dphi1*np.tanh(Y_ub) + phi1_out 
  phi2_0 = 0.5*dphi2*np.tanh(Y_lb) + 0.5*dphi2*np.tanh(Y_ub) + phi2_out

  # get random perturbation
  phi1_rand = sig*(np.random.rand(*Y.shape)-0.5)
  phi2_rand = sig*(np.random.rand(*Y.shape)-0.5)

  # mask for high k-modes
  dx = float(Lx)/Nx;
  dy = float(Ly)/Ny;
  kx = 2*np.pi*np.fft.fftfreq(X.shape[0])/dx
  ky = 2*np.pi*np.fft.fftfreq(X.shape[1])/dy
  Nyquist_kx = np.pi/dx; # maximum wavevector in x
  Nyquist_ky = np.pi/dy; # maximum wavevector in y
  x_mask = np.less_equal(abs(kx), Nyquist_kx/2.);
  y_mask = np.less_equal(abs(ky), Nyquist_ky/2.);
  mask = np.outer(x_mask, y_mask);

  # apply low-pass filter to noise
  phi1_rand_k = np.fft.fftn(phi1_rand)
  phi2_rand_k = np.fft.fftn(phi2_rand)
  phi1_rand_k = phi1_rand_k*mask;
  phi2_rand_k = phi2_rand_k*mask;
  phi1_rand = np.fft.ifftn(phi1_rand_k).real
  phi2_rand = np.fft.ifftn(phi2_rand_k).real
  # ensure phi1 and phi2 are symmetric in y after adding noise
  phi1_rand = 0.5*(phi1_rand + phi1_rand[:, ::-1])
  phi2_rand = 0.5*(phi2_rand + phi2_rand[:, ::-1])

  # apply a velocity perturbation to advect phi one step
  A=0.5; # amplitude of velocity perturbation
  lam_x = 2.*Lx/(np.linspace(1,8,8)[runid]);
  lam_y = Ly/4.;
  xc = 0. #Lx/2.;
  yc = 0. #Ly/2.;

  #U =  1j*Ky[n,m]*np.exp(1j*Kx[n,m]*X+1j*Ky[n,m]*Y)
  #V = -1j*Kx[n,m]*np.exp(1j*Kx[n,m]*X+1j*Ky[n,m]*Y)

  # in x and y
  vx = A*np.cos(2*np.pi*(X-xc)/lam_x)*np.cos(2*np.pi*(Y-yc)/lam_y)
  vy = A*np.sin(2*np.pi*(X-xc)/lam_x)*np.sin(2*np.pi*(Y-yc)/lam_y)
  vz = A*np.zeros(X.shape)

  # in x only
  #vx = A*np.cos(2*np.pi*(X-xc)/lam_x)
  #vy = A*np.zeros(Y.shape)
  #vz = A*np.zeros(X.shape)

  # now, advect phi_1 and phi_2 one step
  (Kx, Ky) = np.meshgrid(kx, ky, indexing='ij')
  phi1_k = np.fft.fftn(phi1_0)
  phi2_k = np.fft.fftn(phi2_0)
  phi1_adv = - vx*np.fft.ifftn(1j*Kx*phi1_k).real - vy*np.fft.ifftn(1j*Ky*phi1_k).real
  phi2_adv = - vx*np.fft.ifftn(1j*Kx*phi2_k).real - vy*np.fft.ifftn(1j*Ky*phi2_k).real

  #phi1 = phi1_0+phi1_rand+phi1_adv
  #phi2 = phi2_0+phi2_rand+phi2_adv
  phi[0] = phi1_0+phi1_adv
  phi[1] = phi2_0+phi2_adv

  # i/o
  ic_dict = {'IC_flag':IC_flag, \
             'phi1_in':phi1_in, \
             'phi2_in':phi2_in, \
             'phi1_out':phi1_out, \
             'phi2_out':phi2_out, \
             'sig':sig, \
             'w':w, \
             'f':f, \
             'lam_x':lam_x, \
             'lam_y':lam_y, \
             'my_seed':my_seed}
# }}}
elif (IC_flag == 8): # read from file, add gradient
# {{{

  eq_dim=1
  if (eq_dim == 1): # (for 1D equilibrium result)

    #phi1_file='../equil/phi1.out'
    #phi2_file='../equil/phi2.out'
    phi1_file='equil/phi1.out'
    phi2_file='equil/phi2.out'
  
    #phi1_1D = np.loadtxt(phi1_file)[:,1];
    #phi2_1D = np.loadtxt(phi2_file)[:,1];
    phi1_1D = np.load(phi1_file)[:,1];
    phi2_1D = np.load(phi2_file)[:,1];
    phi1_0 = np.tile(phi1_1D, (Nx, 1))
    phi2_0 = np.tile(phi2_1D, (Nx, 1))

  elif (eq_dim == 2): # (for 2D equilibrium result)

    #phi1_file='../equil/phi1.out'
    #phi2_file='../equil/phi2.out'

    phi1_file='equil/phi1.out'
    phi2_file='equil/phi2.out'

    #tmp = np.loadtxt(phi1_file);
    tmp = np.load(phi1_file);
    #X = tmp[:, 0].reshape(Nx, Ny);
    #Y = tmp[:, 1].reshape(Nx, Ny);
    phi1_0 = tmp[:, 2].reshape(Nx, Ny);
    phi2_0 = np.load(phi2_file)[:,2].reshape(Nx,Ny);
    #phi2_0 = np.loadtxt(phi2_file)[:, 2].reshape(Nx, Ny);

  # -- paramters --
  f = 1./50. # fraction of inner region area/total area
  sig = 0.e-3 # noise strength
  # -- for piecewise gradient --
  dy = 0.50
  dphi = 0.15
  #dphi = -np.linspace(0.05, 0.2, 4)[runid];
  if (runid < 5):
    m = dphi/dy
  else:
    m = -dphi/dy
  # -- for step down --
  #dphi = 0.18 #np.linspace(0.05, 0.2, 4)[runid]; # a step of greater than 0.1 should cross the spinodal
  #dy1 = 0.02
  #dy2 = 0.08
  #dy3 = 0.20
  #m1 = -dphi/dy1
  #m3 = dphi/dy3

  y_half = np.linspace(0, 1, Ny/2)

  # simple linear gradient
  phi2_grad = np.zeros((Nx, Ny));
  phi2_grad[:, 0:Ny/2] = -m*(Y[:, 0:Ny/2]/Ly - 0.5 + f/2.)
  phi2_grad[:, Ny/2:] = m*(Y[:,Ny/2:]/Ly- 0.5 - f/2.)

  # outer domain (piecewise gradient)
  #phi_half = np.piecewise(y_half, \
  #  [y_half<f, (y_half>=f) & (y_half < (f+dy)), y_half >= (f+dy)], \
  #  [0,        lambda y_half: m*(y_half-f),     m*dy])

  # inner domain (piecewise gradient)
  #phi_half = np.piecewise(y_half, \
  #  [y_half<f-dy, ((y_half >= (f-dy)) & (y_half<f)), y_half >= f], 
  #  [m*dy, lambda y_half: -m*(y_half-f), 0])

  # inner and outer domain (piecewise gradient)
  #phi_half = np.piecewise(y_half, \
  #  [#
  #    y_half<f-dy,#
  #    (y_half >= (f-dy)) & (y_half<f),#
  #    y_half >= f,#
  #    y_half<f,#
  #    (y_half>=f) & (y_half < (f+dy)),#
  #    y_half >= (f+dy)#
  #  ],[#
  #    m*dy,#
  #    lambda y_half: -m*(y_half-f),#
  #    0,#
  #    0,#
  #    lambda y_half: m*(y_half-f),#
  #    m*dy#
  #  ])

  # inner domain (piecewise step down)
  #phi_half = np.piecewise(y_half, #
  #  [#  
  #    y_half<(f-dy1-dy2-dy3), #
  #    ((y_half >= (f-dy1-dy2-dy3)) & (y_half<(f-dy1-dy2))), #
  #    ((y_half >= (f-dy1-dy2)) & (y_half<(f-dy1))), #
  #    ((y_half >= (f-dy1)) & (y_half<f)), #
  #    y_half >= f # 
  #  ], # 
  #  [#  
  #    0,# 
  #    lambda y_half: m3*(y_half-f+dy1+dy2+dy3),#
  #    dphi,#
  #    lambda y_half: m1*(y_half-f),#
  #     0#  
  #  ])
  #phi_whole = np.concatenate((phi_half[::-1], phi_half))
  #phi2_grad = np.tile(phi_whole, (Nx, 1))

  # inner domain, new interfaces
  #phi1_Leq = phi1_0[0, 0]
  #phi1_Req = phi1_0[0, -1]
  #phi2_Leq = phi2_0[0, 0]
  #phi2_Req = phi2_0[0, -1]
  #y_L = abs(y_half-0.11).argmin()
  #y_R = abs(y_half-0.09).argmin()
  #phi1_L =  np.concatenate( ( phi1_0[:, y_L::]-phi1_Leq,  np.ones((X.shape[0], y_L))*(phi1_Req-phi1_Leq)), axis=1)
  #phi1_R =  np.concatenate( (-phi1_0[:, y_R::]+phi1_Leq, -np.ones((X.shape[0], y_R))*(phi1_Req-phi1_Leq)), axis=1)
  #phi1_int = phi1_L + phi1_R
  #phi2_L =  np.concatenate( ( phi2_0[:, y_L::]-phi2_Leq,  np.ones((X.shape[0], y_L))*(phi2_Req-phi2_Leq)), axis=1)
  #phi2_R =  np.concatenate( (-phi2_0[:, y_R::]+phi2_Leq, -np.ones((X.shape[0], y_R))*(phi2_Req-phi2_Leq)), axis=1)
  #phi2_int = phi2_L + phi2_R

  phi1_0 = phi1_0 #+ phi1_int
  phi2_0 = phi2_0 + phi2_grad #+ phi2_int

  # random perturbation
  phi1_rand = sig*(np.random.rand(*Y.shape)-0.5)
  phi2_rand = sig*(np.random.rand(*Y.shape)-0.5)

  # mask for high k-modes
  dx = float(Lx)/Nx;
  dy = float(Ly)/Ny;
  kx = 2*np.pi*np.fft.fftfreq(X.shape[0])/dx
  ky = 2*np.pi*np.fft.fftfreq(X.shape[1])/dy
  Nyquist_kx = np.pi/dx; # maximum wavevector in x
  Nyquist_ky = np.pi/dy; # maximum wavevector in y
  x_mask = np.less_equal(abs(kx), Nyquist_kx/2.);
  y_mask = np.less_equal(abs(ky), Nyquist_ky/2.);
  mask = np.outer(x_mask, y_mask);

  # apply low-pass filter to noise
  phi1_rand_k = np.fft.fftn(phi1_rand)
  phi2_rand_k = np.fft.fftn(phi2_rand)
  phi1_rand_k = phi1_rand_k*mask;
  phi2_rand_k = phi2_rand_k*mask;
  phi1_rand = np.fft.ifftn(phi1_rand_k).real
  phi2_rand = np.fft.ifftn(phi2_rand_k).real

  # ensure phi1 and phi2 are symmetric in y after adding noise
  phi1_rand = 0.5*(phi1_rand + phi1_rand[:, ::-1])
  phi2_rand = 0.5*(phi2_rand + phi2_rand[:, ::-1])

  # apply a velocity perturbation to advect phi one step
  A=0.5; # amplitude of velocity perturbation
  lam_x = 2.*Lx/(np.array([1., 2., 4., 8., 16.])[runid%5]);
  lam_y = Ly/4.;
  xc = 0. #Lx/2.;
  yc = 0. #Ly/2.;

  #U =  1j*Ky[n,m]*np.exp(1j*Kx[n,m]*X+1j*Ky[n,m]*Y)
  #V = -1j*Kx[n,m]*np.exp(1j*Kx[n,m]*X+1j*Ky[n,m]*Y)

  # in x and y
  vx = A*np.cos(2*np.pi*(X-xc)/lam_x)*np.cos(2*np.pi*(Y-yc)/lam_y)
  vy = A*np.sin(2*np.pi*(X-xc)/lam_x)*np.sin(2*np.pi*(Y-yc)/lam_y)
  vz = A*np.zeros(X.shape)

  # in x only
  #vx = A*np.cos(2*np.pi*(X-xc)/lam_x)
  #vy = A*np.zeros(Y.shape)
  #vz = A*np.zeros(X.shape)

  # now, advect phi_1 and phi_2 one step
  (Kx, Ky) = np.meshgrid(kx, ky, indexing='ij')
  phi1_k = np.fft.fftn(phi1_0)
  phi2_k = np.fft.fftn(phi2_0)
  phi1_adv = - vx*np.fft.ifftn(1j*Kx*phi1_k).real - vy*np.fft.ifftn(1j*Ky*phi1_k).real
  phi2_adv = - vx*np.fft.ifftn(1j*Kx*phi2_k).real - vy*np.fft.ifftn(1j*Ky*phi2_k).real

  phi[0] = phi1_0 + phi1_adv + phi1_rand
  phi[1] = phi2_0 + phi2_adv + phi2_rand

  # i/o
  ic_dict = {'IC_flag':IC_flag, \
             'phi1_file':phi1_file, \
             'phi2_file':phi2_file, \
             'f':f, \
             'sig':sig}
  #f1.write('%.16e\t'%m + '# m (slope for gradient)\n')
# }}}
elif (IC_flag == 9): # 1D sinusoid in x or y
# {{{

  #phi_array = np.linspace(0.49, 0.42, 5); #np.array([0.05, 0.15, 0.25, 0.35, 0.45])
  #phi1_array = np.loadtxt('params_PhiAv.in')[:, 0];
  #phi2_array = np.loadtxt('params_PhiAv.in')[:, 1];

  #params
  # note: phi1_0 and phi2_0 cannot be *exactly* equal, add 1e-3 or something
  #phi1_0 = phi1_array[1] # avg concentration of polymer 
  #phi2_0 = phi2_array[1] # avg concentration of non-solvent
  #phi1_0 = phi1_array[runid] # avg concentration of polymer
  #phi2_0 = phi2_array[runid] # avg concentration of non-solvent
  phi1_0 = 0.25
  phi2_0 = 0.15
  T = 2 # periods for sinusoid
  A = 0.1 # sinusoid amplitude

  # phi
  phi[0] = A*np.sin(2*T*np.pi*X/Lx) + phi1_0;
  phi[1] = A*np.sin(2*T*np.pi*X/Lx) + phi2_0;
  #phi1 = A*np.sin(2*T*np.pi*Y/Ly) + phi1_0;
  #phi2 = A*np.sin(2*T*np.pi*Y/Ly) + phi2_0;

  # i/o
  ic_dict = {'IC_flag':IC_flag, \
             'phi1_0':phi1_0, \
             'phi2_0':phi2_0, \
             'T':T, \
             'A':A}
# }}}
elif (IC_flag == 10): # 2D Gaussian pulse
# {{{
  #params
  phi1_0 = 0.25 # avg concentration of polymer
  phi2_0 = 0.1 # avg concentration of non-solvent
  sig2_1 = 1. # variance for phi1
  sig2_2 = 1. # variance for phi2
  A = 0.2 # gaussian amplitude

  # phi
  # normalized (unbounded)
  #phi[0] = A/(2*np.pi*sig2_1)*np.exp((X**2+Y**2)/(2*sig2_1)) + phi1_0;
  #phi[1] = A/(2*np.pi*sig2_2)*np.exp((X**2+Y**2)/(2*sig2_2)) + phi2_0;

  # not-normalized ("A" dictates scale)
  #phi[0] = -A*np.exp(-((X-3*Lx/4)**2+(Y-Ly/2)**2)/(2*sig2_1)) + phi1_0;
  #phi[1] = A*np.exp(-((X-1*Lx/4)**2+(Y-Ly/2)**2)/(2*sig2_2)) + phi2_0;

  # 1D gaussian (X)
  phi[0] = -A*np.exp(-((X-3*Lx/4)**2)/(2*sig2_1)) + phi1_0;
  phi[1] = A*np.exp(-((X-1*Lx/4)**2)/(2*sig2_2)) + phi2_0;

  # 1D gaussian (Y)
  #phi[0] = -A*np.exp(-((Y-1*Ly/4)**2)/(2*sig2_1)) + phi1_0;
  #phi[1] = A*np.exp(-((Y-3*Ly/4)**2)/(2*sig2_2)) + phi2_0;

  # i/o
  ic_dict = {'IC_flag':IC_flag, \
             'phi1_0':phi1_0, \
             'phi2_0':phi2_0, \
             'sig2_1':sig2_1, \
             'sig2_2':sig2_2, \
             'A':A}
# }}}
elif (IC_flag == 11): # 3D sinusoid
# {{{
  #params
  phi1_0 = 0.25 # avg concentration of polymer
  phi2_0 = 0.60 # avg concentration of non-solvent
  Tx = 1 # periods for sinusoid in x-dir
  Ty = 3 # periods for sinusoid in y-dir
  Tz = 5 # periods for sinusoid in y-dir
  A = 0.05 # sinusoid amplitude

  # phi
  phi[0] = A*np.sin(2*Tx*np.pi*X/Lx)*np.sin(2*Ty*np.pi*Y/Ly)*np.sin(2*Tz*np.pi*Z/Lz) + phi1_0;
  phi[1] = A*np.sin(2*Tx*np.pi*X/Lx)*np.sin(2*Ty*np.pi*Y/Ly)*np.sin(2*Tz*np.pi*Z/Lz) + phi2_0;

  # i/o
  ic_dict = {'IC_flag':IC_flag, \
             'phi1_0':phi1_0, \
             'phi2_0':phi2_0, \
             'Tx':Tx, \
             'Ty':Ty, \
             'Tz':Tz, \
             'A':A}
# }}}
elif (IC_flag == 12): # 2D circle/ellipse
# {{{

  #eq_phi = np.loadtxt("phase_diagram_bin.in");
  #alpha = eq_phi[:, :3]
  #beta  = eq_phi[::-1, 3:]

  #params
  phi1_in = 0.98 #0.60 #alpha[79, 0]#0.01 #0.4 # avg concentration of polymer
  phi2_in = 0.01 #0.31 #alpha[79, 1]#0.98 #0.5 # avg concentration of non-solvent
  phi1_out = 0.001 #0.03 #beta[79, 0]#0.25 #0.01 # avg concentration of polymer
  phi2_out = 0.499 #0.86 #beta[79, 1]#0.6 #0.98 # avg concentration of non-solvent

  sig2 = 1.e-3 # noise variance
  w = 0.5 # thickness of tanh boundary
  Rx = 0.05*Lx # ellipse x-axis
  Ry = 0.05*Ly # ellipse y-axis 
  my_seed = np.random.randint(0, 4294967295) #4294... is the maximum seed value

  # basic phi step parallel to the x-axis
  X_c = (X - Lx/2)/Rx
  Y_c = (Y - Ly/2)/Ry
  dphi1 = phi1_out - phi1_in
  dphi2 = phi2_out - phi2_in
  phi1_0 = 0.5*dphi1*(np.tanh((X_c**2 + Y_c**2 - 1)/w) + 1) + phi1_in
  phi2_0 = 0.5*dphi2*(np.tanh((X_c**2 + Y_c**2 - 1)/w) + 1) + phi2_in
  # perhaps make this elliptical? X^2/Rx^2 + Y^2/Ry^2 - 1 == 0

  # get random perturbation
  np.random.seed(my_seed)
  phi1_rand = sig2*(np.random.rand(*X.shape)-0.5)
  phi2_rand = sig2*(np.random.rand(*X.shape)-0.5)
  phi[0] = phi1_0+phi1_rand
  phi[1] = phi2_0+phi2_rand

  # i/o
  ic_dict = {'IC_flag':IC_flag, \
             'phi1_in':phi1_in, \
             'phi2_in':phi2_in, \
             'phi1_out':phi1_out, \
             'phi2_out':phi2_out, \
             'sig2':sig2, \
             'w':w, \
             'Rx':Rx, \
             'Ry':Ry, \
             'my_seed':my_seed}
# }}}
elif (IC_flag == 13): # Multiple spheres at random locations
# {{{
  #params
  phi1_in = 0.98     # avg concentration of polymer
  phi2_in = 0.01    # avg concentration of non-solvent
  phi1_out = 0.001   # avg concentration of polymer
  phi2_out = 0.499  # avg concentration of non-solvent
  sig2 = 1.e-3      # noise variance
  w = 0.3           # thickness of tanh boundary
  loading = 0.2
  volbox = Lx*Ly*Lz

  Rx = 1
  Ry = 1
  Rz = 1
  volsphere = Rx*Ry*Rz*4/3*np.pi
  my_seed  = np.random.randint(0, 2**32-1)
  np.random.seed(my_seed) # set random seed

  num_spheres = int(loading*volbox/volsphere) # number of spheres

  # basic phi step parallel to the x-axis
  dphi1  = phi1_out - phi1_in
  dphi2  = phi2_out - phi2_in

  xs = np.zeros(num_spheres)  #array of x-axis locations of spheres
  ys = np.zeros(num_spheres)  #array of y-axis locations of spheres
  zs = np.zeros(num_spheres)  #array of z-axis locations of spheres
  def checkOverlap(i,newx,newy,newz):
    
    for j in range(i):
      dx = abs(newx - xs[j])
      dy = abs(newy - ys[j])
      dz = abs(newz - zs[j])
      if dx > 0.5*Lx:
        dx = Lx-dx
      if dy > 0.5*Ly:
        dy = Ly-dy
      if dz > 0.5*Lz:
        dz = Lz-dz

      dist = np.sqrt( (dx)**2 + (dy)**2 + (dz)**2 )
      if(dist < 2.1*Rx):
        return(1)

    return(0)

  res = Nx/Lx
  
  for i in range(num_spheres):
    overlap = 1
    while(overlap):
      my_seed  = np.random.randint(0, 2**32-1)
      np.random.seed(my_seed)

      newx = np.random.randint(0, Nx+1)/res
      newy = np.random.randint(Ny/4, 3*Ny/4+1)/res
      newz = np.random.randint(0, Nz+1)/res

      overlap = checkOverlap(i,newx,newy,newz)
      xs[i] = newx
      ys[i] = newy
      zs[i] = newz  
  tanh_sum = 0
  def make_sphere(i):
    py = np.where((Y-ys[i]) <= -2*Ry, (Y-ys[i]+Ly), (Y-ys[i]))
    py = np.where((py) > 2*Ry, (py-Ly), py)
    px = np.where((X-xs[i]) <= -2*Rx, (X-xs[i]+Lx), (X-xs[i]))
    px = np.where((px) > 2*Rx, (px-Lx), px)
    pz = np.where((Z-zs[i]) <= -2*Rz, (Z-zs[i]+Lz), (Z-zs[i]))
    pz = np.where((pz) > 2*Rz, (pz-Lz), pz)
    return (np.tanh(((py)**2/Ry**2 + (px)**2/Rx**2 + (pz)**2/Rz**2 - 1)/w) + 1)
  for i in range(num_spheres):
    tanh_sum = tanh_sum + make_sphere(i)

  lb = 0.2; # lower bound
  ub = 0.8; # upper bound
  Y_lb = (Y-lb*Ly)/w
  Y_ub = (ub*Ly-Y)/w
  phi1_0 = 0.5*dphi1*tanh_sum - (num_spheres-1)*dphi1
  phi2_0 = (0.5*dphi2*tanh_sum - (num_spheres-1)*dphi2)*0.5*(np.tanh(Y_ub)+np.tanh(Y_lb))

  #phi1_0 = 0.5*dphi1*np.tanh(Y_lb) + 0.5*dphi1*np.tanh(Y_ub) + phi1_out
  #phi2_0 -= 0.5*dphi2*np.tanh(Y_ub)# + 0.5*dphi2*np.tanh(Y_lb)
  lb = 0.8; # lower bound
  ub = 1; # upper bound
  Y_lb = (Y-lb*Ly)/w
  Y_ub = (ub*Ly-Y)/w
  #phi2_0 -= 0.5*dphi2*np.tanh(Y_ub) + 0.5*dphi2*np.tanh(Y_lb)

  # get random perturbation
  phi1_rand = sig2*(np.random.rand(*X.shape)-0.5)
  phi2_rand = sig2*(np.random.rand(*X.shape)-0.5)
  
  phi[0] = phi1_0 + phi1_rand + phi1_in
  phi[1] = phi2_0 + phi2_rand + phi2_in

  # i/o
  ic_dict = {'IC_flag':IC_flag, \
             'phi1_in':phi1_in, \
             'phi2_in':phi2_in, \
             'phi1_out':phi1_out, \
             'phi2_out':phi2_out, \
             'sig2':sig2, \
             'w':w, \
             'num_spheres':num_spheres, \
             'Rx':Rx, \
             'Ry':Ry, \
             'Rz':Rz, \
             'my_seed':my_seed}
# }}}
elif (IC_flag == 14): # Retracting droplet
# {{{  
  #params
  phi1_in  = 0.02 #0.60 #alpha[79, 0]#0.01 #0.4 # avg concentration of polymer
  phi2_in  = 0.952 #0.31 #alpha[79, 1]#0.98 #0.5 # avg concentration of non-solvent
  phi1_ins = 0.85
  phi2_ins = 0.05
  phi1_out = 0.02 #0.03 #beta[79, 0]#0.25 #0.01 # avg concentration of polymer
  phi2_out = 0.048 #0.86 #beta[79, 1]#0.6 #0.98 # avg concentration of non-solvent
  sig2     = 1.e-3  # noise variance
  w        = 0.3    # thickness of tanh boundary
  Rx       = 10     # ellipse x-axis
  Ry       = 30     # ellipse y-axis
  Rxs      = 2      # droplet x-axis
  Rys      = 2      # droplet y-axis
  my_seed  = np.random.randint(0, 4294967295) #4294... is the maximum seed value

  # basic phi step parallel to the x-axis
  dphi1  = phi1_out - phi1_in
  dphi2  = phi2_out - phi2_in
  dphi1s = phi1_in  - phi1_ins
  dphi2s = phi2_in  - phi2_ins # this is arbitrary (might need to change for different conditions or radius values or positions)

  xs = np.zeros(2)
  ys = np.zeros(2)
  xs[0] = Lx/2-1
  ys[0] = Ly/2-1
  xs[1] = Lx/2-1
  ys[1] = Ly/2-1-Ry/2
  
  tanh_sum  = 0
  tanh_sum2 = -2
  tanh_sum  = tanh_sum  + (np.tanh(((Y - ys[0])**2/Ry**2  + (X - xs[0])**2/Rx**2  - 1)/(w)) + 1)
  tanh_sum2 = tanh_sum2 + (np.tanh(((Y - ys[1])**2/Rys**2 + (X - xs[1])**2/Rxs**2 - 1)/(w)) + 1)

  phi1_0 = 0.5*dphi1*(tanh_sum )
  phi2_0 = 0.5*dphi2*(tanh_sum )
  phi1_s = 0.5*dphi1s*(tanh_sum2 )
  phi2_s = 0.5*dphi2s*(tanh_sum2 )
  
  # get random perturbation
  np.random.seed(my_seed)
  phi1_rand = sig2*(np.random.rand(*X.shape)-0.5)
  phi2_rand = sig2*(np.random.rand(*X.shape)-0.5)
  phi[0] = phi1_0 + phi1_rand + phi1_in + phi1_s
  phi[1] = phi2_0 + phi2_rand + phi2_in + phi2_s

  # i/o
  ic_dict = {'IC_flag':IC_flag, \
             'phi1_in':phi1_in, \
             'phi2_in':phi2_in, \
             'phi1_out':phi1_out, \
             'phi2_out':phi2_out, \
             'sig2':sig2, \
             'w':w, \
             'Rx':Rx, \
             'Ry':Ry, \
             'my_seed':my_seed}

# }}}
elif (IC_flag == 15): # Interfacial droplet
# {{{  
  #params
  my_seed  = np.random.randint(0, 4294967295) #4294... is the maximum seed value
  phi1_in = 0.9      #0.60 #alpha[79, 0]#0.01 #0.4 # avg concentration of polymer
  phi2_in = 0.05      #0.31 #alpha[79, 1]#0.98 #0.5 # avg concentration of non-solvent
  phi1_out = 0.05     #0.03 #beta[79, 0]#0.25 #0.01 # avg concentration of polymer
  phi2_out1 = 0.9
  phi2_out2 = 0.1
  dphi1  = phi1_out - phi1_in
  dphi21 = phi2_out1 - phi2_in
  dphi22 = phi2_out2 - phi2_in
  sig2 = 1.e-3        # noise variance
  w = 0.3             # thickness of tanh boundary
  
  Rx = 5          # droplet x-axis
  Ry = 5          # droplet y-axis 
  
  X_c = np.zeros([Nx,Nx])
  Y_c = np.zeros([Ny,Ny])
  X_c = (X - Lx/2)/Rx
  Y_c = (Y - Ly/2)/Ry
  
  phi1_0  = 0.5*dphi1  * (np.tanh((X_c**2 + Y_c**2 - 1)/w) + 1) 
  phi2_01 = 0.5*dphi21 * (np.tanh((X_c**2 + Y_c**2 - 1)/w) + 1)
  phi2_02 = 0.5*dphi22 * (np.tanh((X_c**2 + Y_c**2 - 1)/w) + 1)
  phi2_03 = np.zeros([Ny,Ny])
  phi2_03[:,:int(Nx/2)] = phi2_01[:,:int(Nx/2)]
  phi2_03[:,int(Nx/2):] = phi2_02[:,int(Nx/2):]
  
  # get random perturbation
  np.random.seed(my_seed)
  phi1_rand = sig2*(np.random.rand(*X.shape)-0.5)
  phi2_rand = sig2*(np.random.rand(*X.shape)-0.5)
  phi[0] = phi1_0  + phi1_rand + phi1_in
  phi[1] = phi2_03 + phi2_rand + phi2_in

  # i/o
  ic_dict = {'IC_flag':IC_flag, \
             'phi1_in':phi1_in, \
             'phi2_in':phi2_in, \
             'phi1_out':phi1_out, \
             'sig2':sig2, \
             'w':w, \
             'Rx':Rx, \
             'Ry':Ry, \
             'my_seed':my_seed}
# }}}
elif (IC_flag == 16): # VIPS
# {{{
  ### We introduce a diffuse interface between the left domain (film)
  ### and the right domain (bath/atmosphere) where the concentration
  ### is set to a constant. This IC is useful for DIRICHLET boundary
  ### conditions for VIPS. To use Neumann with VIPS, use IC_flag = 0
  
  # params
  bc_params = np.loadtxt('params_BCs_phi.in')
  phi2_ord = bc_params[9]
  phi2_BC = bc_params[17]
 
  if (phi2_ord == 0): #check it's Dirichlet
    #film compositions
    phi1_0 = 0.15
    phi2_0 = 0.20 
    w = 0.8 # boundary thickness
    f = 0.85 #5./12. # fraction of inner region area/total area

    # basic phi step parallel to the x-axis
    dphi2 = phi2_BC - phi2_0
    phi2_1 = 0.5*dphi2*np.tanh(w*(X-f*Lx)) + (phi2_0 + phi2_BC)/2

    # phi
    phi[0] = np.zeros(X.shape) + phi1_0
    phi[1] = phi2_1

    # i/o
    ic_dict = {'IC_flag':IC_flag, \
               'phi1_0':phi1_0, \
               'phi2_0':phi2_0, \
               'phi2_BC':phi2_BC}

  else:
    print("Error: this BC is for Dirichlet VIPS. For Neumann VIPS, use option 0")
    sys.exit();
# }}}
elif (IC_flag == 17): # Multiple spheres at random locations 2D
#{{{
  #params
  phi1_in = 0.98     # avg concentration of polymer
  phi2_in = 0.01    # avg concentration of non-solvent
  phi1_out = 0.001   # avg concentration of polymer
  phi2_out = 0.499  # avg concentration of non-solvent
  sig2 = 1.e-3      # noise variance
  w = 0.5           # thickness of tanh boundary
  loading = 0.2
  volbox = Lx*Ly

  Rx = 1
  Ry = 1

  volsphere = Rx*Ry*np.pi #change to circle
  my_seed  = np.random.randint(0, 2**32-1)
  np.random.seed(my_seed) # set random seed

  num_spheres = int(loading*volbox/volsphere) # number of spheres

  # basic phi step parallel to the x-axis
  dphi1  = phi1_out - phi1_in
  dphi2  = phi2_out - phi2_in

  xs = np.zeros(num_spheres)  #array of x-axis locations of spheres
  ys = np.zeros(num_spheres)  #array of y-axis locations of spheres
  def checkOverlap(i,newx,newy):
    
    for j in range(i):
      dx = abs(newx - xs[j])
      dy = abs(newy - ys[j])
      if dx > 0.5*Lx:
        dx = Lx-dx
      if dy > 0.5*Ly:
        dy = Ly-dy
      #adjusting the dist < makes spheres generate futher apart from eachother (2.0 means they are touching)
      dist = np.sqrt( (dx)**2 + (dy)**2 )
      if(dist < 2.5*Rx):
        return(1)

    return(0)

  res = Nx/Lx
  
  for i in range(num_spheres):
    overlap = 1
    while(overlap):
      my_seed  = np.random.randint(0, 2**32-1)
      np.random.seed(my_seed)

      #this is where the spheres can be placed (zero is the bottom of the box 1 is the top)
      newx = np.random.randint(0, Nx+1)/res
      newy = np.random.randint(0, Ny+1)/res

      overlap = checkOverlap(i,newx,newy)
      xs[i] = newx
      ys[i] = newy
        
  tanh_sum = 0
  def make_sphere(i):
    py = np.where((Y-ys[i]) <= -2*Ry, (Y-ys[i]+Ly), (Y-ys[i]))
    py = np.where((py) > 2*Ry, (py-Ly), py)
    px = np.where((X-xs[i]) <= -2*Rx, (X-xs[i]+Lx), (X-xs[i]))
    px = np.where((px) > 2*Rx, (px-Lx), px)
    return (np.tanh(((py)**2/Ry**2 + (px)**2/Rx**2 - 1)/w) + 1)
  for i in range(num_spheres):
    tanh_sum = tanh_sum + make_sphere(i)

  lb = 0.0; # lower bound
  ub = 1.0; # upper bound
  Y_lb = (Y-lb*Ly)/w
  Y_ub = (ub*Ly-Y)/w
  phi1_0 = 0.5*dphi1*tanh_sum - (num_spheres-1)*dphi1
  phi2_0 = (0.5*dphi2*tanh_sum - (num_spheres-1)*dphi2) * 0.5*(np.tanh(Y_ub)+np.tanh(Y_lb))

  #phi1_0 = 0.5*dphi1*np.tanh(Y_lb) + 0.5*dphi1*np.tanh(Y_ub) + phi1_out
  #phi2_0 -= 0.5*dphi2*np.tanh(Y_ub)# + 0.5*dphi2*np.tanh(Y_lb)
  lb = 0.8; # lower bound
  ub = 1; # upper bound
  Y_lb = (Y-lb*Ly)/w
  Y_ub = (ub*Ly-Y)/w
  #phi2_0 -= 0.5*dphi2*np.tanh(Y_ub) + 0.5*dphi2*np.tanh(Y_lb)

  # get random perturbation
  phi1_rand = sig2*(np.random.rand(*X.shape)-0.5)
  phi2_rand = sig2*(np.random.rand(*X.shape)-0.5)
  
  phi[0] = phi1_0 + phi1_rand + phi1_in
  phi[1] = phi2_0 + phi2_rand + phi2_in

  # i/o
  ic_dict = {'IC_flag':IC_flag, \
             'phi1_in':phi1_in, \
             'phi2_in':phi2_in, \
             'phi1_out':phi1_out, \
             'phi2_out':phi2_out, \
             'sig2':sig2, \
             'w':w, \
             'num_spheres':num_spheres, \
             'Rx':Rx, \
             'Ry':Ry, \
             'my_seed':my_seed}
# }}}
# --- (100 Series) Quaternary NIPS ---
elif (IC_flag == 100): # homogeneous initial condition
# {{{
  # params
  phi_0 = np.zeros(NIcomp);
  phi_0[0] = 0.25 # avg conc. for polymer
  phi_0[1] = 0.6 # avg conc. for nonsolvent
  phi_0[2] = 0.1 # avg conc. for solvent

  # phi
  for i in range(NIcomp):
    phi[i] = np.zeros(X.shape) + phi_0[i]
  phi_0list = phi_0.tolist()
  # i/o
  ic_dict = {'IC_flag':IC_flag, 
             'phi_0':phi_0list}
# }}}
elif (IC_flag == 101): # uniform noise about homogeneous initial conditions
# {{{
  # params
  phi_0 = np.zeros(NIcomp);
  phi_0[0] = 0.25 #0.25
  phi_0[1] = 0.4  #0.6
  phi_0[2] = 0.2  #0.

  sig = 1e-2 # std. deviation for random noise
  my_seed = 1234 #np.random.randint(0, 4294967295) #1234
  np.random.seed(my_seed)

  # phi
  for i in range(NIcomp):
    phi[i] = sig*np.random.rand(*X.shape)
    phi_k = np.fft.fftn(phi[i])
    phi_k[np.unravel_index(0, X.shape)] = phi_0[i]*Ngrid # set 0 element to mean
    phi[i] = np.fft.ifftn(phi_k).real

  phi_0list = phi_0.tolist()
  # i/o
  ic_dict = {'IC_flag':IC_flag, 
             'phi_0':phi_0list,
             'sig':sig,
             'my_seed':my_seed}
# }}}
elif (IC_flag == 102): # guassian random initial conditions
# {{{ 
  # params
  phi_0 = np.zeros(NIcomp);
  phi_0[0] = 0.25 # avg conc. for polymer
  phi_0[1] = 0.6 # avg conc. for nonsolvent
  phi_0[2] = 0.1 # avg conc. for solvent
  
  sig = 1e-2 # std. deviation for random noise (IC_flag == 1, 2, 4 only)
  my_seed = 1234 #np.random.randint(0, sys.maxint) #1234

  # mask for high k-modes
  dx = Lx/Nx;
  dy = Ly/Ny;
  kx = 2*np.pi*np.fft.fftfreq(X.shape[0])/dx
  ky = 2*np.pi*np.fft.fftfreq(X.shape[1])/dy
  Nyquist_kx = np.pi/dx; # maximum wavevector in x
  Nyquist_ky = np.pi/dy; # maximum wavevector in y
  x_mask = np.less_equal(abs(kx), Nyquist_kx/2.);
  y_mask = np.less_equal(abs(ky), Nyquist_ky/2.);
  mask = np.outer(x_mask, y_mask);

  # phi
  np.random.seed(my_seed)
  for i in range(NIcomp):
    phi[i] = sig*np.random.randn(*X.shape)
    phi_k = np.fft.fftn(phi[i])
    phi_k[np.unravel_index(0, X.shape)] = phi_0[i]*Ngrid
    phi_k = phi_k*mask;
    phi[i] = np.fft.ifftn(phi_k).real

  phi_0list = phi_0.tolist()
  # i/o
  ic_dict = {'IC_flag':IC_flag, 
             'phi_0':phi_0list,
             'sig':sig,
             'my_seed':my_seed}
# }}}
elif (IC_flag == 103): # tanh profile in y
# {{{

  #eq_phi = np.loadtxt("phase_diagram_bin.in");
  #alpha = eq_phi[:, :3]
  #beta  = eq_phi[::-1, 3:]

  #params
  #phi1_in = alpha[60, 0]#0.01 #0.4 # avg concentration of polymer
  #phi2_in = alpha[60, 1]#0.98 #0.5 # avg concentration of non-solvent
  #phi1_out = beta[60, 0]#0.25 #0.01 # avg concentration of polymer
  #phi2_out = beta[60, 1]#0.6 #0.98 # avg concentration of non-solvent

  #params
  phi_in = np.zeros(NIcomp);
  phi_out = np.zeros(NIcomp);
  phi_in[0] = 0.28 # polymer volume fraction in inner region
  phi_in[1] = 0.50 # non-solvent volume fraction in inner region
  phi_in[2] = 0.10 # solvent volume fraction in inner region
  phi_out[0] = 0.02 # polymer volume fraction in outer region
  phi_out[1] = 0.80 # non-solvent volume fraction in outer region
  phi_out[2] = 0.10 # solvent volume fraction in outer region
  sig = 1e-3; #9.9e-3 #0.99e-2 # noise strength
  w = 3. # boundary thickness
  f = 1./2. #5./12. # fraction of inner region area/total area
  my_seed = np.random.randint(0, 4294967295) #4294... is the maximum seed value

  # basic phi step parallel to the x-axis
  lb = 0.5 - 0.5*f; # lower bound
  ub = 0.5 + 0.5*f; # upper bound
  Y_lb = (Y-lb*Ly)/w
  Y_ub = (ub*Ly-Y)/w
  dphi = np.zeros(NIcomp);
  for i in range(NIcomp):
    dphi[i] = phi_in[i] - phi_out[i]
  np.random.seed(my_seed)
  phi_0 = np.zeros(NIcomp, dtype=object);
  for i in range(NIcomp):
    phi_0[i] = 0.5*dphi[i]*np.tanh(Y_lb) + 0.5*dphi[i]*np.tanh(Y_ub) + phi_out[i] 

  # get random perturbation
  phi_rand = np.zeros(NIcomp, dtype=object);
  for i in range(NIcomp):
    phi_rand[i] = sig*(np.random.rand(*Y.shape)-0.5)

  # mask for high k-modes
  dx = float(Lx)/Nx;
  dy = float(Ly)/Ny;
  kx = 2*np.pi*np.fft.fftfreq(X.shape[0])/dx
  ky = 2*np.pi*np.fft.fftfreq(X.shape[1])/dy
  Nyquist_kx = np.pi/dx; # maximum wavevector in x
  Nyquist_ky = np.pi/dy; # maximum wavevector in y
  x_mask = np.less_equal(abs(kx), Nyquist_kx/2.);
  y_mask = np.less_equal(abs(ky), Nyquist_ky/2.);
  mask = np.outer(x_mask, y_mask);

  # apply low-pass filter to noise
  for i in range(NIcomp):
    phi_rand_k = np.fft.fftn(phi_rand[i])
    phi_rand_k = phi_rand_k*mask;
    phi_rand = np.fft.ifftn(phi_rand_k).real
    # ensure phi[i] are symmetric in y after adding noise
    phi_rand = 0.5*(phi_rand + phi_rand[:, ::-1])
    phi[i] = phi_0[i] #+phi_rand

  phi_inlist = phi_in.tolist()
  phi_outlist = phi_out.tolist()
  # i/o
  ic_dict = {'IC_flag':IC_flag, 
             'phi_in':phi_inlist,
             'phi_out':phi_outlist,
             'sig':sig,
             'w':w,
             'f':f,
             'my_seed':my_seed}
# }}}
elif (IC_flag == 104): # sinusoids in x and y
# {{{
  #params
  phi_0 = np.zeros(NIcomp);
  phi_0[0] = 0.25 # avg concentration of polymer
  phi_0[1] = 0.60 # avg concentration of non-solvent
  phi_0[2] = 0.10 # avg concentration of solvent
  Tx = 1 # periods for sinusoid in x-dir
  Ty = 3 # periods for sinusoid in y-dir
  A = 0.05 # sinusoid amplitude

  # phi
  for i in range(NIcomp):
    phi[i] = A*np.sin(2*Tx*np.pi*X/Lx)*np.sin(2*Ty*np.pi*Y/Ly) + phi_0[i];

  phi_0list = phi_0.tolist()
  # i/o
  ic_dict = {'IC_flag':IC_flag, 
             'phi_0':phi_0list,
             'Tx':Tx,
             'Ty':Ty,
             'A':A}
# }}}
elif (IC_flag == 105): # gradient in y-direction
# {{{
  #params
  phi_N = np.zeros(NIcomp);
  phi_S = np.zeros(NIcomp);
  phi_N[0] = 0.30 # concentration of polymer on north end
  phi_S[0] = 0.30 # concentration of polymer on south end
  phi_N[1] = 0.65 # concentration of non-solvent on north end
  phi_S[1] = 0.50 # concentration of non-solvent on south end
  phi_N[2] = 0.65 # concentration of solvent on north end
  phi_S[2] = 0.50 # concentration of solvent on south end
  sig = 5e-3 # std. deviation for random noise
  my_seed = np.random.randint(0, 4294967295) #1234

  # phi
  np.random.seed(my_seed)
  for i in range(NIcomp):
    phi[i] = sig*np.random.rand(*X.shape)
    phi_k = np.fft.fftn(phi[i])
    phi_k[np.unravel_index(0, X.shape)] = 0 # set mean to 0
    phi[i] = np.fft.ifftn(phi_k).real
    phi[i][:, :int(Ny/2)] = phi[i][:, :int(Ny/2)] - 2*(phi_N[i]-phi_S[i])*Y[:, :int(Ny/2)]/Ly + phi_N[i];
    phi[i][:, int(Ny/2):] = phi[i][:, int(Ny/2):] + 2*(phi_N[i]-phi_S[i])*(Y[:, int(Ny/2):]/Ly-0.5) + phi_S[i];

  phi_Nlist = phi_N.tolist()
  phi_Slist = phi_S.tolist()
  # i/o
  ic_dict = {'IC_flag':IC_flag, 
             'phi_N':phi_Nlist,
             'phi_S':phi_Slist,
             'sig':sig,
             'my_seed':my_seed}
# }}}
elif (IC_flag == 106): # 1D tanh profile in x
# {{{

  #phi_params = np.loadtxt("../../setup/run"+str(runid)+"/params_Equilibrium.in")

  #params
  phi_in = np.zeros(NIcomp);
  phi_out = np.zeros(NIcomp);
  phi_in[0] = 0.54 #phi_params[0] # polymer volume fraction in inner region
  phi_in[1] = 0.35 #phi_params[1] # non-solvent volume fraction in inner region
  phi_in[2] = 0.25 # solvent volume fraction in inner region
  phi_out[0] = 0.01 #phi_params[2] # polymer volume fraction in outer region
  phi_out[1] = 0.85 #phi_params[3] # non-solvent volume fraction in outer region
  phi_out[2] = 0.85 # solvent volume fraction in outer region
  f = 0.5 #phi_params[4] #1./2. # fraction of inner region area/total area
  sig = 0.e-3 #0.5e-3 # noise strength
  w = 2 # boundary thickness
  my_seed = np.random.randint(0, 4294967295) #4294... is the maximum seed value

  # mask for high k-modes
  dx = float(Lx)/Nx;
  kx = 2*np.pi*np.fft.fftfreq(X.shape[0])/dx
  Nyquist_kx = np.pi/dx; # maximum wavevector in x
  x_mask = np.less_equal(abs(kx), Nyquist_kx/2.);

  # phi
  lb = 0.5 - 0.5*f; # lower bound
  ub = 0.5 + 0.5*f; # upper bound
  X_lb = (X-lb*Lx)/w
  X_ub = (ub*Lx-X)/w
  dphi = np.zeros(NIcomp);
  for i in range(NIcomp):
    dphi[i] = phi_in[i] - phi_out[i]
  np.random.seed(my_seed)
  for i in range(NIcomp):
    phi[i] = 0.5*dphi[i]*np.tanh(X_lb) + 0.5*dphi[i]*np.tanh(X_ub) + phi_out[i] + sig*(np.random.rand(*X.shape)-0.5)

  # apply low-pass filter
  for i in range(NIcomp):
    phi_k = np.fft.fft(phi[i])
    phi_k = phi_k*x_mask;
    phi[i] = np.fft.ifft(phi_k).real

  # ensure phi1 and phi2 are symmetric in x after adding noise
  for i in range(NIcomp):
    phi[i] = 0.5*(phi[i] + phi[i][::-1])

  phi_inlist = phi_in.tolist()
  phi_outlist = phi_out.tolist()
  # i/o
  ic_dict = {'IC_flag':IC_flag, 
             'phi_in':phi_inlist,
             'phi_out':phi_outlist,
             'sig':sig,
             'w':w,
             'f':f,
             'my_seed':my_seed}
# }}}
elif (IC_flag == 107): # tanh profile in y with roll cells
# {{{
  #params
  phi_in = np.zeros(NIcomp);
  phi_out = np.zeros(NIcomp);
  phi_in[0] = 0.35 # polymer volume fraction in inner region
  phi_in[1] = 0.45 # non-solvent volume fraction in inner region
  phi_in[2] = 0.45 # solvent volume fraction in inner region
  phi_out[0] = 0.02 # polymer volume fraction in outer region
  phi_out[1] = 0.97 # non-solvent volume fraction in outer region
  phi_out[2] = 0.97 # solvent volume fraction in outer region
  sig = 9.9e-3 #0.99e-2 # noise strength
  w = 4 # boundary thickness
  f = 1./2. #5./12. # fraction of inner region area/total area
  my_seed = np.random.randint(0, 4294967295) #4294... is the maximum seed value

  # basic phi step parallel to the x-axis
  lb = 0.5 - 0.5*f; # lower bound
  ub = 0.5 + 0.5*f; # upper bound
  Y_lb = (Y-lb*Ly)/w
  Y_ub = (ub*Ly-Y)/w
  dphi = np.zeros(NIcomp);
  for i in range(NIcomp):
    dphi[i] = phi_in[i] - phi_out[i]
  np.random.seed(my_seed)
  phi_0 = np.zeros(NIcomp, dtype=object); 
  for i in range(NIcomp):
    phi_0[i] = 0.5*dphi[i]*np.tanh(Y_lb) + 0.5*dphi[i]*np.tanh(Y_ub) + phi_out[i] 

  # get random perturbation
  phi_rand = np.zeros(NIcomp, dtype=object);
  for i in range(NIcomp):
    phi_rand[i] = sig*(np.random.rand(*Y.shape)-0.5)

  # mask for high k-modes
  dx = float(Lx)/Nx;
  dy = float(Ly)/Ny;
  kx = 2*np.pi*np.fft.fftfreq(X.shape[0])/dx
  ky = 2*np.pi*np.fft.fftfreq(X.shape[1])/dy
  Nyquist_kx = np.pi/dx; # maximum wavevector in x
  Nyquist_ky = np.pi/dy; # maximum wavevector in y
  x_mask = np.less_equal(abs(kx), Nyquist_kx/2.);
  y_mask = np.less_equal(abs(ky), Nyquist_ky/2.);
  mask = np.outer(x_mask, y_mask);

  # apply low-pass filter to noise
  for i in range(NIcomp):
    phi_rand_k = np.fft.fftn(phi_rand[i])
    phi_rand_k = phi_rand_k*mask;
    phi_rand[i] = np.fft.ifftn(phi_rand_k).real
    # ensure phi1 and phi2 are symmetric in y after adding noise
    phi_rand[i] = 0.5*(phi_rand[i] + phi_rand[i][:, ::-1])

  # apply a velocity perturbation to advect phi one step
  A=0.5; # amplitude of velocity perturbation
  lam_x = 2.*Lx/(np.linspace(1,8,8)[runid]);
  lam_y = Ly/4.;
  xc = 0. #Lx/2.;
  yc = 0. #Ly/2.;

  #U =  1j*Ky[n,m]*np.exp(1j*Kx[n,m]*X+1j*Ky[n,m]*Y)
  #V = -1j*Kx[n,m]*np.exp(1j*Kx[n,m]*X+1j*Ky[n,m]*Y)

  # in x and y
  vx = A*np.cos(2*np.pi*(X-xc)/lam_x)*np.cos(2*np.pi*(Y-yc)/lam_y)
  vy = A*np.sin(2*np.pi*(X-xc)/lam_x)*np.sin(2*np.pi*(Y-yc)/lam_y)
  vz = A*np.zeros(X.shape)

  # in x only
  #vx = A*np.cos(2*np.pi*(X-xc)/lam_x)
  #vy = A*np.zeros(Y.shape)
  #vz = A*np.zeros(X.shape)

  # now, advect phi_1 and phi_2 one step
  (Kx, Ky) = np.meshgrid(kx, ky, indexing='ij')
  for i in range(NIcomp):
    phi_k = np.fft.fftn(phi_0[i])
    phi_adv = - vx*np.fft.ifftn(1j*Kx*phi_k).real - vy*np.fft.ifftn(1j*Ky*phi_k).real
    #phi[i] = phi_0[i]+phi_rand[i]+phi_adv
    phi[i] = phi_0[i]+phi_adv

  phi_inlist = phi_in.tolist()
  phi_outlist = phi_out.tolist()
  # i/o
  ic_dict = {'IC_flag':IC_flag, 
             'phi_in':phi_inlist,
             'phi_out':phi_outlist,
             'sig':sig,
             'w':w,
             'f':f,
             'lam_x':lam_x,
             'lam_y':lam_y,
             'my_seed':my_seed}
# }}}
elif (IC_flag == 108): # read from file, add gradient
# {{{

  eq_dim=1

  phi_file = np.zeros(NIcomp, dtype=object)
  phi_0 = np.zeros(NIcomp, dtype=object)
  if (eq_dim == 1): # (for 1D equilibrium result)

    phi_file[0]='../equil/phi1.out'
    phi_file[1]='../equil/phi2.out'
    phi_file[2]='../equil/phi3.out'
  
    phi_1D = np.zeros(NIcomp, dtype=object)
    for i in range(NIcomp):
      phi_1D[i] = np.loadtxt(phi_file[i])[:,1];

    for i in range(NIcomp):
      phi_0[i] = np.tile(phi_1D[i], (Nx, 1))

  elif (eq_dim == 2): # (for 2D equilibrium result)

    phi_file[0]='../equil/phi1.out'
    phi_file[1]='../equil/phi2.out'
    phi_file[2]='../equil/phi3.out'

    tmp = np.loadtxt(phi_file[0]);
    #X = tmp[:, 0].reshape(Nx, Ny);
    #Y = tmp[:, 1].reshape(Nx, Ny);
    for i in range(NIcomp):
      phi_0[i] = np.loadtxt(phi_file[i])[:, 2].reshape(Nx, Ny);

  # -- paramters --
  f = 1./2. # fraction of inner region area/total area
  sig = 0.e-3 # noise strength
  # -- for piecewise gradient --
  dy = 0.50
  dphi = 0.15
  #dphi = -np.linspace(0.05, 0.2, 4)[runid];
  if (runid < 5):
    m = dphi/dy
  else:
    m = -dphi/dy
  # -- for step down --
  #dphi = 0.18 #np.linspace(0.05, 0.2, 4)[runid]; # a step of greater than 0.1 should cross the spinodal
  #dy1 = 0.02
  #dy2 = 0.08
  #dy3 = 0.20
  #m1 = -dphi/dy1
  #m3 = dphi/dy3

  y_half = np.linspace(0, 1, Ny/2)

  # simple linear gradient
  phi2_grad = np.zeros((Nx, Ny));
  phi2_grad[:, 0:Ny/2] = -m*(Y[:, 0:Ny/2]/Ly - 0.5 + f/2.)
  phi2_grad[:, Ny/2:] = m*(Y[:,Ny/2:]/Ly- 0.5 - f/2.)

  # outer domain (piecewise gradient)
  #phi_half = np.piecewise(y_half, \
  #  [y_half<f, (y_half>=f) & (y_half < (f+dy)), y_half >= (f+dy)], \
  #  [0,        lambda y_half: m*(y_half-f),     m*dy])

  # inner domain (piecewise gradient)
  #phi_half = np.piecewise(y_half, \
  #  [y_half<f-dy, ((y_half >= (f-dy)) & (y_half<f)), y_half >= f], 
  #  [m*dy, lambda y_half: -m*(y_half-f), 0])

  # inner and outer domain (piecewise gradient)
  #phi_half = np.piecewise(y_half, \
  #  [#
  #    y_half<f-dy,#
  #    (y_half >= (f-dy)) & (y_half<f),#
  #    y_half >= f,#
  #    y_half<f,#
  #    (y_half>=f) & (y_half < (f+dy)),#
  #    y_half >= (f+dy)#
  #  ],[#
  #    m*dy,#
  #    lambda y_half: -m*(y_half-f),#
  #    0,#
  #    0,#
  #    lambda y_half: m*(y_half-f),#
  #    m*dy#
  #  ])

  # inner domain (piecewise step down)
  #phi_half = np.piecewise(y_half, #
  #  [#  
  #    y_half<(f-dy1-dy2-dy3), #
  #    ((y_half >= (f-dy1-dy2-dy3)) & (y_half<(f-dy1-dy2))), #
  #    ((y_half >= (f-dy1-dy2)) & (y_half<(f-dy1))), #
  #    ((y_half >= (f-dy1)) & (y_half<f)), #
  #    y_half >= f # 
  #  ], # 
  #  [#  
  #    0,# 
  #    lambda y_half: m3*(y_half-f+dy1+dy2+dy3),#
  #    dphi,#
  #    lambda y_half: m1*(y_half-f),#
  #     0#  
  #  ])
  #phi_whole = np.concatenate((phi_half[::-1], phi_half))
  #phi2_grad = np.tile(phi_whole, (Nx, 1))

  # inner domain, new interfaces
  #phi_Leq = np.zeros(NIcomp);
  #phi_Req = np.zeros(NIcomp);
  #for i in range(NIcomp):
  #  phi_Leq[i] = phi_0[i][0, 0]
  #  phi_Req[i] = phi_0[i][0, -1]
  #y_L = abs(y_half-0.11).argmin()
  #y_R = abs(y_half-0.09).argmin()
  #phi_L = np.zeros(NIcomp);
  #phi_R = np.zeros(NIcomp);
  #for i in range(NIcomp):
  #  phi_L[i] =  np.concatenate( ( phi_0[i][:, y_L::]-phi_Leq[i],  np.ones((X.shape[0], y_L))*(phi_Req[i]-phi_Leq[i])), axis=1)
  #  phi_R[i] =  np.concatenate( (-phi_0[i][:, y_R::]+phi_Leq[i], -np.ones((X.shape[0], y_R))*(phi_Req[i]-phi_Leq[i])), axis=1)
  #phi_int = np.zeros(NIcomp);
  #for i in range(NIcomp):
  #  phi_int[i] = phi_L[i] + phi_R[i]

  for i in range(NIcomp):
    phi_0[i] = phi_0[i] #+ phi_int[i]
    if i==1:
      phi_0[i] = phi_0[i] + phi2_grad

  # random perturbation
  phi_rand = np.zeros(NIcomp);
  for i in range(NIcomp):
    phi_rand[i] = sig*(np.random.rand(*Y.shape)-0.5)

  # mask for high k-modes
  dx = float(Lx)/Nx;
  dy = float(Ly)/Ny;
  kx = 2*np.pi*np.fft.fftfreq(X.shape[0])/dx
  ky = 2*np.pi*np.fft.fftfreq(X.shape[1])/dy
  Nyquist_kx = np.pi/dx; # maximum wavevector in x
  Nyquist_ky = np.pi/dy; # maximum wavevector in y
  x_mask = np.less_equal(abs(kx), Nyquist_kx/2.);
  y_mask = np.less_equal(abs(ky), Nyquist_ky/2.);
  mask = np.outer(x_mask, y_mask);

  # apply low-pass filter to noise
  for i in range(NIcomp):
    phi_rand_k = np.fft.fftn(phi_rand[i])
    phi_rand_k = phi_rand_k*mask;
    phi_rand[i] = np.fft.ifftn(phi_rand_k).real

  # ensure phi1 and phi2 are symmetric in y after adding noise
  for i in range(NIcomp):
    phi_rand[i] = 0.5*(phi_rand[i] + phi_rand[i][:, ::-1])

  # apply a velocity perturbation to advect phi one step
  A=0.5; # amplitude of velocity perturbation
  lam_x = 2.*Lx/(np.array([1., 2., 4., 8., 16.])[runid%5]);
  lam_y = Ly/4.;
  xc = 0. #Lx/2.;
  yc = 0. #Ly/2.;

  #U =  1j*Ky[n,m]*np.exp(1j*Kx[n,m]*X+1j*Ky[n,m]*Y)
  #V = -1j*Kx[n,m]*np.exp(1j*Kx[n,m]*X+1j*Ky[n,m]*Y)

  # in x and y
  vx = A*np.cos(2*np.pi*(X-xc)/lam_x)*np.cos(2*np.pi*(Y-yc)/lam_y)
  vy = A*np.sin(2*np.pi*(X-xc)/lam_x)*np.sin(2*np.pi*(Y-yc)/lam_y)
  vz = A*np.zeros(X.shape)

  # in x only
  #vx = A*np.cos(2*np.pi*(X-xc)/lam_x)
  #vy = A*np.zeros(Y.shape)
  #vz = A*np.zeros(X.shape)

  # now, advect phi_1 and phi_2 one step
  (Kx, Ky) = np.meshgrid(kx, ky, indexing='ij')
  for i in range(NIcomp):
    phi_k = np.fft.fftn(phi_0[i])
    phi_adv = - vx*np.fft.ifftn(1j*Kx*phi_k).real - vy*np.fft.ifftn(1j*Ky*phi_k).real
    phi[i] = phi_0[i] + phi_adv[i] + phi_rand[i]

  phi_filelist = phi_file.tolist()
  # i/o
  ic_dict = {'IC_flag':IC_flag, 
             'phi_file':phi_filelist,
             'sig':sig,
             'f':f}
  #f1.write('%.16e\t'%m + '# m (slope for gradient)\n')
# }}}
elif (IC_flag == 109): # 1D sinusoid in x or y
# {{{

  #phi_array = np.linspace(0.49, 0.42, 5); #np.array([0.05, 0.15, 0.25, 0.35, 0.45])
  #phi1_array = np.loadtxt('params_PhiAv.in')[:, 0];
  #phi2_array = np.loadtxt('params_PhiAv.in')[:, 1];

  #params
  # note: phi1_0 and phi2_0 cannot be *exactly* equal, add 1e-3 or something
  #phi1_0 = phi1_array[1] # avg concentration of polymer 
  #phi2_0 = phi2_array[1] # avg concentration of non-solvent
  #phi1_0 = phi1_array[runid] # avg concentration of polymer
  #phi2_0 = phi2_array[runid] # avg concentration of non-solvent
  phi_0 = np.zeros(NIcomp);
  phi_0[0] = 0.25
  phi_0[1] = 0.15
  phi_0[2] = 0.10
  T = 2 # periods for sinusoid
  A = 0.1 # sinusoid amplitude

  # phi
  for i in range(NIcomp):
    phi[i] = A*np.sin(2*T*np.pi*X/Lx) + phi_0[i];
    #phi[i] = A*np.sin(2*T*np.pi*Y/Ly) + phi_0[i];

  phi_0list = phi_0.tolist()
  # i/o
  ic_dict = {'IC_flag':IC_flag, 
             'phi_0':phi_0list,
             'T':T,
             'A':A}
# }}}
elif (IC_flag == 110): # 2D Gaussian pulse
# {{{
  #params
  phi_0 = np.zeros(NIcomp);
  phi_0[0] = 0.25 # avg concentration of polymer
  phi_0[1] = 0.1 # avg concentration of non-solvent
  phi_0[2] = 0.1 # avg concentration of solvent
  sig2 = np.zeros(NIcomp);
  sig2[0] = 1. # variance for phi1
  sig2[1] = 1. # variance for phi2
  sig2[2] = 1. # variance for phi3
  A = 0.2 # gaussian amplitude

  # phi
  # normalized (unbounded)
  #for i in range(NIcomp):
    #phi[i] = A/(2*np.pi*sig2[i])*np.exp((X**2+Y**2)/(2*sig2[i])) + phi_0[i];

  # not-normalized ("A" dictates scale)
  #for i in range(NIcomp):
    #phi[i] = -A*np.exp(-((X-3*Lx/4)**2+(Y-Ly/2)**2)/(2*sig2[i])) + phi_0[i];

  # 1D gaussian (X)
  for i in range(NIcomp):
    phi[i] = -A*np.exp(-((X-3*Lx/4)**2)/(2*sig2[i])) + phi_0[i];

  # 1D gaussian (Y)
  #for i in range(NIcomp):
  #phi[i] = -A*np.exp(-((Y-1*Ly/4)**2)/(2*sig2[i])) + phi_0[i];

  phi_0list = phi_0.tolist()
  sig2list = sig2.tolist()
  # i/o
  ic_dict = {'IC_flag':IC_flag, 
             'phi_0':phi_0list,
             'sig2':sig2list,
             'A':A}
# }}}
elif (IC_flag == 111): # 3D sinusoid
# {{{
  #params
  phi_0 = np.zeros(NIcomp);
  phi_0[0] = 0.25 # avg concentration of polymer
  phi_0[1] = 0.60 # avg concentration of non-solvent
  phi_0[2] = 0.60 # avg concentration of solvent
  Tx = 1 # periods for sinusoid in x-dir
  Ty = 3 # periods for sinusoid in y-dir
  Tz = 5 # periods for sinusoid in y-dir
  A = 0.05 # sinusoid amplitude

  # phi
  for i in range(NIcomp):
    phi[i] = A*np.sin(2*Tx*np.pi*X/Lx)*np.sin(2*Ty*np.pi*Y/Ly)*np.sin(2*Tz*np.pi*Z/Lz) + phi_0[i];

  phi_0list = phi_0.tolist()
  # i/o
  ic_dict = {'IC_flag':IC_flag, 
             'phi_0':phi_0list,
             'Tx':Tx,
             'Ty':Ty,
             'Tz':Tz,
             'A':A}
# }}}
elif (IC_flag == 112): # 2D circle/ellipse
# {{{

  #eq_phi = np.loadtxt("phase_diagram_bin.in");
  #alpha = eq_phi[:, :3]
  #beta  = eq_phi[::-1, 3:]

  #params
  phi_in = np.zeros(NIcomp);
  phi_out = np.zeros(NIcomp);
  phi_in[0] = 0.03 #0.60 #alpha[79, 0]#0.01 #0.4 # avg concentration of polymer
  phi_in[1] = 0.86 #0.31 #alpha[79, 1]#0.98 #0.5 # avg concentration of non-solvent
  phi_in[2] = 0.86 # avg concentration of solvent
  phi_out[0] = 0.60 #0.03 #beta[79, 0]#0.25 #0.01 # avg concentration of polymer
  phi_out[1] = 0.03 #0.86 #beta[79, 1]#0.6 #0.98 # avg concentration of non-solvent
  phi_out[2] = 0.03 # avg concentration of solvent
  sig2 = 0.e-3 # noise variance
  w = 0.5 # thickness of tanh boundary
  Rx = 0.2*Lx # ellipse x-axis
  Ry = 0.2*Ly # ellipse y-axis 
  my_seed = np.random.randint(0, 4294967295) #4294... is the maximum seed value

  # basic phi step parallel to the x-axis
  X_c = (X - Lx/2)/Rx
  Y_c = (Y - Ly/2)/Ry
  dphi = np.zeros(NIcomp);
  phi_0 = np.zeros(NIcomp, dtype=object);
  for i in range(NIcomp):
    dphi[i] = phi_out[i] - phi_in[i]
    phi_0[i] = 0.5*dphi[i]*(np.tanh((X_c**2 + Y_c**2 - 1)/w) + 1) + phi_in[i]
  # perhaps make this elliptical? X^2/Rx^2 + Y^2/Ry^2 - 1 == 0

  # get random perturbation
  np.random.seed(my_seed)
  for i in range(NIcomp):
    phi_rand = sig2*(np.random.rand(*X.shape)-0.5)
    phi[i] = phi_0[i]+phi_rand

  phi_inlist = phi_in.tolist()
  phi_outlist = phi_out.tolist()
  # i/o
  ic_dict = {'IC_flag':IC_flag, 
             'phi_in':phi_inlist,
             'phi_out':phi_outlist,
             'sig2':sig2,
             'w':w,
             'Rx':Rx,
             'Ry':Ry,
             'my_seed':my_seed}
# }}}
# --- (200 Series) Block Polymers ---
elif (IC_flag == 200): # homogeneous initial condition
# {{{
  # params
  phi1_0 = 0.4     # homopolymer A
  phi2_0 = 0.1     # block A of diblock
  phi3_0 = phi2_0  # block B of diblock

  # phi
  phi[0] = np.zeros(X.shape) + phi1_0
  phi[1] = np.zeros(X.shape) + phi2_0
  phi[2] = np.zeros(X.shape) + phi3_0

  # i/o
  ic_dict = {'IC_flag':IC_flag,
             'phi1_0':phi1_0,
             'phi2_0':phi2_0,
             'phi3_0':phi3_0}
# }}}
elif (IC_flag == 201): # uniform noise about homogeneous initial conditions
# {{{
  # params
  phi1_0 = 0.4 #0.46
  phi2_0 = 0.5-phi1_0 #0.44 #0.45 + 0.05*runid
  phi3_0 = phi2_0
  sig = 2e-2 # std. deviation for random noise
  my_seed = 1234  #np.random.randint(0, 4294967295) #1234

  # phi
  np.random.seed(my_seed)
  phi1 = phi1_0*sig*np.random.rand(*X.shape)
  phi2 = phi2_0*sig*np.random.rand(*X.shape)
  phi3 = phi3_0*sig*np.random.rand(*X.shape)
  phi1_k = np.fft.fftn(phi1)
  phi2_k = np.fft.fftn(phi2)
  phi3_k = np.fft.fftn(phi3)
  phi1_k[np.unravel_index(0, X.shape)] = phi1_0*Ngrid # set 0 element to mean
  phi2_k[np.unravel_index(0, X.shape)] = phi2_0*Ngrid
  phi3_k[np.unravel_index(0, X.shape)] = phi3_0*Ngrid
  phi[0] = np.fft.ifftn(phi1_k).real
  phi[1] = np.fft.ifftn(phi2_k).real
  phi[2] = np.fft.ifftn(phi3_k).real

  # i/o
  ic_dict = {'IC_flag':IC_flag,
             'phi1_0':phi1_0,
             'phi2_0':phi2_0,
             'phi3_0':phi3_0,
             'sig':sig,
             'my_seed':my_seed}
# }}}
elif (IC_flag == 202): # guassian random initial conditions
# {{{ 
  pass
# }}}
elif (IC_flag == 203): # tanh profile in y
# {{{
  pass
# }}}
elif (IC_flag == 204): # sinusoids in x and y
# {{{
  pass
# }}}
elif (IC_flag == 205): # gradient in y-direction
# {{{
  pass
# }}}
elif (IC_flag == 206): # 1D tanh profile in x
# {{{

  #phi_params = np.loadtxt("../../setup/run"+str(runid)+"/params_Equilibrium.in")

  #params
  phi1_in = 0.94  # Homopolymer A volume fraction in inner region
  phi2_in = 0.5 - (0.02+phi1_in)/2  # Diblock's A volume fraction in inner region
  phi3_in = phi2_in  # Diblock's B volume fraction in inner region
  phi1_out = 0.02 # Homopolymer A volume fraction in outer region
  phi2_out = phi2_in # Diblock's A volume fraction in outer region
  phi3_out = phi3_in # Diblock's B volume fraction in outer region
  f = 0.5 #phi_params[4] #1./2. # fraction of inner region area/total area
  sig = 2.e-3 #0.5e-3 # noise strength
  w = 4 # boundary thickness
  my_seed = np.random.randint(0, 4294967295) #4294... is the maximum seed value

  # mask for high k-modes
  dx = float(Lx)/Nx;
  kx = 2*np.pi*np.fft.fftfreq(X.shape[0])/dx
  Nyquist_kx = np.pi/dx; # maximum wavevector in x
  x_mask = np.less_equal(abs(kx), Nyquist_kx/2.);

  # phi
  lb = 0.5 - 0.5*f; # lower bound
  ub = 0.5 + 0.5*f; # upper bound
  X_lb = (X-lb*Lx)/w
  X_ub = (ub*Lx-X)/w
  dphi1 = phi1_in - phi1_out
  dphi2 = phi2_in - phi2_out
  dphi3 = phi3_in - phi3_out
  np.random.seed(my_seed)

  phi1_noise = sig*np.random.rand(*X.shape)
  phi2_noise = sig*np.random.rand(*X.shape)
  phi3_noise = sig*np.random.rand(*X.shape)
   
  phi1 = 0.5*dphi1*np.tanh(X_lb) + 0.5*dphi1*np.tanh(X_ub) + phi1_out + phi1_noise
  phi2 = 0.5*dphi2*np.tanh(X_lb) + 0.5*dphi2*np.tanh(X_ub) + phi2_out + phi2_noise
  phi3 = 0.5*dphi3*np.tanh(X_lb) + 0.5*dphi3*np.tanh(X_ub) + phi3_out + phi3_noise
 
  # set mean to (phi_in + phi_out)/2 and apply low-pass filter
  phi1_k = np.fft.fft(phi1)
  phi2_k = np.fft.fft(phi2)
  phi3_k = np.fft.fft(phi3)
  phi1_k[np.unravel_index(0, X.shape)] = Ngrid*(phi1_in + phi1_out)/2
  phi2_k[np.unravel_index(0, X.shape)] = Ngrid*(phi2_in + phi2_out)/2
  phi3_k[np.unravel_index(0, X.shape)] = Ngrid*(phi3_in + phi3_out)/2
  phi1_k = phi1_k*x_mask;
  phi2_k = phi2_k*x_mask;
  phi3_k = phi3_k*x_mask;
  phi1 = np.fft.ifft(phi1_k).real
  phi2 = np.fft.ifft(phi2_k).real
  phi3 = np.fft.ifft(phi3_k).real
  # ensure phi1 and phi2 are symmetric in x after adding noise
  phi[0] = 0.5*(phi1 + phi1[::-1])
  phi[1] = 0.5*(phi2 + phi2[::-1])
  phi[2] = 0.5*(phi3 + phi3[::-1])

  # i/o
  ic_dict = {'IC_flag':IC_flag,
             'phi1_in':phi1_in,
             'phi2_in':phi2_in,
             'phi3_in':phi3_in,
             'phi1_out':phi1_out,
             'phi2_out':phi2_out,
             'phi3_out':phi3_out,
             'sig':sig,
             'w':w,
             'f':f,
             'my_seed':my_seed}
# }}}
elif (IC_flag == 207): # tanh profile in y with roll cells
# {{{
  pass
# }}}
elif (IC_flag == 208): # read from file, add velocity pertubation
# {{{

  eq_dim=1
  if (eq_dim == 1): # (for 1D equilibrium result)

    phi1_file='../equil/phi1.out'
    phi2_file='../equil/phi2.out'
    phi3_file='../equil/phi3.out'
     
    phi1_1D = np.loadtxt(phi1_file)[:,1];
    phi2_1D = np.loadtxt(phi2_file)[:,1];
    phi3_1D = np.loadtxt(phi3_file)[:,1];

    phi1_0 = np.tile(phi1_1D, (Nx, 1))
    phi2_0 = np.tile(phi2_1D, (Nx, 1))
    phi3_0 = np.tile(phi3_1D, (Nx, 1))

  elif (eq_dim == 2): # (for 2D equilibrium result)

    phi1_file='../equil/phi1.out'
    phi2_file='../equil/phi2.out'
    phi3_file='../equil/phi3.out'
    #tmp = np.loadtxt(phi1_file);
    #X = tmp[:, 0].reshape(Nx, Ny);
    #Y = tmp[:, 1].reshape(Nx, Ny);
    phi1_0 = np.loadtxt(phi1_file)[:, 2].reshape(Nx, Ny);
    phi2_0 = np.loadtxt(phi2_file)[:, 2].reshape(Nx, Ny);
    phi3_0 = np.loadtxt(phi3_file)[:, 2].reshape(Nx, Ny);

  # -- parameters --
  f = 1./2. # fraction of inner region area/total area
  sig = 0.e-3 # noise strength
  # -- for piecewise gradient --
  dy = 0.50
  dphi = 0
  #dphi = -np.linspace(0.05, 0.2, 4)[runid];
  if (runid < 5):
    m = dphi/dy
  else:
    m = -dphi/dy

  # -- for step down --
  #dphi = 0.18 #np.linspace(0.05, 0.2, 4)[runid]; # a step of greater than 0.1 should cross the spinodal
  #dy1 = 0.02
  #dy2 = 0.08
  #dy3 = 0.20
  #m1 = -dphi/dy1
  #m3 = dphi/dy3

  y_half = np.linspace(0, 1, Ny/2)

  # simple linear gradient
  phi2_grad = np.zeros((Nx, Ny));
  phi2_grad[:, 0:Ny/2] = -m*(Y[:, 0:Ny/2]/Ly - 0.5 + f/2.)
  phi2_grad[:, Ny/2:] = m*(Y[:,Ny/2:]/Ly- 0.5 - f/2.)

  # random perturbation
  phi1_rand = sig*(np.random.rand(*Y.shape)-0.5)
  phi2_rand = sig*(np.random.rand(*Y.shape)-0.5)
  phi3_rand = sig*(np.random.rand(*Y.shape)-0.5)

  # mask for high k-modes
  dx = float(Lx)/Nx;
  dy = float(Ly)/Ny;
  kx = 2*np.pi*np.fft.fftfreq(X.shape[0])/dx
  ky = 2*np.pi*np.fft.fftfreq(X.shape[1])/dy
  Nyquist_kx = np.pi/dx; # maximum wavevector in x
  Nyquist_ky = np.pi/dy; # maximum wavevector in y
  x_mask = np.less_equal(abs(kx), Nyquist_kx/2.);
  y_mask = np.less_equal(abs(ky), Nyquist_ky/2.);
  mask = np.outer(x_mask, y_mask);

  # apply low-pass filter to noise
  phi1_rand_k = np.fft.fftn(phi1_rand)
  phi2_rand_k = np.fft.fftn(phi2_rand)
  phi3_rand_k = np.fft.fftn(phi3_rand)
  phi1_rand_k = phi1_rand_k*mask;
  phi2_rand_k = phi2_rand_k*mask;
  phi3_rand_k = phi3_rand_k*mask;
  phi1_rand = np.fft.ifftn(phi1_rand_k).real
  phi2_rand = np.fft.ifftn(phi2_rand_k).real
  phi3_rand = np.fft.ifftn(phi3_rand_k).real

  # ensure phi1 and phi2 are symmetric in y after adding noise
  phi1_rand = 0.5*(phi1_rand + phi1_rand[:, ::-1])
  phi2_rand = 0.5*(phi2_rand + phi2_rand[:, ::-1])
  phi3_rand = 0.5*(phi3_rand + phi3_rand[:, ::-1])

  # apply a velocity perturbation to advect phi one step
  A=0; # amplitude of velocity perturbation
  lam_x = 2.*Lx/(np.array([1., 2., 4., 8., 16.])[runid%5]);
  lam_y = Ly/4.;
  xc = 0. #Lx/2.;
  yc = 0. #Ly/2.;

  #U =  1j*Ky[n,m]*np.exp(1j*Kx[n,m]*X+1j*Ky[n,m]*Y)
  #V = -1j*Kx[n,m]*np.exp(1j*Kx[n,m]*X+1j*Ky[n,m]*Y)

  # in x and y
  vx = A*np.cos(2*np.pi*(X-xc)/lam_x)*np.cos(2*np.pi*(Y-yc)/lam_y)
  vy = A*np.sin(2*np.pi*(X-xc)/lam_x)*np.sin(2*np.pi*(Y-yc)/lam_y)
  vz = A*np.zeros(X.shape)

  # in x only
  #vx = A*np.cos(2*np.pi*(X-xc)/lam_x)
  #vy = A*np.zeros(Y.shape)
  #vz = A*np.zeros(X.shape)

  # now, advect phi_1 and phi_2 one step
  (Kx, Ky) = np.meshgrid(kx, ky, indexing='ij')
  phi1_k = np.fft.fftn(phi1_0)
  phi2_k = np.fft.fftn(phi2_0)
  phi3_k = np.fft.fftn(phi3_0)  
  phi1_adv = - vx*np.fft.ifftn(1j*Kx*phi1_k).real - vy*np.fft.ifftn(1j*Ky*phi1_k).real
  phi2_adv = - vx*np.fft.ifftn(1j*Kx*phi2_k).real - vy*np.fft.ifftn(1j*Ky*phi2_k).real
  phi3_adv = - vx*np.fft.ifftn(1j*Kx*phi3_k).real - vy*np.fft.ifftn(1j*Ky*phi3_k).real
  phi[0] = phi1_0 + phi1_adv + phi1_rand
  phi[1] = phi2_0 + phi2_adv + phi2_rand
  phi[2] = phi3_0 + phi3_adv + phi3_rand

  # i/o
  ic_dict = {'IC_flag':IC_flag,
             'phi1_file':phi1_file,
             'phi2_file':phi2_file,
             'phi3_file':phi3_file,
             'sig':sig,
             'm':m,
             'f':f}
# }}}
elif (IC_flag == 209): # 1D sinusoid in x or y
# {{{

  #params
  phi1_0 = 0.4
  phi2_0 = 0.5 - phi1_0
  phi3_0 = phi2_0

  T = 2 # periods for sinusoid
  A = 0.3 # sinusoid amplitude

  # phi
  phi[0] = (A*np.sin(2*T*np.pi*X/Lx) + 1)*phi1_0;
  phi[1] = (A*np.sin(2*T*np.pi*X/Lx) + 1)*phi2_0;
  phi[2] = (A*np.sin(2*T*np.pi*X/Lx) + 1)*phi3_0;
  #phi1 = A*np.sin(2*T*np.pi*Y/Ly) + phi1_0;
  #phi2 = A*np.sin(2*T*np.pi*Y/Ly) + phi2_0;

  # i/o
  ic_dict = {'IC_flag':IC_flag,
             'phi1_0':phi1_0,
             'phi2_0':phi2_0,
             'phi3_0':phi3_0,
             'T':T,
             'A':A}
# }}}
elif (IC_flag == 210): # 2D Gaussian pulse
# {{{
  pass
# }}}
elif (IC_flag == 211): # 3D sinusoid
# {{{
  pass
# }}}
elif (IC_flag == 212): # 2D circle/ellipse
# {{{

  #params
  phi1_in = 0.97 #0.60 #alpha[79, 0]#0.01 #0.4 # avg concentration of polymer
  phi2_in = 0.01 #0.31 #alpha[79, 1]#0.98 #0.5 # avg concentration of non-solvent
  phi3_in = phi2_in
  phi1_out = 1-phi1_in-phi2_in-phi3_in #0.03 #beta[79, 0]#0.25 #0.01 # avg concentration of polymer
  phi2_out = phi2_in #0.86 #beta[79, 1]#0.6 #0.98 # avg concentration of non-solvent
  phi3_out = phi2_out
  sig2 = 0.e-3 # noise variance
  w = 0.75 # thickness of tanh boundary
  Rx = Lx/np.sqrt(2*np.pi)/1.5 #ellipse x-axis
  Ry = Rx #Lx/np.sqrt(2*np.pi)/2 # ellipse y-axis 
  my_seed = np.random.randint(0, 4294967295) #4294... is the maximum seed value

  # basic phi step parallel to the x-axis
  X_c = (X - Lx/2)/Rx
  Y_c = (Y - Ly/2)/Ry
  dphi1 = phi1_out - phi1_in
  dphi2 = phi2_out - phi2_in
  dphi3 = phi3_out - phi3_in
  phi1_0 = 0.5*dphi1*(np.tanh((X_c**2 + (Y_c-Ly/4/Ry)**2 - 1)/w) + np.tanh((X_c**2 + (Y_c+Ly/4/Ry)**2 - 1)/w)) + phi1_in
  phi2_0 = 0.5*dphi2*(np.tanh((X_c**2 + Y_c**2 - 1)/w) + 1) + phi2_in
  phi3_0 = 0.5*dphi3*(np.tanh((X_c**2 + Y_c**2 - 1)/w) + 1) + phi3_in
  # perhaps make this elliptical? X^2/Rx^2 + Y^2/Ry^2 - 1 == 0

  # get random perturbation
  np.random.seed(my_seed)
  phi1_rand = sig2*(np.random.rand(*X.shape)-0.5)
  phi2_rand = sig2*(np.random.rand(*X.shape)-0.5)
  phi3_rand = sig2*(np.random.rand(*X.shape)-0.5)
  phi[0] = phi1_0+phi1_rand
  phi[1] = phi2_0+phi2_rand
  phi[2] = phi3_0+phi3_rand

  # i/o
  ic_dict = {'IC_flag':IC_flag,
             'phi1_in':phi1_in,
             'phi2_in':phi2_in,
             'phi3_in':phi3_in,
             'phi1_out':phi1_out,
             'phi2_out':phi2_out,
             'phi3_out':phi3_out,
             'sig2':sig2,
             'w':w,
             'Rx':Rx,
             'Ry':Ry,
             'my_seed':my_seed}
# }}}

# ----------------------------------------------
#  Pick which type of velocity IC
# ----------------------------------------------
# 0: zero velocity
# 1: 2D sinusoids

VelIC_flag = 0

if (VelIC_flag == 0): # zero velocity
# {{{

  # vel
  vx = np.zeros(X.shape)
  if (dim > 1): vy = np.zeros(X.shape)
  if (dim > 2): vz = np.zeros(X.shape)

  # i/o
  vel_dict = {'VelIC_flag':VelIC_flag}

# }}}
elif (VelIC_flag == 1): # sinusoids
# {{{

  A = np.sqrt(2)/2*1e-1 # sinusoid amplitude
  Tx = 1 # periods for sinusoid in x-dir
  Ty = 3 # periods for sinusoid in y-dir

  # vel
  vx = A*np.sin(2*Tx*np.pi*X/Lx)*np.sin(2*Ty*np.pi*Y/Ly)
  if (dim > 1): vy = A*np.sin(2*Tx*np.pi*X/Lx)*np.sin(2*Ty*np.pi*Y/Ly)
  if (dim > 2): vz = A*np.sin(2*Tx*np.pi*X/Lx)*np.sin(2*Ty*np.pi*Y/Ly)

  # i/o
  vel_dict = {'VelIC_flag':VelIC_flag,
              'A':A,
              'Tx':Tx,
              'Ty':Ty}
# }}}

# ----------------------------------------------
#  Define distribution of body force for flow
# ----------------------------------------------
# Pick the type of flow based on the body force
# (Only valid for Model H and ModelH_gradmu)
# 0: shear flow
# 1: planar extensional flow

BF_FlowType_flag = 0

if (BF_FlowType_flag == 0): # shear flow
# {{{

  # shear rate = bf * delta_y / (2 * viscosity)
  bf = 1.;   # magnitude of body force
  bfx = np.zeros(X.shape)
  if (dim > 1): bfx = bfx + bf*((Y==0.))-bf*((Y==Ly/2))
  if (dim > 1): bfy = np.zeros(X.shape)
  if (dim > 2): bfz = np.zeros(X.shape)

# }}}
elif (BF_FlowType_Flag == 1): # extensional flow
# {{{

  # ext rate   = bf * delta_y / viscosity
  bf = 1.;   # magnitude of body force
  bfx = np.zeros(X.shape)

  if (dim > 1):
    bfx = bfx - bf*(np.abs(X-Lx/2)-Lx/4)*( \
             ((Y==Ly/Ny)+(Y==(Ly/2-Ly/Ny))).astype(int) \
            -((Y==(Ly/2+Ly/Ny))+(Y==(Ly-Ly/Ny))).astype(int))

    bfy = np.zeros(X.shape)
    bfy = bfy + bf*((Y==0.)*(X<Lx/2)*(X!=0.)+(Y==Ly/2)*(X>Lx/2)) \
              - bf*((Y==0)*(X>Lx/2)+(Y==Ly/2)*(X<Lx/2)*(X!=0.))
    bfy = bfy + bf*(np.abs(Y-Ly/2)-Ly/4)*( \
             ((X==Lx/Nx)+(X==(Lx/2-Lx/Nx))).astype(int) \
            -((X==(Lx/2+Lx/Nx))+(X==(Lx-Lx/Nx))).astype(int))

  if (dim > 2): bfz = np.zeros(X.shape)

# }}}

invert = np.ones(Y.shape)
#invert[:,:Ny//2,:] *=1
# ------------------------
# write ICs
# ------------------------

# {{{

json_dict.update(ic_dict)
json_dict.update(vel_dict)
f1 = open('params_ICs.in', 'w')
json.dump(json_dict,f1,indent=4)
f1.close()

# verify inital conditions are valid
phi_N = np.ones(X.shape)
for i in range(NIcomp):
  phi_N -= phi[i]

if ( np.any(phi_N < 0.) or np.any(phi_N > 1.) ):
  print('  (warning) I.C.s are not valid')

for i in range(NIcomp):
  if ( np.any(phi[i] < 0.) or np.any(phi[i] > 1.) ):
    print('  (warning) I.C.s are not valid')

def write_fmt_ps_field(field, name):
# {{{
  f1 = open(name, 'wb')
  f1.write(b'# Format version 3\n')
  f1.write(b'# nfields = 1\n')
  f1.write(b'# NDim = %d\n'%dim)
  if (dim == 1): f1.write(b'# PW grid = %d\n'%Nx)
  if (dim == 2): f1.write(b'# PW grid = %d'%Nx+b' %d\n'%Ny)
  if (dim == 3): f1.write(b'# PW grid = %d'%Nx+b' %d'%Ny+b' %d\n'%Nz)
  f1.write(b'# k-space data = 0 , complex data = 0\n')
  if (dim == 1): f1.write(b'# Columns: x field.real\n')
  if (dim == 2): f1.write(b'# Columns: x y field.real\n')
  if (dim == 3): f1.write(b'# Columns: x y z field.real\n')

  if (dim == 1):
    np.savetxt(f1, np.vstack([X, field]).T, fmt='%+0.4e %+0.10e')
  
  if (dim == 2):
    for i in range(Nx):
      np.savetxt(f1, np.vstack([X[i, :], Y[i, :], field[i, :]]).T, fmt='%+0.4e %+0.4e %+0.10e')
      if (i < Nx-1): f1.write(b'\n')

  if (dim == 3):
    for i in range(Nx):
      for j in range(Ny):
        np.savetxt(f1, np.vstack([X[i, j, :], Y[i, j, :], Z[i, j, :], field[i, j, :]]).T, fmt='%+0.4e %+0.4e %+0.4e %+0.10e')
      if (i < Nx-1): f1.write(b'\n')
  f1.close()
# }}}

def write_bin_ps_field(field, name):
# {{{
  version = 51 #3
  nfields = 1
  kspacedata=False
  IsComplex=True
  griddim = [Nx, Ny, Nz]
  boxtensor = np.zeros((3,3))
  boxtensor[0,0] = Lx
  boxtensor[1,1] = Ly
  boxtensor[2,2] = Lz
  elementsize = 16 # element size, number of bytes per element (8 for double * 2 for complex)

  outfile = open(name, "wb")
  outfile.write(struct.pack("@9s", b"FieldBin\0")) # need null character at end for c-style string
  outfile.write(struct.pack("@I", version))
  outfile.write(struct.pack("@i", nfields))
  outfile.write(struct.pack("@I", dim))
  for i in range(dim): 
    outfile.write(struct.pack("@L", griddim[i]))
  outfile.write(struct.pack("@?", kspacedata))
  outfile.write(struct.pack("@?", IsComplex))
  for i in range(dim):
    for j in range(dim):
      outfile.write(struct.pack("@d", boxtensor[i, j]))
  outfile.write(struct.pack("@L", elementsize))

  if dim==1:
    for i in range(Nx):
      outfile.write(struct.pack("@2d", field[i], 0.))
  if dim==2:
    for i in range(Nx):
      for j in range(Ny):
        outfile.write(struct.pack("@2d", field[i,j], 0.))
  if dim==3:
    for i in range(Nx):
      for j in range(Ny):
        for k in range(Nz):
          outfile.write(struct.pack("@2d", field[i,j,k], 0.))

  outfile.close()
# }}}

def write_fmt_fd_field(field, name):
# {{{
  f1 = open(name, 'wb')
  f1.write(b'# Format version 3\n')
  f1.write(b'# nfields = %d\n'%Nx)
  f1.write(b'# NDim = %d\n'%(dim-1))
  if (dim == 2): f1.write(b'# PW grid = %d\n'%Ny)
  if (dim == 3): f1.write(b'# PW grid = %d'%Ny+b' %d\n'%Nz)
  f1.write(b'# k-space data = 0 , complex data = 0\n')
  if (dim == 2): f1.write(b'# Columns: x fielddata\n')
  if (dim == 3): f1.write(b'# Columns: x y fielddata\n')

  if (dim == 2):
    myfmt = '%+0.4e' + ' %+0.16e'*Nx
    np.savetxt(f1, np.vstack([Y[0, :], field]).T, fmt=myfmt) 

  if (dim == 3):
    for j in range(Ny):
      myfmt = '%+0.4e %+0.4e' + ' %+0.16e'*Nx
      np.savetxt(f1, np.vstack([Y[0, j, :], Z[0, j, :], field[:, j, :]]).T, fmt=myfmt)
      if (j < Ny-1): f1.write(b'\n')
  f1.close()
# }}}

def write_bin_fd_field(field, name):
# {{{
  version = 51 #3
  nfields = Nx
  kspacedata=False
  IsComplex=True
  fddim = dim-1
  griddim = [Ny, Nz]
  boxtensor = np.zeros((2,2))
  boxtensor[0,0] = Ly
  boxtensor[1,1] = Lz
  elementsize = 16 # element size, number of bytes per element (8 for double * 2 for complex)

  outfile = open(name, "wb")
  outfile.write(struct.pack("@9s", b"FieldBin\0")) # need null character at end for c-style string
  outfile.write(struct.pack("@I", version))
  outfile.write(struct.pack("@i", nfields))
  outfile.write(struct.pack("@I", fddim))
  for i in range(fddim): 
    outfile.write(struct.pack("@L", griddim[i]))
  outfile.write(struct.pack("@?", kspacedata))
  outfile.write(struct.pack("@?", IsComplex))
  for i in range(fddim):
    for j in range(fddim):
      outfile.write(struct.pack("@d", boxtensor[i, j]))
  outfile.write(struct.pack("@L", elementsize))

  if dim==1:
    for i in range(Nx):
      outfile.write(struct.pack("@2d", field[i], 0.))
  if dim==2:
    for i in range(Nx):
      for j in range(Ny):
        outfile.write(struct.pack("@2d", field[i,j], 0.))
  if dim==3:
    for i in range(Nx):
      for j in range(Ny):
        for k in range(Nz):
          outfile.write(struct.pack("@2d", field[i,j,k], 0.))

  outfile.close()
# }}}

if (PS_flag == 0 and BinaryIn_flag == 1):
  for i in range(NIcomp):
    write_bin_ps_field(phi[i], 'phi%d.in'%(i+1))
  if (WriteVel_flag == 1):
    write_bin_ps_field(vx, 'vx.in')
    if (dim > 1): write_bin_ps_field(vy, 'vy.in')
    if (dim > 2): write_bin_ps_field(vz, 'vz.in')
  if (BodyForceFlag == 1):
    write_bin_ps_field(bfx, 'bfx.in')
    if (dim > 1): write_bin_ps_field(bfy, 'bfy.in')
    if (dim > 2): write_bin_ps_field(bfz, 'bfz.in')
  if (BodyForceFlag == 2):
    write_bin_ps_field(invert,'invert.in')

elif (PS_flag == 0 and BinaryIn_flag == 0):
  for i in range(NIcomp):
    write_fmt_ps_field(phi[i], 'phi%d.in'%(i+1))
  if (WriteVel_flag == 1):
    write_fmt_ps_field(vx, 'vx.in')
    if (dim > 1): write_fmt_ps_field(vy, 'vy.in')
    if (dim > 2): write_fmt_ps_field(vz, 'vz.in')
  if (BodyForceFlag == 1):
    write_fmt_ps_field(bfx, 'bfx.in')
    if (dim > 1): write_fmt_ps_field(bfy, 'bfy.in')
    if (dim > 2): write_fmt_ps_field(bfz, 'bfz.in')
  if (BodyForceFlag == 2):
    write_fmt_ps_field(invert,'invert.in')

elif (PS_flag == 1 and BinaryIn_flag == 1):
  for i in range(NIcomp):
    write_bin_fd_field(phi[i], 'phi%d.in'%(i+1))
  if (WriteVel_flag == 1):
    write_bin_fd_field(vx, 'vx.in')
    if (dim > 1): write_bin_fd_field(vy, 'vy.in')
    if (dim > 2): write_bin_fd_field(vz, 'vz.in')
  if (BodyForceFlag == 1):
    write_bin_fd_field(bfx, 'bfx.in')
    if (dim > 1): write_bin_fd_field(bfy, 'bfy.in')
    if (dim > 2): write_bin_fd_field(bfz, 'bfz.in')
  if (BodyForceFlag == 2):
    write_bin_fd_field(invert,'invert.in')

elif (PS_flag == 1 and BinaryIn_flag == 0):
  for i in range(NIcomp):
    write_fmt_fd_field(phi[i], 'phi%d.in'%(i+1))
  if (WriteVel_flag == 1):
    write_fmt_fd_field(vx, 'vx.in')
    if (dim > 1): write_fmt_fd_field(vy, 'vy.in')
    if (dim > 2): write_fmt_fd_field(vz, 'vz.in')
  if (BodyForceFlag == 1):
    write_fmt_fd_field(bfx, 'bfx.in')
    if (dim > 1): write_fmt_fd_field(bfy, 'bfy.in')
    if (dim > 2): write_fmt_fd_field(bfz, 'bfz.in')
  if (BodyForceFlag == 2):
    write_fmt_fd_field(invert,'invert.in')

# }}}

