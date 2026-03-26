#!/usr/bin/env python3

import numpy as np
import sys
import os
import json

if (len(sys.argv) >= 2):
  runid = int(sys.argv[1]); # should be an integer
else:
  runid = 0 # if nothing supplied, run it for zero

print('Writing paramters for runid: ', runid)

# ----------------------------------------------
#  No. of components (for all Models) and BCs
# ----------------------------------------------
# (Written in params_Keys.in at the bottom)

Ncomp = 3 
NIcomp = Ncomp-1
NIcomp_max = 4 # largest number of independent components possible (currently corresponds to 4-component systems)

# ----------------------------------------------
#  Parallel Parameters
# ----------------------------------------------
# {{{

# MPI
MPInnodes = 1 # Set to 1. Number of MPI nodes is determined at submission time.

# OpenMP
OMPnthreads = 1 # Set to 1. Number of OMP threads is determined at submission time.

# GPU
CUDAdeviceID = -1  # Set to -1. GPU device ID is determined at submission time.
CUDAnthreads = 128 # number of threads per block (block size)
maxnumelements = int(1e8) # max. total number of threads

#i/o
print('  - writing Parallel parameters (Nnodes = ', \
  MPInnodes, ', Nthreads = ', OMPnthreads, ', GPU = ', CUDAdeviceID, ')')
f1 = open('params_Parallel.in', 'w')
json.dump({'MPInnodes':MPInnodes, \
           'OMPnthreads':OMPnthreads, \
           'CUDAdeviceID':CUDAdeviceID, \
           'CUDAnthreads':CUDAnthreads, \
           'maxnumelements':maxnumelements},f1,indent=4)
f1.close()
# }}}

# ----------------------------------------------
#  Grid Parameters
# ----------------------------------------------
# {{{

# Grid flag
# 0: Pseudo-spectral
# 1: Hybrid finite-difference/pseudospectral
GridFlag = 0

# number of dimensions (trumps all below)
dim = 1

# box length
Lx = 32
Ly = 32
Lz = 32

# grid points
# For PS methods, one should try to use grids that are a factor of 2.
# This increases the efficiency of the FFTs.
#
# For a FD method (Nx only), one can use any size. However, remember 
# that dx = Lx/(Nx-1) for a non-periodic BC. So, often one might want to use:
# Nx_{FD} = Nx_{PS} + 1

Nx = int(6*Lx)
Ny = int(6*Ly)
Nz = int(6*Lz)

#i/o
print('  - writing Grid parameters (Dim = ', dim, \
      ': ', Nx, 'x', Ny, 'x', Nz, ')')
f1 = open('params_Grid.in', 'w')
json.dump({'GridFlag':GridFlag, \
           'dim':dim, \
           'Lx':Lx, \
           'Ly':Ly, \
           'Lz':Lz, \
           'Nx':Nx, \
           'Ny':Ny, \
           'Nz':Nz},f1,indent=4)
f1.close()
# }}}

# ----------------------------------------------
#  Free Energy Model Parameters
# ----------------------------------------------
# 0: N component, Linear chemical potential
# 1: N component, Flory--Huggins--de Gennes (FHG)
# 2: N > 4 component, Uneyama-Doi for A + sym. C-D + B + ...
# 3: N > 4 component, Uneyama-Doi for A + non-sym. C-D + B ...
# 4: N component, FHG - Kappa3 fix
#
# In 2 and 3: component 0 is a homopolymer
#             components 1 & 2 form a linear diblock
#             all other components including implicit are homopolymers

EnergyModelFlag = 1

print('  - writing Energy Model parameters (', EnergyModelFlag, ')')

if (EnergyModelFlag == 0):
# {{{

  # Hessian matrix
  Hess = np.zeros((NIcomp, NIcomp));
  for i in range(NIcomp):
    Hess[i,i] = 2.;
  
  # Gradient coefficients
  Kappa = np.zeros((NIcomp));
  for i in range(NIcomp):
    Kappa[i] = 1.;
  Hesslist = Hess.tolist()
  Kappalist = Kappa.tolist()
  f1 = open('params_EnergyModel.in', 'w')
  json.dump({'EnergyModelFlag':EnergyModelFlag, \
  #for i in range(NIcomp):
   # for j in range(NIcomp):
  'Hess':Hesslist, \
  #for i in range(NIcomp):
  'Kappa':Kappalist},f1,indent=4)
  f1.close()

# }}}
if (EnergyModelFlag == 1):
# {{{

  if (Ncomp < 2):
    print('***       (Error)        ***')
    print('*** Must have Ncomp >= 2 ***')

  Nr = 1 # reference N
  C_reg = 0.e-3 # regularization strength parameter
  delta = 1.e-3 # regularization scale parameter

  N = np.zeros(NIcomp_max);
  Chi = np.zeros((NIcomp_max, NIcomp_max));
  Kappa = np.zeros((NIcomp));

  # Add degree of polymerization here
  N[0] = 8. # polymer
  N[1] = 1. # non-solvent
  N[2] = 1. # solvent
  #N[3] = 1. # additive

  # Add Flory Huggins Interaction Parameters here
  Chi_c = 0.5*(np.sqrt(1./N[0])+np.sqrt(1./N[1]))**2 # binary critical ChiN for polymer/non-solvent

  Chi[0,1] = 1.01#1.2*Chi_c # polymer/non-solvent
  Chi[0,2] = 0.04 # polymer/solvent
  Chi[1,2] = 0.84 # non-solvent/solvent
  #Chi[1,3] = 0. # non-solvent/additive
  #Chi[2,3] = 0. # solvent/additive

  # Modify gradient coefficients here
  w = 2.0 # approx width in Rg
  prefactor = 2.4

  for i in range(NIcomp):
    Kappa[i] = 2.; #Chi[0,1]*Nr * (w/prefactor)**2 * (Chi[0,1]/Chi_c-1)
  
  Nlist = N.tolist()
  Chilist = Chi.tolist()
  Kappalist = Kappa.tolist()
  f1 = open('params_EnergyModel.in', 'w')
  json.dump({'EnergyModelFlag':EnergyModelFlag, \
  # (Flory Huggins + Square Gradient)\n')
             'Nr':Nr, \
  #(degree of polymerization for reference)\n')
             'C_reg':C_reg, \
  # (Regularization strength parameter)\n')
             'delta':delta, \
  # (Regularization scale parameter)\n')
  #for i in range(Ncomp):
             'N':Nlist, \
  #%d (deg. of polymerization for %d)\n'%(i,i))
  #for i in range(Ncomp):
    #for j in range(i+1, Ncomp):
             'Chi':Chilist, \
  #%d%d (Flory parameter between %d and %d)\n'%(i,j,i,j))
  #for i in range(NIcomp):
             'Kappa':Kappalist},f1,indent=4) 
  #grad coeff (%d)\n'%(i))
  f1.close()

# }}}
if ((EnergyModelFlag == 2) or (EnergyModelFlag == 3)):
# {{{

  if (Ncomp < 2):
    print('***       (Error)        ***')
    print('*** Must have Ncomp >= 2 ***')

  Nr = 40 # reference N
  C_reg = 0.e-5 # regularization strength parameter
  delta = 1.e-3 # regularization scale parameter

  N = np.zeros(Ncomp);
  Chi = np.zeros((Ncomp, Ncomp));
  Kappa = np.zeros((Ncomp));
  C = np.zeros((Ncomp, Ncomp));

  # Add degree of polymerization here
  N[0] = Nr # homopolymer A
  N[1] = Nr/2 # diblock's A
  N[2] = Nr/2 # diblock's B
  N[3] = Nr # homopolymer B

  # Add Flory Huggins Interaction Parameters here
  Chi_c = 0.5*(np.sqrt(1./N[0])+np.sqrt(1./N[1]))**2 # ignore this --  binary critical ChiN for polymer/non-solvent

  Chi[0,1] = 0    # A-h/A-d
  Chi[0,2] = 0.15  # A-h/B-d
  Chi[0,3] = Chi[0,2]  # A-h/B-h
  Chi[1,2] = Chi[0,2]  # A-d/B-d
  Chi[1,3] = Chi[0,2]  # A-d/B-h
  Chi[2,3] = 0    # B-h/B-d

  # Modify gradient coefficients here
  w = 2.0 # approx width in Rg
  prefactor = 2.4

  #for i in range(Ncomp):
  Kappa[0] = 0.25; #Chi[0,1]*Nr * (w/prefactor)**2 * (Chi[0,1]/Chi_c-1)
  Kappa[1] = 0.25;
  Kappa[2] = 0.25;
  Kappa[3] = 0.25;

  # Modify Coulomb coefficients here
  C[1,1] = 3 
  C[1,2] = -C[1,1]
  C[2,2] = C[1,1]
    
  
  Nlist = N.tolist()
  Chilist = Chi.tolist()
  Kappalist = Kappa.tolist()
  Clist = C.tolist()
  f1 = open('params_EnergyModel.in', 'w')
  json.dump({'EnergyModelFlag':EnergyModelFlag, \
   #(Flory Huggins + Square Gradient)\n')
             'Nr':Nr, \
   # (degree of polymerization for reference)\n')
             'C_reg':C_reg, \
   # (Regularization strength parameter)\n')
             'delta':delta, \
   # (Regularization scale parameter)\n')
  #for i in range(Ncomp):
             'N':Nlist, \
   #%d (deg. of polymerization for %d)\n'%(i,i))
  #for i in range(Ncomp):
   # for j in range(i+1, Ncomp):
             'Chi':Chilist, \
   #%d%d (Flory parameter between %d and %d)\n'%(i,j,i,j))
  #for i in range(Ncomp):
             'Kappa':Kappalist, \
   #grad coeff (%d)\n'%(i))
  #for i in range(Ncomp):
    #for j in range(i,Ncomp):
             'C':Clist},f1,indent=4) # C%d%d (Coulomb coefficient between %d and %d)\n'%(i,j,i,j))
  f1.close()
# }}}

if (EnergyModelFlag == 4):
# {{{

  if (Ncomp < 2):
    print('***       (Error)        ***')
    print('*** Must have Ncomp >= 2 ***')

  Nr = 1 # reference N
  C_reg = 0.e-3 # regularization strength parameter
  delta = 0.e-3 # regularization scale parameter

  N = np.zeros(NIcomp_max);
  Chi = np.zeros((NIcomp_max, NIcomp_max));
  Kappa = np.ones((NIcomp, NIcomp))

  # Add degree of polymerization here
  N[0] = 50. # polymer
  N[1] = 1. # non-solvent
  N[2] = 1. # solvent
  #N[3] = 1. # additive

  # Add Flory Huggins Interaction Parameters here
  Chi_c = 0.5*(np.sqrt(1./N[0])+np.sqrt(1./N[1]))**2 # binary critical ChiN for polymer/non-solvent
#1.0084848301906106 0.04164476303364435 0.8445049657627101
  Chi[0,1] = 1.01*Chi_c # polymer/non-solvent
  Chi[0,2] = 0.04 # polymer/solvent
  Chi[1,2] = 0.84 # non-solvent/solvent
  #Chi[1,3] = 0. # non-solvent/additive
  #Chi[2,3] = 0. # solvent/additive

  # Modify gradient coefficients here
  w = 2.0 # approx width in Rg
  prefactor = 2.4
  # Kappa = Chi[0,1]*Nr * (w/prefactor)**2 * (Chi[0,1]/Chi_c-1)
  k0 = 0.25; #polymer gradient coefficient 
  k1 = 0.25; #non-solvent gradient coefficient 
  k2 = 0.25; #solvent gradient coefficient 
  # Gradient Coefficient Matrix (off-diagonal are assumed to be the same!)
  # We evaluate a binary Kappa since phi2 = 1 - phi0 - phi1
  Kappa[0,0] = k0 + k2
  Kappa[0,1] = k2
  Kappa[1,0] = Kappa[0,1]
  Kappa[1,1] = k1 + k2

  
  Nlist = N.tolist()
  Chilist = Chi.tolist()
  Kappalist = Kappa.tolist()
  f1 = open('params_EnergyModel.in', 'w')
  json.dump({'EnergyModelFlag':EnergyModelFlag, \
   #((FHG - Kappa3 fix)\n')
             'Nr':Nr, \
   # (degree of polymerization for reference)\n')
             'C_reg':C_reg, \
   # (Regularization strength parameter)\n')
             'delta':delta, \
   # (Regularization scale parameter)\n')
  #for i in range(Ncomp):
             'N':Nlist, \
   #%d (deg. of polymerization for %d)\n'%(i,i))
  #for i in range(Ncomp):
   # for j in range(i+1, Ncomp):
             'Chi':Chilist, \
   #%d%d (Flory parameter between %d and %d)\n'%(i,j,i,j))
  #for i in range(Ncomp):
             'Kappa':Kappalist},f1,indent=4)
   #grad coeff (%d)\n'%(i))
  f1.close()
# }}}

# ----------------------------------------------
#  Mobility Model Parameters
# ----------------------------------------------
# 0: Constant
# 1: Rouse model
# 2: Rouse model scaled by the viscosity

MobilityModelFlag = 2

print('  - writing Mobility Model parameters (', MobilityModelFlag, ')')

if (MobilityModelFlag == 0):
# {{{

  # Mobility matrix
  Mob = np.zeros((NIcomp, NIcomp));
  for i in range(NIcomp):
    Mob[i,i] = 1.;
  Moblist = Mob.tolist()
  f1 = open('params_MobilityModel.in', 'w')
  json.dump({'MobilityModelFlag':MobilityModelFlag, \
  #for i in range(NIcomp):
   # for j in range(NIcomp):
             'Mob':Moblist},f1,indent=4) #mobility coeff (%d,%d)\n'%(i,j))
  f1.close()

# }}}
if (MobilityModelFlag == 1):
# {{{

  f1 = open('params_MobilityModel.in', 'w')
  json.dump({'MobilityModelFlag':MobilityModelFlag},f1,indent=4)# (Rouse model)\n')
  f1.close()

# }}}
if (MobilityModelFlag == 2):
# {{{

  # Mobility Prefactor Matrix (Should be symmetric!)
  C = np.ones((NIcomp, NIcomp)) # should be ones unless you know what you are doing
  C[0,0] = 1.0
  C[0,1] = 1.0
  C[1,0] = 1.0
  C[1,1] = 1.0
 
  Clist = C.tolist()
  f1 = open('params_MobilityModel.in', 'w')
  json.dump({'MobilityModelFlag':MobilityModelFlag, \
# (Rouse model scaled by viscosity.)\n')
 # for i in range(NIcomp):
  #  for j in range(NIcomp):
             'C':Clist},f1,indent=4)# mobility coeff (%d,%d)\n'%(i,j))
  f1.close()

# }}}

# ----------------------------------------------
#  Viscosity Model Parameters
# ----------------------------------------------
# 0: Linear viscosity
# 1: Exponential viscosity
# 2: Sigmoidal viscosity
# 3: VFTH viscosity model

ViscModelFlag=1

print('  - writing Viscosity Model parameters (', ViscModelFlag, ')')
#if (MobilityModelFlag == 2 and ViscModelFlag != 1):
#  print('   ***Error: MobilityModelFlag 2 only works with ViscModelFlag 1***   ');
#  sys.ext();

if (ViscModelFlag == 0):
# {{{

  eta_r = 1.0 # reference viscosity, typically = eta_3, the solvent

  eta = np.zeros(NIcomp_max);
  eta[0] = 1.0 # polymer viscosity
  eta[1] = 1.0 # non-solvent viscosity
  eta[2] = 1.0 # solvent viscosity
  eta[3] = 1.0 # additive viscosity
  etalist = eta.tolist()
  f1 = open('params_ViscModel.in', 'w')
  json.dump({'ViscModelFlag':ViscModelFlag, \
# (Linear Variable Viscosity)\n')
             'eta_r':eta_r, \
# (reference viscosity)\n')
#  for i in range(Ncomp):
             'eta':etalist},f1,indent=4)# eta_%d\n'%i)
  f1.close()

# }}}
elif (ViscModelFlag == 1):
# {{{

  eta_r = 1.0 # reference viscosity, typically the smallest (e.g. solvent)

  eta = np.zeros(NIcomp_max);
  eta[0] = 1.0 # polymer viscosity
  eta[1] = 1.0 # non-solvent viscosity
  eta[2] = 1.0 # solvent viscosity
  eta[3] = 1.0 # additive viscosity

  phi_T = np.zeros(NIcomp);
  w = np.zeros(NIcomp);

  phi_T[0] = 0.33; # threshold concentration of polymer
  w[0] = 0.005; # width of transition for polymer

  expon= np.zeros(NIcomp,dtype=int); # == 1 for component i to be exponential, else linear
  expon[0] = 1

  etalist = eta.tolist()
  phi_Tlist = phi_T.tolist()
  wlist = w.tolist()
  exponlist = expon.tolist()
  f1 = open('params_ViscModel.in', 'w')
  json.dump({'ViscModelFlag':ViscModelFlag, \
#(Exponential Variable Viscosity)\n')
             'eta_r':eta_r, \
# (reference viscosity)\n')
             'eta':etalist, \
              # eta_%d\n'%i)
             'phi_T':phi_Tlist, \
              # (threshold concentration of %d)\n'%i)
             'w':wlist, \
              # transition width for comp %d\n'%i)
             'expon':exponlist},f1,indent=4) # ==1 if comp %d is exponential (else linear)\n'%i)
  f1.close()
# }}}
elif (ViscModelFlag == 2):
# {{{

  eta_r = 1.0 # reference viscosity, typically the smallest (e.g. solvent)

  eta = np.zeros(NIcomp_max);
  eta[0] = 1.0 # polymer viscosity
  eta[1] = 1.0 # non-solvent viscosity
  eta[2] = 1.0 # solvent viscosity
  eta[3] = 1.0 # additive viscosity

  phi_T = np.zeros(NIcomp);
  w = np.zeros(NIcomp);

  phi_T[0] = 0.08; # threshold concentration of polymer
  w[0] = 0.008; # width of transition for polymer

  sigon = np.zeros(NIcomp,dtype=int); # == 1 for component i to be sigmoidal, else linear
  sigon[0] = 1

  etalist = eta.tolist()
  phi_Tlist = phi_T.tolist()
  wlist = w.tolist()
  sigonlist = sigon.tolist()
  f1 = open('params_ViscModel.in', 'w')
  json.dump({'ViscModelFlag':ViscModelFlag, \
#(Sigmoidal Variable Viscosity)\n')
             'eta_r':eta_r, \
# (reference viscosity)\n')
             'eta':etalist, \
              # eta_%d\n'%i)
             'phi_T':phi_Tlist, \
              # (threshold concentration of %d)\n'%i)
             'w':wlist, \
              # transition width for comp %d\n'%i)
             'sigon':sigonlist},f1,indent=4) # ==1 if comp %d is exponential (else linear)\n'%i)
  f1.close()

# }}}
elif (ViscModelFlag == 3):
# {{{

  eta_r = 1.0 # reference viscosity, typically the smallest (e.g. solvent)
  eta_s = 1.0 # solvent/nonsolvent viscosity, typically, eta_s = eta_r.
  gamma = 0.5 # between 0 and 1 is a good operating point
  phiDagger = 0.4  #glass transition concentration.
  etaCapFlag = 1 #switch for etaCap
  etaCap = 1.e+300 #cap for viscosity, i.e., computational infinity. 

  f1 = open('params_ViscModel.in', 'w')
  json.dump({'ViscModelFlag':ViscModelFlag, \
# (VFTH Viscosity Model)\n')
  'etaCapFlag':etaCapFlag, \
# (VFTH Viscosity Model Cap Flag)\n')
  'eta_r':eta_r, \
# (reference viscosity)\n')
  'eta_s':eta_s, \
# (solvent/nonsolvent viscosity)\n')
  'gamma':gamma, \
# (exponent coefficient)\n')
  'phiDagger':phiDagger, \
# (glass transition concentration)\n')
  'etaCap':etaCap},f1,indent=4) # (viscosity cap [numerical infinity])\n')
  f1.close()

# }}}

# ----------------------------------------------
#  Time Integration Parameters 
# ----------------------------------------------
# 0: Model B
# 1: Model H
# 2: Model H with deforming grid
# 3: Model B using the del_mu formulation 
# 4: Model H using the del_mu formulation (body force included)

TimeIntFlag = 1
print('  - writing Time Integration parameters (', TimeIntFlag, ')' )

# parameters
# {{{
output_interval_flag = 0; # 0: output on linear interval, 1: output on a logarithmic interval
t_max = 2e5 # final time
Nt_disp = 501  # number of timepoints to write to file

var_dt_flag = 0 # 0: No variable timestep, 1: variable timestep
dt0 = 1.e-6 # initial timestep (remains const if no step-doubling)
max_dt = 5.e-1 # be careful that this is smaller than dt_disp (below)
min_dt = 1.e-9 # this should be bigger than 1e-12

CapPhi0Flag = 0 # if flag=1, artifically stops cross-over above glassMax
phi0_max = 0.400000*0.99 #1 # maximum value of phi[0] allowed in simulation

phi_err_tol_flag = 0 # adaptive phi_err_tol flag (== 0 if constant, == 1 if adaptive)
phi_err_tol_max = 1.e-6 # maximum (or constant) tolerance for truncation error in step-doubling (only one that matters if phi_err_tol_flag == 0)
phi_err_tol_min = 1.e-9 # minimum tolerance for truncation error near glassMax
phi_err_tol_w = 3.3333e-3 # width of phi0 crossover between max and min phi_err_tol
phi_err_tol_mid = 0.02 # delta of shift in phi0 from (glassmax-phi0max) with tol = 0.5*(tol_max-tol_min)+tol_min

N_step_disp = 1 # number of timesteps between writes to job.log
Nt_max = int(1e7) # maximum number of timesteps
phi_iter_max = 100; # max number of iterations when phi doesn't converge

# ** velocity only **
v_err_tol = 1.e-6 # tolerance for error in velocities (model H only)
BodyForce_flag = 0 #0 = none, 1 = from file (for TimeInt = 1 or 4 only) 2 = component force. Only TimeInt = 1
bf_phi_T = phi_T[0] #Bodyforce phi threshold for component body force
bf_w = w[0] #Bodyforce width for component body force
bf_scale = [0,1.0,0] #scale of the bodyforce for compnent body force. [x,y,z]
bf_comp = 0 #component for the bf to act upon. Only works with bf flag 2
if bf_comp > NIcomp or bf_comp < 0:
  print('bf component must be a valid compoent')
  sys.exit()

G = np.zeros((dim, dim)) # shear rate for deforming grid only (TimeInt = 2)
if (dim > 1): G[1,0] = 2 # script breaks on this line when dim=1 
advection_flag = 0 # advection formulation flag, 1: v\cdot\grad\phi (Tree 2017); 0: \div(phi*v)

# Calculate dt_disp. write to file every dt_disp
if ( output_interval_flag == 0 ):
  dt_disp = t_max/float(Nt_disp-1)
else:
  dt_disp = np.log10(t_max)/float(Nt_disp-1)

if (output_interval_flag == 0 and max_dt > dt_disp): 
  max_dt = dt_disp;
  print("Warning: max_dt was larger than dt_disp and has been changed.")
  print("max_dt = ", max_dt)

if ((dt_disp-np.round(dt_disp/dt0)*dt0) > 1e-10):
  dt0 = dt_disp/np.ceil(dt_disp/dt0);
  print("Warning: dt0 is not an integer multiple of dt_disp and has been changed.")
  print("dt0 = ", dt0)

# check if variable viscosity or not (model H only)
var_visc_flag = 0
if (ViscModelFlag == 0 or ViscModelFlag == 1 or ViscModelFlag == 2):
  for i in range(NIcomp):
    if ( abs(eta[i]-eta[Ncomp-1])/eta_r > 1e-12 ):
      var_visc_flag = 1
      break
elif (ViscModelFlag == 3):
  var_visc_flag = 1
# }}}

# thermal fluctuations
enablefluctuations =1;
reductionFactor = 1.;

# i/o
if (TimeIntFlag == 0 or TimeIntFlag == 3):
# {{{

  f1 = open('params_TimeInt.in', 'w')
  json.dump({'TimeIntFlag':TimeIntFlag, \
  # (Model B)\n')
  't_max':t_max, \
  'dt0':dt0, \
  'max_dt':max_dt, \
  'min_dt':min_dt, \
  'dt_disp':dt_disp, \
  'output_interval_flag':output_interval_flag, \
  'Nt_max':Nt_max, \
  'N_step_disp':N_step_disp, \
  'phi_iter_max':phi_iter_max, \
  'var_dt_flag':var_dt_flag, \
  'CapPhi0Flag':CapPhi0Flag, \
  'phi0_max':phi0_max, \
  'phi_err_tol_flag':phi_err_tol_flag, \
  'phi_err_tol_max':phi_err_tol_max, \
  'phi_err_tol_min':phi_err_tol_min, \
  'phi_err_tol_w':phi_err_tol_w, \
  'phi_err_tol_mid':phi_err_tol_mid, \
  'enablefluctuations':enablefluctuations, \
  'reductionFactor':reductionFactor},f1,indent=4) 
  f1.close()

# }}}
if (TimeIntFlag == 1):
# {{{

  # i/o
  f1 = open('params_TimeInt.in', 'w')
  json.dump({'TimeIntFlag':TimeIntFlag, \
  't_max':t_max, \
  'dt0':dt0, \
  'max_dt':max_dt, \
  'min_dt':min_dt, \
  'dt_disp':dt_disp, \
  'phi_err_tol_max':phi_err_tol_max, \
  'v_err_tol':v_err_tol, \
  'var_visc_flag':var_visc_flag, \
  'output_interval_flag':output_interval_flag, \
  'Nt_max':Nt_max, \
  'N_step_disp':N_step_disp, \
  'phi_iter_max':phi_iter_max, \
  'var_dt_flag':var_dt_flag, \
  'CapPhi0Flag':CapPhi0Flag, \
  'phi0_max':phi0_max, \
  'BodyForce_flag':BodyForce_flag, \
  'bf_phiT':bf_phi_T, \
  'bf_w':bf_w, \
  'bf_scale':bf_scale, \
  'bf_comp':bf_comp, \
  'enablefluctuations':enablefluctuations, \
  'reductionFactor':reductionFactor, \
  'advection_flag':advection_flag},f1,indent=4) 
  f1.close()

# }}}
if (TimeIntFlag == 4):
# {{{

  # i/o
  f1 = open('params_TimeInt.in', 'w')
  json.dump({'TimeIntFlag':TimeIntFlag, \
  't_max':t_max, \
  'dt0':dt0, \
  'max_dt':max_dt, \
  'min_dt':min_dt, \
  'dt_disp':dt_disp, \
  'phi_err_tol_max':phi_err_tol_max, \
  'v_err_tol':v_err_tol, \
  'var_visc_flag':var_visc_flag, \
  'output_interval_flag':output_interval_flag, \
  'Nt_max':Nt_max, \
  'N_step_disp':N_step_disp, \
  'phi_iter_max':phi_iter_max, \
  'var_dt_flag':var_dt_flag, \
  'CapPhi0Flag':CapPhi0Flag, \
  'phi0_max':phi0_max, \
  'BodyForce_flag':BodyForce_flag, \
  'enablefluctuations':enablefluctuations, \
  'reductionFactor':reductionFactor, \
  'advection_flag':advection_flag},f1,indent=4) 
  f1.close()

# }}}
if (TimeIntFlag == 2):
# {{{
  Glist = G.tolist()
  # i/o
  f1 = open('params_TimeInt.in', 'w')
  json.dump({'TimeIntFlag':TimeIntFlag, \
  't_max':t_max, \
  'dt0':dt0, \
  'max_dt':max_dt, \
  'min_dt':min_dt, \
  'dt_disp':dt_disp, \
  'phi_err_tol_max':phi_err_tol_max, \
  'v_err_tol':v_err_tol, \
  'var_visc_flag':var_visc_flag, \
  'output_interval_flag':output_interval_flag, \
  'Nt_max':Nt_max, \
  'N_step_disp':N_step_disp, \
  'phi_iter_max':phi_iter_max, \
  'var_dt_flag':var_dt_flag, \
  'CapPhi0Flag':CapPhi0Flag, \
  'phi0_max':phi0_max, \
  'G':Glist},f1,indent=4)# (Velocity gradient tensor element G %d %d)\n'%(i,j,i,j))
  f1.close()

# }}}

# ----------------------------------------------
#  Reaction Model Parameters
# ----------------------------------------------
# 0: No Reactions
# 1: Binary Reaction (Reactive Blending)

ReactionsModelFlag=0

print('  - writing Reaction Model parameters (', ReactionsModelFlag, ')')

if (ReactionsModelFlag==0):
#{{{
  f1 = open('params_Reactions.in', 'w')
  json.dump({'ReactionsModelFlag':ReactionsModelFlag},f1,indent = 4) # (No Reactions)\n')
  f1.close()
#}}}
elif (ReactionsModelFlag==1):
#{{{
  n = int(N[0]/N[1]);
  kf = 1
  kb = 0/(Nr*1.5);

  f1 = open('params_Reactions.in', 'w')
  json.dump({'ReactionsModelFlag':ReactionsModelFlag, \
             'kf':kf, \
             'kb':kb, \
             'Nr':Nr, \
             'N':Nlist, \
             'n':n},f1,indent=4)
  f1.close()
#}}}


# ----------------------------------------------
#  Flags and # components to input file
# ----------------------------------------------

# {{{
f1 = open('params_Keys.in', 'w')
json.dump({'Ncomp':Ncomp, \
          'NIcomp':NIcomp, \
          'GridFlag':GridFlag, \
          'EnergyModelFlag':EnergyModelFlag, \
          'MobilityModelFlag':MobilityModelFlag, \
          'ViscModelFlag':ViscModelFlag, \
          'TimeIntFlag':TimeIntFlag, \
          'ReactionsModelFlag':ReactionsModelFlag},f1,indent=4)
f1.close()
# }}}

