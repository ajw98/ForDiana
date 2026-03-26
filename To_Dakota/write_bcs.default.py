#!/usr/bin/env python3

import numpy as np
import sys
import os
import json

if (len(sys.argv) == 2):
  runid = int(sys.argv[1]); # should be an integer
else:
  runid = 0 # if nothing supplied, run it for zero

print('Writing boundary conditions for runid: ', runid )

# Number of BCs for each variable at each boundary
# * 2nd order gradients: NBCs = 1
# * 4th order gradients: NBCs = 2
# TODO: This should be tied to the Energy model at
#       some point, but for now I'm leaving it manual.
NBCs=2;

# read in some flags
# {{{
if (os.path.isfile("params_Keys.in")):
  # 0 = periodic, 1 = non-periodic
  f = open('params_Keys.in')
  keys = json.load(f)
  Ncomp = int(keys["Ncomp"])
  NIcomp = int(keys["NIcomp"])
  GridFlag = keys["GridFlag"]
else:
  print('*** Error in write_bcs.py ***')
  print('Run write_params.py first')

if (os.path.isfile("params_Grid.in")):
  g = open('params_Grid.in')
  jsongrid = json.load(g)
  Dim = int(jsongrid["dim"])
else:
  print('*** Error in write_bcs.py ***')
  print('Cannot open params_Grid.in')
  print('Run write_params.py first')

WriteVelFlag=0;
if (os.path.isfile('params_TimeInt.in')):
  fopen = open('params_TimeInt.in')
  TimeInt_params = json.load(fopen)
  TimeIntFlag = TimeInt_params["TimeIntFlag"]
  if (TimeIntFlag == 0 or TimeIntFlag == 3): # Model B
    WriteVelFlag=0;
  if (TimeIntFlag == 1 or TimeIntFlag == 2 or TimeIntFlag == 4): # Model H
    WriteVelFlag=1;

# }}}

# ----------------------------------------------
#  Boundary Conditions 
# ----------------------------------------------
# Note: Choices apply to both phi_i and 
#       velocities if using Model H
#
# 0: Periodic 
# 1: Constant BCs
# 2: Time-dependent Phi BCs
# 3: VIPS BCs 
# 4: Two walls (not yet implemented)

BCFlag = 0;

# Check Grid versus BCs
# {{{
print('  - writing Boundary conditions of type: (', BCFlag, ')')
print('  - writing velocities? (', WriteVelFlag==1, ')')
print('  - Number of BCs for each variable at each boundary: (', NBCs, ')')

if (BCFlag != 0 and GridFlag == 0):
  print('*** Error: Inconsistent Flags. BCFlag is set to a non-periodic condition,' \
        ' but the grid is set to PS. ***')
  sys.exit();

if (BCFlag == 0 and GridFlag != 0):
  print(' *** Warning: Flags may be inconsistent. BCFlag is set to periodic conditions,' \
        ' but the grid is set to hybrid PS/FD. ***')
# }}}

if (BCFlag == 0):
# {{{

  f1 = open('params_BCs_phi.in', 'w')
  json.dump({'BCFlag':BCFlag},f1,indent=4)
  f1.close()

  if (WriteVelFlag):
    f1 = open('params_BCs_vel.in', 'w')
    json.dump({'BCFlag':BCFlag},f1,indent=4)
    f1.close()

# }}}
elif (BCFlag == 1):
# {{{

  # GhostNode flag is an extra flag determining
  # the order of accuracy of the FD approximation 
  # for the second order derivative
  # - GhostNodeFlag = 0, 1st order FD terms. This is necessary when the
  #                      2nd order is the lowest order derivative
  # - GhostNodeFlag = 1, use 2nd order all other times
  GhostNodeFlag = 1

  # --- set up BCs by hand, loop will write them to file ---
  # 2 BCs at x = 0 and 2 BCs at x = Lx,

  # (1) pick the order of derivatives of each BC
  # Set the derivative order of the BCs
  # order = 0: f = val
  # order = 1: f' = val
  # order = 2: f'' = val
  # etc
  # (note 1: no two orders can be the same)
  # (note 2: the lower order condition must come first)
  # e.g. possible combos: 0101, 0202, 1212, 0303, 1313, 2323

  Order_xeq0 = np.zeros((NBCs, NIcomp),dtype=int);
  Order_xeqL = np.zeros((NBCs, NIcomp),dtype=int);

  # matrices are rows = # BC, cols = component
  Order_xeq0[0][0] = 0
  Order_xeq0[0][1] = 0
  Order_xeq0[1][0] = 1
  Order_xeq0[1][1] = 1

  Order_xeqL[0][0] = 0
  Order_xeqL[0][1] = 0
  Order_xeqL[1][0] = 1
  Order_xeqL[1][1] = 1

  # check assumptions
  if (NBCs > 1):
    for i in range(1, NBCs):
      for j in range(NIcomp):
        if ( Order_xeq0[i-1][j] >= Order_xeq0[i][j] or
             Order_xeqL[i-1][j] >= Order_xeqL[i][j] ):
          print('  *** Error: Order[%d,%d] must be greater' + \
                'than Order[%d,%d] ***'%(i,j,i-1,j))

  # (2) pick the value of the BCs of each component

  BC_xeq0 = np.zeros((NBCs, NIcomp));
  BC_xeqL = np.zeros((NBCs, NIcomp));

  # matrices are rows = # BC, cols = component
  BC_xeq0[0][0] = 0.25 # 0.0 #
  BC_xeq0[0][1] = 0.6  # 0.0 #
  BC_xeq0[1][0] = 0.0  # 0.05
  BC_xeq0[1][1] = 0.0  #-0.02

  BC_xeqL[0][0] = 0.25 # 0.0 #
  BC_xeqL[0][1] = 0.6  # 0.0 #
  BC_xeqL[1][0] = 0.0  #-0.05
  BC_xeqL[1][1] = 0.0  # 0.02
  Order_xeq0list = Order_xeq0.tolist()
  Order_xeqLlist = Order_xeqL.tolist()
  BC_xeq0list = BC_xeq0.tolist()
  BC_xeqLlist = BC_xeqL.tolist()    
  f1 = open('params_BCs_phi.in', 'w')
  json.dump({'BCFlag':BCFlag, \
             'NBCs':NBCs, \
             'GhostNodeFlag':GhostNodeFlag, \
           # GhostNodeFlag (0 = 1st order, 1 = 2nd order)\n')
             'Order_xeq0':Order_xeq0list, \
           # Derivative Order @ x=0 for bc=%d, comp=%d\n'%(i,j))
             'Order_xeqL':Order_xeqLlist, \
           # Derivative Order @ x=L for bc=%d, comp=%d\n'%(i,j))
             'BC_xeq0':BC_xeq0list, \
           # BC @ x=0 for bc=%d, comp=%d\n'%(i,j))
             'BC_xeqL':BC_xeqLlist},f1,indent=4)
           # BC @ x=L for bc=%d, comp=%d\n'%(i,j))
  f1.close()

  if (WriteVelFlag):

    # (3) set the velocity BCs if applicable

    vel_0 = np.zeros(Dim);
    vel_L = np.zeros(Dim);

    # vel = [vx, vy, vz]
    vel_0[0] = 0 # vx
    if (Dim > 1): vel_0[1] = 0.1 # vy
    if (Dim > 2): vel_0[2] = 0 # vz

    vel_L[0] = 0 # vx
    if (Dim > 1): vel_L[1] = 0 # vy
    if (Dim > 2): vel_L[2] = 0 # vz
    vel_0list = vel_0.tolist()
    vel_Llist = vel_L.tolist()
    f1 = open('params_BCs_vel.in', 'w')
    json.dump({'BCFlag':BCFlag, \
             'vel_0':vel_0list, \
             'vel_L':vel_Llist},f1,indent=4)
    f1.close()

# }}}
elif (BCFlag == 2):
# {{{

  # GhostNode flag is an extra flag determining
  # the order of accuracy of the FD approximation 
  # for the second order derivative
  # - GhostNodeFlag = 0, 1st order FD terms. This is necessary when the
  #                      2nd order is the lowest order derivative
  # - GhostNodeFlag = 1, use 2nd order all other times
  GhostNodeFlag = 1

  # --- set up BCs by hand, loop will write them to file ---
  # 2 BCs at x = 0 and 2 BCs at x = Lx,

  # (1) pick the order of derivatives of each BC
  # Set the derivative order of the BCs
  # order = 0: f = val
  # order = 1: f' = val
  # order = 2: f'' = val
  # etc
  # (note 1: no two orders can be the same)
  # (note 2: the lower order condition must come first)
  # e.g. possible combos: 0101, 0202, 1212, 0303, 1313, 2323

  Order_xeq0 = np.zeros((NBCs, NIcomp),dtype=int);
  Order_xeqL = np.zeros((NBCs, NIcomp),dtype=int);

  # Order[BC number, component]
  # rows = BC number, cols = component
  # X = 0
  Order_xeq0[0][0] = 0
  Order_xeq0[0][1] = 0
  Order_xeq0[1][0] = 1
  Order_xeq0[1][1] = 1

  # X = L
  Order_xeqL[0][0] = 0
  Order_xeqL[0][1] = 0
  Order_xeqL[1][0] = 1
  Order_xeqL[1][1] = 1

  # check assumptions
  if (NBCs > 1):
    for i in range(1, NBCs):
      for j in range(NIcomp):
        if ( Order_xeq0[i-1][j] >= Order_xeq0[i][j] or
             Order_xeqL[i-1][j] >= Order_xeqL[i][j] ):
          print('  *** Error: Order[%d,%d] must be greater' + \
                'than Order[%d,%d] ***'%(i,j,i-1,j))

  # (2) Time series for the BCs of each component

  # Read in or set BC time data
  Nt = 51
  t = np.linspace(0, 1e2, Nt)
  BC_xeq0 = np.zeros((NBCs, NIcomp, Nt)); # BC @ x=0[BC number, component]
  BC_xeqL = np.zeros((NBCs, NIcomp, Nt)); # BC @ x=L[BC number, component]

  # X = 0
  BC_xeq0[0][0] = 0.25 + 0*t  # BC 0, component 0
  BC_xeq0[0][1] = 0.6 + 0*t   # BC 0, component 1
  BC_xeq0[1][0] = 0*t         # BC 1, component 0
  BC_xeq0[1][1] = 0*t         # BC 1, component 1

  # X = L
  BC_xeqL[0][0] = 0.20*(1. - np.exp(-t/1.5e1)) + 0.05
  BC_xeqL[0][1] = 0.55*(1. - np.exp(-t/5.e0)) + 0.05
  BC_xeqL[1][0] = 0*t
  BC_xeqL[1][1] = 0*t
  Order_xeq0list = Order_xeq0.tolist()
  Order_xeqLlist = Order_xeqL.tolist()
  # params file
  f1 = open('params_BCs_phi.in', 'w')
  json.dump({'BCFlag':BCFlag, \
             'NBCs':NBCs, \
             'GhostNodeFlag':GhostNodeFlag, \
           # GhostNodeFlag (0 = 1st order, 1 = 2nd order)\n')
             'Order_xeq0':Order_xeq0list, \
           # Derivative Order @ x=0 for bc=%d, comp=%d\n'%(i,j))
             'Order_xeqL':Order_xeqLlist, \
           # Derivative Order @ x=L for bc=%d, comp=%d\n'%(i,j))
             'Nt':Nt},f1,indent=4) # no of time points in BC file\n')
  f1.close()

  # time dependent BCs
  f1 = open('params_BCs_phi_time.in', 'wb')

  f1.write(b'# t (times)\n')
  np.savetxt(f1, t, fmt='%+22.15e')
  f1.write(b'\n')

  for i in range(NBCs):
    for j in range(NIcomp):
      f1.write(b'# BC = %d, Comp = %d, x = 0\n'%(i, j))
      np.savetxt(f1, BC_xeq0[i, j], fmt="%+22.15e")
      f1.write(b'\n')

  for i in range(NBCs):
    for j in range(NIcomp):
      f1.write(b'# BC = %d, Comp = %d, x = L\n'%(i, j))
      np.savetxt(f1, BC_xeqL[i, j], fmt='%+22.15e')
      f1.write(b'\n')

  f1.close()

  if (WriteVelFlag):

    # Read in or set data
    vel_0 = np.zeros((Dim, Nt));
    vel_L = np.zeros((Dim, Nt));

    # vel = [vx, vy, vz]
    vel_0[0] = 0. + 0.*t # vx
    if (Dim > 1): vel_0[1] = 0.+5.e-4*t # vy
    if (Dim > 2): vel_0[2] = 0.+0*t # vz

    vel_L[0] = 0.0 + 0.*t # vx
    if (Dim > 1): vel_L[1] = 0.0+0*t # vy
    if (Dim > 2): vel_L[2] = 0.0+0*t # vz

    # params file
    f1 = open('params_BCs_vel.in', 'w')
    json.dump({'BCFlag':BCFlag, \
               'Nt':Nt},f1,indent=4)
    f1.close()

    # time dependent BCs
    f1 = open('params_BCs_vel_time.in', 'wb')

    f1.write(b'# t (times)\n')
    np.savetxt(f1, t, fmt='%+22.15e')
    f1.write(b'\n')

    for i in range(Dim):
      if(i==0):
        f1.write(b'#vx at x=0\n')
        np.savetxt(f1, vel_0[0], fmt='%+22.15e')
        f1.write(b'\n')
      elif(i==1):
        f1.write(b'#vy at x=0\n')
        np.savetxt(f1, vel_0[1], fmt='%+22.15e')
        f1.write(b'\n')
      elif(i==2):
        f1.write(b'#vz at x=0\n')
        np.savetxt(f1, vel_0[2], fmt='%+22.15e')
        f1.write(b'\n')

    for i in range(Dim):
      if(i==0):
        f1.write(b'#vx at x=L\n')
        np.savetxt(f1, vel_L[0], fmt='%+22.15e')
        f1.write(b'\n')
      elif(i==1):
        f1.write(b'#vy at x=L\n')
        np.savetxt(f1, vel_L[1], fmt='%+22.15e')
        f1.write(b'\n')
      elif(i==2):
        f1.write(b'#vz at x=L\n')
        np.savetxt(f1, vel_L[2], fmt='%+22.15e')
        f1.write(b'\n')

    f1.close()

# }}}
elif (BCFlag == 3):
# {{{
  #### VIPS Boundary Conditions ####

  GhostNodeFlag = 1 #VIPS always requires this to be 1.
  ddxAcc = 2 #choice: 2, 4, 6
  #accuracy order of dphi/dx at x = Lx when Order_xeqL[0][1] = 0
  #ddxAcc does not matter when Order_xeqL[0][1] = 1

  # --- set up BCs by hand, loop will write them to file ---
  # 2 BCs at x = 0 and 2 BCs at x = Lx,

  # (1) pick the order of derivatives of each BC
  # Set the derivative order of the BCs
  # order = 0: f = val
  # order = 1: f' = val
  # order = 2: f'' = val
  # etc
  # (note 1: no two orders can be the same)
  # (note 2: the lower order condition must come first)
  # e.g. possible combos: 0101, 0202, 1212, 0303, 1313, 2323

  Order_xeq0 = np.zeros((NBCs, NIcomp),dtype=int);
  Order_xeqL = np.zeros((NBCs, NIcomp),dtype=int);

  # matrices are rows = # BC, cols = component
  
  # Bottom of the film
  Order_xeq0[0][0] = 1
  Order_xeq0[0][1] = 1
  Order_xeq0[1][0] = 3
  Order_xeq0[1][1] = 3
  #x=0 notes:
  #First and third order derivatives are usually
  #set to zero because it's the bottom of the film.
  #You are free to change these BCs though.

  # Top of the film
  Order_xeqL[0][0] = 1
  Order_xeqL[0][1] = 0 #Set to either 0 or 1.
  Order_xeqL[1][0] = 3
  Order_xeqL[1][1] = 2 #VIPS ignores. 
  #x=Lx notes:
  #First BC for nonsolvent must be either zero or
  #first order. Second nonsolvent BC is ignored.
 
  # check assumptions
  if (NBCs > 1):
    for i in range(1, NBCs):
      for j in range(NIcomp):
        if ( Order_xeq0[i-1][j] >= Order_xeq0[i][j] or
             Order_xeqL[i-1][j] >= Order_xeqL[i][j] ):
          print('  *** Error: Order[%d,%d] must be greater' + \
                'than Order[%d,%d] ***'%(i,j,i-1,j))
 
  if (ddxAcc != 2 and ddxAcc != 4 and ddxAcc != 6 and Order_xeqL[0][1] == 0):
    print ("Pick a valid accuracy order");
    sys.exit();

  # (2) pick the value of the BCs of each component

  BC_xeq0 = np.zeros((NBCs, NIcomp));
  BC_xeqL = np.zeros((NBCs, NIcomp));

  # matrices are rows = # BC, cols = component
  
  # Bottom of the film 
  BC_xeq0[0][0] = 0.0 
  BC_xeq0[0][1] = 0.0
  BC_xeq0[1][0] = 0.0 
  BC_xeq0[1][1] = 0.0
  #x=0 notes:
  #First and third order derivatives are
  #usually set to zero because it's the bottom of the film.
  #You are free to change these BCs though.
  
  # Top of the film
  BC_xeqL[0][0] = 0.0 #set to zero 
  BC_xeqL[0][1] = 0.25 #controls nonsolvent entry
  BC_xeqL[1][0] = 0.0 #set to zero 
  BC_xeqL[1][1] = 0.0 #VIPS ignores
  #x=Lx notes:
  #Polymer 1st and 3rd order derivatives must be 
  #set to zero. Second nonsolvent BC is ignored.
    
  Order_xeq0list = Order_xeq0.tolist()
  Order_xeqLlist = Order_xeqL.tolist()
  BC_xeq0list = BC_xeq0.tolist()
  BC_xeqLlist = BC_xeqL.tolist()    
  f1 = open('params_BCs_phi.in', 'w')
  json.dump({'BCFlag':BCFlag, \
             'NBCs':NBCs, \
             'GhostNodeFlag':GhostNodeFlag, \
           # GhostNodeFlag (0 = 1st order, 1 = 2nd order)\n')
             'Order_xeq0':Order_xeq0list, \
           # Derivative Order @ x=0 for bc=%d, comp=%d\n'%(i,j))
             'Order_xeqL':Order_xeqLlist, \
           # Derivative Order @ x=L for bc=%d, comp=%d\n'%(i,j))
             'BC_xeq0':BC_xeq0list, \
           # BC @ x=0 for bc=%d, comp=%d\n'%(i,j))
             'BC_xeqL':BC_xeqLlist, \
           # BC @ x=L for bc=%d, comp=%d\n'%(i,j))
             'ddxAcc':ddxAcc},f1,indent=4)
  f1.close()

  if (WriteVelFlag):
     print('  *** Error: VIPS and Model H not supported.')
# }}}    

