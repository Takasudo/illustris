# --- create_npz_bh.py :      
# (1) read protocluster subhalo data from npz file
# (2) read BH data from illustris snapshot within radius d_N
# (3) save all data to "npzdata_protocluster_bh" file 

from math import *
import numpy as np
import illustris_python as il
import sys
args = sys.argv

snap = int(args[1])
redshift = float(args[2])
z0id = int(args[3])
calclabel = args[4]

if snap==33:
    basePath="/Users/taka/Dropbox/works/proto/TNG300-3/output"
else:
    basePath="/Users/taka/Desktop/TNG300-3/output"

flag_bh = 0 # [0: from snapshot, 1: from npz file]  

# Cosmology

hubble = 0.7

# (1) load subhalo map data

npz = np.load('npzdata_protocluster_subhalo_snap'+str(snap)+"_z0id"+str(z0id)+"_"+str(calclabel)+".npz")
center = npz['CoM']
xmin = npz['xmin']
ymin = npz['ymin']
zmin = npz['zmin']
xmax = npz['xmax']
ymax = npz['ymax']
zmax = npz['zmax']
d_20 = npz['d_20']
d_50 = npz['d_50']
d_80 = npz['d_80']
d_100 = npz['d_100']

# (2) bh particles

# (2a) load all bh data within a box of (d_100)^3

if (flag_bh == 0): # from snapshot
    bh_prop = il.snapshot.loadSubset(basePath,snap,'bh',['Coordinates','BH_Mdot','BH_Mass'])
    print("Load all BH data")
    bh_x = bh_prop['Coordinates'][:,0][(zmin<bh_prop['Coordinates'][:,2])
                                       &(bh_prop['Coordinates'][:,2]<zmax)
                                       &(ymin<bh_prop['Coordinates'][:,1])
                                       &(bh_prop['Coordinates'][:,1]<ymax)
                                       &(xmin<bh_prop['Coordinates'][:,0])
                                       &(bh_prop['Coordinates'][:,0]<xmax)]
    print("set bh_x")
    bh_y = bh_prop['Coordinates'][:,1][(zmin<bh_prop['Coordinates'][:,2])
                                       &(bh_prop['Coordinates'][:,2]<zmax)
                                       &(ymin<bh_prop['Coordinates'][:,1])
                                       &(bh_prop['Coordinates'][:,1]<ymax)
                                       &(xmin<bh_prop['Coordinates'][:,0])
                                       &(bh_prop['Coordinates'][:,0]<xmax)]
    print("set bh_y")
    bh_z = bh_prop['Coordinates'][:,2][(zmin<bh_prop['Coordinates'][:,2])
                                       &(bh_prop['Coordinates'][:,2]<zmax)
                                       &(ymin<bh_prop['Coordinates'][:,1])
                                       &(bh_prop['Coordinates'][:,1]<ymax)
                                       &(xmin<bh_prop['Coordinates'][:,0])
                                       &(bh_prop['Coordinates'][:,0]<xmax)]
    print("set bh_z")
    bh_Mdot = bh_prop['BH_Mdot'][(zmin<bh_prop['Coordinates'][:,2])
                                 &(bh_prop['Coordinates'][:,2]<zmax)
                                 &(ymin<bh_prop['Coordinates'][:,1])
                                 &(bh_prop['Coordinates'][:,1]<ymax)
                                 &(xmin<bh_prop['Coordinates'][:,0])
                                 &(bh_prop['Coordinates'][:,0]<xmax)]
    print("set bh_Mdot")
    del bh_prop
    np.savez('bh_snap'+str(snap)+"_z0id"+str(z0id)+"_"+str(calclabel),
             bh_x = bh_x,
             bh_y = bh_y,
             bh_z = bh_z,
             bh_Mdot = bh_Mdot
             )
else: # from npz file
    npz = np.load("bh_snap"+str(snap)+"_z0id"+str(z0id)+"_"+str(calclabel)+".npz")
    bh_x = npz['bh_x']
    bh_y = npz['bh_y']
    bh_z = npz['bh_z']
    bh_Mdot = npz['bh_Mdot']

print("Load gas particles (#="+str(len(bh_Mdot))+")")

# (2b) select gas particles within spheres of d50, d80, d100

bh_Mdot_d20 = bh_Mdot[((bh_x-center[0])**2.0 + (bh_y-center[1])**2.0 + (bh_z-center[2])**2.0)**0.5 < d_20]
bh_Mdot_d50 = bh_Mdot[((bh_x-center[0])**2.0 + (bh_y-center[1])**2.0 + (bh_z-center[2])**2.0)**0.5 < d_50]
bh_Mdot_d80 = bh_Mdot[((bh_x-center[0])**2.0 + (bh_y-center[1])**2.0 + (bh_z-center[2])**2.0)**0.5 < d_80]
bh_Mdot_d100 = bh_Mdot[((bh_x-center[0])**2.0 + (bh_y-center[1])**2.0 + (bh_z-center[2])**2.0)**0.5 < d_100]
del bh_Mdot

bh_x_d20 = bh_x[((bh_x-center[0])**2.0 + (bh_y-center[1])**2.0 + (bh_z-center[2])**2.0)**0.5 < d_20]
bh_x_d50 = bh_x[((bh_x-center[0])**2.0 + (bh_y-center[1])**2.0 + (bh_z-center[2])**2.0)**0.5 < d_50]
bh_x_d80 = bh_x[((bh_x-center[0])**2.0 + (bh_y-center[1])**2.0 + (bh_z-center[2])**2.0)**0.5 < d_80]
bh_x_d100 = bh_x[((bh_x-center[0])**2.0 + (bh_y-center[1])**2.0 + (bh_z-center[2])**2.0)**0.5 < d_100]

bh_y_d20 = bh_y[((bh_x-center[0])**2.0 + (bh_y-center[1])**2.0 + (bh_z-center[2])**2.0)**0.5 < d_20]
bh_y_d50 = bh_y[((bh_x-center[0])**2.0 + (bh_y-center[1])**2.0 + (bh_z-center[2])**2.0)**0.5 < d_50]
bh_y_d80 = bh_y[((bh_x-center[0])**2.0 + (bh_y-center[1])**2.0 + (bh_z-center[2])**2.0)**0.5 < d_80]
bh_y_d100 = bh_y[((bh_x-center[0])**2.0 + (bh_y-center[1])**2.0 + (bh_z-center[2])**2.0)**0.5 < d_100]

bh_z_d20 = bh_z[((bh_x-center[0])**2.0 + (bh_y-center[1])**2.0 + (bh_z-center[2])**2.0)**0.5 < d_20]
bh_z_d50 = bh_z[((bh_x-center[0])**2.0 + (bh_y-center[1])**2.0 + (bh_z-center[2])**2.0)**0.5 < d_50]
bh_z_d80 = bh_z[((bh_x-center[0])**2.0 + (bh_y-center[1])**2.0 + (bh_z-center[2])**2.0)**0.5 < d_80]
bh_z_d100 = bh_z[((bh_x-center[0])**2.0 + (bh_y-center[1])**2.0 + (bh_z-center[2])**2.0)**0.5 < d_100]

del bh_x
del bh_y
del bh_z

# (3) save data

np.savez("npzdata_protocluster_bh_snap"+str(snap)+"_z0id"+str(z0id)+"_"+str(calclabel),
         bh_Mdot_d20 = bh_Mdot_d20,
         bh_Mdot_d50 = bh_Mdot_d50,
         bh_Mdot_d80 = bh_Mdot_d80,
         bh_Mdot_d100 = bh_Mdot_d100)
