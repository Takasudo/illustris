# create_npz_gas.py
# + create npz file for gas particles inside a protocluster (bh_snap...npz)
# + need to run create_npz_subhalo.py in advance   

from math import *
import numpy as np
import illustris_python as il

# basic data

snap = 33
redshift = 2.
z0id = 0
calclabel = 'test'

basePath="/Users/takasudo/illustris/TNG300-3/output"

flag_bh = 0 # [0: from snapshot, 1: from npz file]

# (1) load subhalo data

npz = np.load('npzdata_protocluster_subhalo_snap'+str(snap)+"_z0id"+str(z0id)+"_"+str(calclabel)+".npz")
center = npz['CoM']
xmin = npz['xmin']
ymin = npz['ymin']
zmin = npz['zmin']
xmax = npz['xmax']
ymax = npz['ymax']
zmax = npz['zmax']

# (2) load and save all BH particle data within a protocluster region  

if (flag_bh == 0): # load from snapshot and save
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
else: # load from npz file
    npz = np.load("bh_snap"+str(snap)+"_z0id"+str(z0id)+"_"+str(calclabel)+".npz")
    bh_x = npz['bh_x']
    bh_y = npz['bh_y']
    bh_z = npz['bh_z']
    bh_Mdot = npz['bh_Mdot']

print("Load BH particles (#="+str(len(bh_Mdot))+")")
