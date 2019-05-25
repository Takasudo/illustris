from math import *
import numpy as np
import matplotlib.pyplot as plt
import illustris_python as il

basePath="/Users/taka/Dropbox/works/proto/TNG300-3/output"

flag_gal = 0 # [0: from snapshot, 1: from npz file]
flag_gas = 0 # [0: from snapshot, 1: from npz file]

# (1) load subhalo data

if (flag_gal == 0): # from snapshot
    tree = il.sublink.loadTree(basePath,99,id=0,fields=['SubhaloPos','SnapNum','SubhaloSFR','SubhaloVmaxRad'],onlyMPB=False)
    #--- select subhalos with tree['SnapNum']==33 (z=2)"
    subhalo_pos = []
    subhalo_sfr = []
    subhalo_size = []
    for n in range(len(tree['SnapNum'])):
        snapnum = tree['SnapNum'][n]
        if(snapnum==33):
            subhalo_pos.append(tree['SubhaloPos'][n])
            subhalo_sfr.append(tree['SubhaloSFR'][n])
            subhalo_size.append(tree['SubhaloVmaxRad'][n])
    del tree
    np.savez('subhalo',
             subhalo_pos = subhalo_pos,
             subhalo_sfr = subhalo_sfr,
             subhalo_size = subhalo_size)    
else:
    npz = np.load("subhalo.npz")
    subhalo_pos = npz['subhalo_pos']
    subhalo_sfr = npz['subhalo_sfr']
    subhalo_size = npz['subhalo_size']

# (2) set coordinates

center = subhalo_pos[0]
dx,dy,dz = 15000,15000,15000
xmin,xmax =  center[0]-dx,center[0]+dx
ymin,ymax =  center[1]-dy,center[1]+dy
zmin,zmax =  center[2]-dz,center[2]+dz

# (3) load gas data

if (flag_gas == 0): # from snapshot 
    gas_prop = il.snapshot.loadSubset(basePath,33,'gas',['Coordinates','Density'])
    # -- select gas particles in region of interest
    gas_x = gas_prop['Coordinates'][:,0][(zmin<gas_prop['Coordinates'][:,2])
                                         &(gas_prop['Coordinates'][:,2]<zmax)
                                         &(ymin<gas_prop['Coordinates'][:,1])
                                         &(gas_prop['Coordinates'][:,1]<ymax)
                                         &(xmin<gas_prop['Coordinates'][:,0])
                                         &(gas_prop['Coordinates'][:,0]<xmax)]
    gas_y = gas_prop['Coordinates'][:,0][(zmin<gas_prop['Coordinates'][:,2])
                                         &(gas_prop['Coordinates'][:,2]<zmax)
                                         &(ymin<gas_prop['Coordinates'][:,1])
                                         &(gas_prop['Coordinates'][:,1]<ymax)
                                         &(xmin<gas_prop['Coordinates'][:,0])
                                         &(gas_prop['Coordinates'][:,0]<xmax)]
    gas_density = gas_prop['Density'][(zmin<gas_prop['Coordinates'][:,2])
                                      &(gas_prop['Coordinates'][:,2]<zmax)
                                      &(ymin<gas_prop['Coordinates'][:,1])
                                      &(gas_prop['Coordinates'][:,1]<ymax)
                                      &(xmin<gas_prop['Coordinates'][:,0])
                                      &(gas_prop['Coordinates'][:,0]<xmax)]
    del gas_prop
    np.savez('gas',
             gas_x = gas_x,
             gas_y = gas_y,
             gas_density = gas_density)
else:
    npz = np.load("gas.npz")
    gas_x = npz['gas_x']
    gas_y = npz['gas_y']
    gas_density = npz['gas_density']

# (4) plot

import matplotlib.patches as patches
import matplotlib as mpl

# (4a) plot gas density map (numpy histogram2d)

# define the bin edges

xedges = np.linspace(xmin,xmax,20)
yedges = np.linspace(ymin,ymax,20)

# create a histogram H

H, xedges, yedges = np.histogram2d(gas_x, gas_y, weights=gas_density,bins=(xedges,yedges),normed=True)
H = H.T

# use pcolormesh

X,Y = np.meshgrid(xedges, yedges)
plt.figure()
plt.pcolormesh(X, Y, H)

# (4b) plot galaxies

pos_x = [row[0] for row in subhalo_pos]
pos_y = [row[1] for row in subhalo_pos]
plt.scatter(pos_x,pos_y,color='white',marker='.')

# (4c) save

plt.xlabel('x [ckpc/h]')
plt.ylabel('y [ckpc/h]')
plt.savefig("fig_gas.pdf")
