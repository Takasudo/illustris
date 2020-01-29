# create_npz_halo.py
# + create npz file for subhalos inside a protocluster (subhalo_snap...npz)
# + create npz file for propertoes of the protocluster (npzdata_protocluster_subhalo_snap...npz)

from math import *
import numpy as np
import illustris_python as il

# basic data

snap = 33
redshift = 2.
z0id = 0
calclabel = 'test'

basePath="/Users/takasudo/illustris/TNG300-3/output"

flag_gal = 0 # [0: from snapshot, 1: from npz file]

# parameter

hubble = 0.7

# (1) load all subhalos that belong to a protocluster

if (flag_gal == 0): # load from snapshot and save
    tree = il.sublink.loadTree(basePath,99,z0id,fields=['SubhaloPos','SnapNum','SubhaloMass','SubhaloSFR','SubhaloVmaxRad','SubhaloBHMdot'],onlyMPB=False)
    #--- select subhalos with tree['SnapNum']==33 (z=2)"
    subhalo_pos = []
    subhalo_mass = []
    subhalo_sfr = []
    subhalo_size = []
    subhalo_bh = []
    for n in range(len(tree['SnapNum'])):
        snapnum = tree['SnapNum'][n]
        if(snapnum==snap):
            subhalo_pos.append(tree['SubhaloPos'][n])
            subhalo_sfr.append(tree['SubhaloSFR'][n])
            subhalo_mass.append(tree['SubhaloMass'][n])
            subhalo_size.append(tree['SubhaloVmaxRad'][n])
            subhalo_bh.append(tree['SubhaloBHMdot'][n])
    del tree
    np.savez('subhalo_snap'+str(snap)+"_z0id"+str(z0id)+"_"+str(calclabel),
             subhalo_pos = subhalo_pos,
             subhalo_sfr = subhalo_sfr,
             subhalo_mass = subhalo_mass,
             subhalo_size = subhalo_size,
             subhalo_bh = subhalo_bh)    
else: # load from npz file
    npz = np.load('subhalo_snap'+str(snap)+"_z0id"+str(z0id)+"_"+str(calclabel)+".npz")
    subhalo_pos = npz['subhalo_pos']
    subhalo_sfr = npz['subhalo_sfr']
    subhalo_mass = npz['subhalo_mass']
    subhalo_size = npz['subhalo_size']
    subhalo_bh = npz['subhalo_bh']

# (2) properties of subhalos within different radii

# (2a) set np arrays

Nhalo = len(subhalo_mass)     # total number of subhalos
mass = np.array(subhalo_mass) # array of subhalo mass
pos_x = np.array([row[0] for row in subhalo_pos]) # x of subhalos
pos_y = np.array([row[1] for row in subhalo_pos]) # y of subhalos
pos_z = np.array([row[2] for row in subhalo_pos]) # z of subhalos    

# (2b) define center and radius (d20, d50, d80, d100)

center = [sum(pos_x*mass)/sum(subhalo_mass), sum(pos_y*mass)/sum(subhalo_mass), sum(pos_z*mass)/sum(subhalo_mass)] # Center of Mass for all subhalos

diff = subhalo_pos - np.array([center for n in range(Nhalo)]) # List of (subhalo position) - (center of mass)

distance_from_center = []
for d in range (len(diff)):
    dist = (diff[d][0]**2.0 + diff[d][1]**2.0 + diff[d][2]**2.0)**0.5
    distance_from_center.append(dist)
distance = np.array(distance_from_center) # List of distances between subhalo and CoM

sorted_distance = np.sort(distance)       # Sorted distances
d_20 = sorted_distance[int(Nhalo*0.2)]    # 20% distance
d_50 = sorted_distance[int(Nhalo*0.5)]    # 50% distance
d_80 = sorted_distance[int(Nhalo*0.8)]    # 80% distance
d_100 = sorted_distance[Nhalo-1]          # 100% distance

if len(subhalo_size)==1:
    d_100 = subhalo_size[0]  # re-define 100% distance

xmin, xmax = center[0] - d_100, center[0] + d_100
ymin, ymax = center[1] - d_100, center[1] + d_100
zmin, zmax = center[2] - d_100, center[2] + d_100 

# (2b) define physical quantities within d20, d50, d80, d100

subhalo_d20_pos = []
subhalo_d20_sfr = []
subhalo_d20_size = []
subhalo_d20_mass = []

subhalo_d50_pos = []
subhalo_d50_sfr = []
subhalo_d50_size = [] 
subhalo_d50_mass = []

subhalo_d80_pos = []
subhalo_d80_sfr = []
subhalo_d80_size = []
subhalo_d80_mass = []

subhalo_d100_pos = []
subhalo_d100_sfr = []
subhalo_d100_size = []
subhalo_d100_mass = []

for i in range (Nhalo):
    if (distance[i] < d_20):
        subhalo_d20_pos.append(subhalo_pos[i])
        subhalo_d20_sfr.append(subhalo_sfr[i])
        subhalo_d20_size.append(subhalo_size[i])
        subhalo_d20_mass.append(subhalo_mass[i])
    if (distance[i] < d_50):
        subhalo_d50_pos.append(subhalo_pos[i])
        subhalo_d50_sfr.append(subhalo_sfr[i])
        subhalo_d50_size.append(subhalo_size[i])
        subhalo_d50_mass.append(subhalo_mass[i])
    if (distance[i] < d_80):
        subhalo_d80_pos.append(subhalo_pos[i])
        subhalo_d80_sfr.append(subhalo_sfr[i])
        subhalo_d80_size.append(subhalo_size[i])
        subhalo_d80_mass.append(subhalo_mass[i])
    if (distance[i] < d_100):
        subhalo_d100_pos.append(subhalo_pos[i])
        subhalo_d100_sfr.append(subhalo_sfr[i])
        subhalo_d100_size.append(subhalo_size[i])
        subhalo_d100_mass.append(subhalo_mass[i])

print("Load subhalo data")

# (2c) find most massive halo

for i in range (len(mass)):
    if mass[i] == np.max(mass):
        most_massive_pos = subhalo_pos[i]
        most_massive_mass = mass[i]
        most_massive_size = subhalo_size[i]

# (2d) from comiving kpc/h to physical Mpc

d_20_pMpc = d_20*1.e-3/(1.+redshift)/hubble    # from [ckpc/h] to [Mpc]
d_50_pMpc = d_50*1.e-3/(1.+redshift)/hubble
d_80_pMpc = d_80*1.e-3/(1.+redshift)/hubble
d_100_pMpc = d_100*1.e-3/(1.+redshift)/hubble

# (3) save data

np.savez("npzdata_protocluster_subhalo_snap"+str(snap)+"_z0id"+str(z0id)+"_"+str(calclabel),
         CoM = center,
         xmin = xmin,
         ymin = ymin,
         zmin = zmin,
         xmax = xmax,
         ymax = ymax,
         zmax = zmax,
         d_20 = d_20,
         d_50 = d_50,
         d_80 = d_80,
         d_100 = d_100,
         d20_physMpc = d_20_pMpc,
         d50_physMpc = d_50_pMpc,
         d80_physMpc = d_80_pMpc,
         d100_physMpc = d_100_pMpc,
         subhalo_d20_sfr = subhalo_d20_sfr, # [M_sun/yr]
         subhalo_d50_sfr = subhalo_d50_sfr,
         subhalo_d80_sfr = subhalo_d80_sfr,
         subhalo_d100_sfr = subhalo_d100_sfr,
         subhalo_d20_mass = subhalo_d20_mass, # [10^10 M_sun/h]
         subhalo_d50_mass = subhalo_d50_mass,
         subhalo_d80_mass = subhalo_d80_mass,
         subhalo_d100_mass = subhalo_d100_mass,
         subhalo_d20_pos = subhalo_d20_pos,
         subhalo_d50_pos = subhalo_d50_pos,
         subhalo_d80_pos = subhalo_d80_pos,
         subhalo_d100_pos = subhalo_d100_pos
         )
