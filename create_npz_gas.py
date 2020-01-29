# create_npz_gas.py
# + create npz file for gas particles inside a protocluster (gas_snap...npz)
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

flag_gas = 0 # [0: from snapshot, 1: from npz file]  

# (1) load subhalo data

npz = np.load('npzdata_protocluster_subhalo_snap'+str(snap)+"_z0id"+str(z0id)+"_"+str(calclabel)+".npz")

center = npz['CoM']
xmin = npz['xmin']
ymin = npz['ymin']
zmin = npz['zmin']
xmax = npz['xmax']
ymax = npz['ymax']
zmax = npz['zmax']

# (2) load and save all gas particle data within a protocluster region

if (flag_gas == 0): # load from snapshot and save
    gas_prop = il.snapshot.loadSubset(basePath,snap,'gas',['Coordinates','Density','Masses','MagneticField','StarFormationRate','InternalEnergy','ElectronAbundance','ParticleIDs'])
    print("Load all gas data")
    gas_x = gas_prop['Coordinates'][:,0][(zmin<gas_prop['Coordinates'][:,2])
                                         &(gas_prop['Coordinates'][:,2]<zmax)
                                         &(ymin<gas_prop['Coordinates'][:,1])
                                         &(gas_prop['Coordinates'][:,1]<ymax)
                                         &(xmin<gas_prop['Coordinates'][:,0])
                                         &(gas_prop['Coordinates'][:,0]<xmax)]
    print("set gas_x")
    gas_y = gas_prop['Coordinates'][:,1][(zmin<gas_prop['Coordinates'][:,2])
                                         &(gas_prop['Coordinates'][:,2]<zmax)
                                         &(ymin<gas_prop['Coordinates'][:,1])
                                         &(gas_prop['Coordinates'][:,1]<ymax)
                                         &(xmin<gas_prop['Coordinates'][:,0])
                                         &(gas_prop['Coordinates'][:,0]<xmax)]
    print("set gas_y")
    gas_z = gas_prop['Coordinates'][:,2][(zmin<gas_prop['Coordinates'][:,2])
                                         &(gas_prop['Coordinates'][:,2]<zmax)
                                         &(ymin<gas_prop['Coordinates'][:,1])
                                         &(gas_prop['Coordinates'][:,1]<ymax)
                                         &(xmin<gas_prop['Coordinates'][:,0])
                                         &(gas_prop['Coordinates'][:,0]<xmax)]
    print("set gas_z")
    gas_density = gas_prop['Density'][(zmin<gas_prop['Coordinates'][:,2])
                                      &(gas_prop['Coordinates'][:,2]<zmax)
                                      &(ymin<gas_prop['Coordinates'][:,1])
                                      &(gas_prop['Coordinates'][:,1]<ymax)
                                      &(xmin<gas_prop['Coordinates'][:,0])
                                      &(gas_prop['Coordinates'][:,0]<xmax)]
    print("set gas_density")
    gas_mass = gas_prop['Masses'][(zmin<gas_prop['Coordinates'][:,2])
                                    &(gas_prop['Coordinates'][:,2]<zmax)
                                    &(ymin<gas_prop['Coordinates'][:,1])
                                    &(gas_prop['Coordinates'][:,1]<ymax)
                                    &(xmin<gas_prop['Coordinates'][:,0])
                                    &(gas_prop['Coordinates'][:,0]<xmax)]
    print("set gas_mass")
    gas_mag = gas_prop['MagneticField'][(zmin<gas_prop['Coordinates'][:,2])
                                     &(gas_prop['Coordinates'][:,2]<zmax)
                                     &(ymin<gas_prop['Coordinates'][:,1])
                                     &(gas_prop['Coordinates'][:,1]<ymax)
                                     &(xmin<gas_prop['Coordinates'][:,0])
                                     &(gas_prop['Coordinates'][:,0]<xmax)]
    print("set gas_mag")
    gas_sfr = gas_prop['StarFormationRate'][(zmin<gas_prop['Coordinates'][:,2])
                                     &(gas_prop['Coordinates'][:,2]<zmax)
                                     &(ymin<gas_prop['Coordinates'][:,1])
                                     &(gas_prop['Coordinates'][:,1]<ymax)
                                     &(xmin<gas_prop['Coordinates'][:,0])
                                     &(gas_prop['Coordinates'][:,0]<xmax)]
    print("set gas_sfr")
    gas_u = gas_prop['InternalEnergy'][(zmin<gas_prop['Coordinates'][:,2])
                                       &(gas_prop['Coordinates'][:,2]<zmax)
                                       &(ymin<gas_prop['Coordinates'][:,1])
                                       &(gas_prop['Coordinates'][:,1]<ymax)
                                       &(xmin<gas_prop['Coordinates'][:,0])
                                       &(gas_prop['Coordinates'][:,0]<xmax)]
    print("set gas_u")
    gas_id = gas_prop['ParticleIDs'][(zmin<gas_prop['Coordinates'][:,2])
                                     &(gas_prop['Coordinates'][:,2]<zmax)
                                     &(ymin<gas_prop['Coordinates'][:,1])
                                     &(gas_prop['Coordinates'][:,1]<ymax)
                                     &(xmin<gas_prop['Coordinates'][:,0])
                                     &(gas_prop['Coordinates'][:,0]<xmax)]
    print("set gas_id")
    gas_xe = gas_prop['ElectronAbundance'][(zmin<gas_prop['Coordinates'][:,2])
                                           &(gas_prop['Coordinates'][:,2]<zmax)
                                           &(ymin<gas_prop['Coordinates'][:,1])
                                           &(gas_prop['Coordinates'][:,1]<ymax)
                                           &(xmin<gas_prop['Coordinates'][:,0])
                                           &(gas_prop['Coordinates'][:,0]<xmax)]
    del gas_prop
    np.savez('gas_snap'+str(snap)+"_z0id"+str(z0id)+"_"+str(calclabel),
             gas_x = gas_x,
             gas_y = gas_y,
             gas_z = gas_z,
             gas_density = gas_density,
             gas_mass = gas_mass,
             gas_mag = gas_mag,
             gas_sfr = gas_sfr,
             gas_u = gas_u,
             gas_id = gas_id,
             gas_xe = gas_xe
             )
else: # load from existing npz file
    npz = np.load("gas_snap"+str(snap)+"_z0id"+str(z0id)+"_"+str(calclabel)+".npz")
    gas_x = npz['gas_x']
    gas_y = npz['gas_y']
    gas_z = npz['gas_z']
    gas_density = npz['gas_density']
    gas_mass = npz['gas_mass']
    gas_mag = npz['gas_mag']
    gas_sfr = npz['gas_sfr']
    gas_u = npz['gas_u']
    gas_id = npz['gas_id']

print("Load gas particles (#="+str(len(gas_density))+")")
