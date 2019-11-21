# --- create_npz_gas.py :      
# (1) read protocluster subhalo data from npz file
# (2) read gas data from illustris snapshot within radius d_N
# (3) save all data to "npzdata_protocluster_gas" file 

from math import *
import numpy as np
import illustris_python as il
import sys
args = sys.argv

snap = int(args[1])
redshift = float(args[2])
# 99 : z=0
# 33 : z=2
# 50 : z=1
# 67 : z=0.5
# 91 : z=0.1
# 98 : z=0.01

z0id = int(args[3])
# id = 0 : most massive at z=0

calclabel = args[4]

if snap==33 or snap==34:
    basePath="/Users/taka/Dropbox/works/proto/TNG300-3/output"
else:
    basePath="/Users/taka/Desktop/TNG300-3/output"

flag_gas = 0 # [0: from snapshot, 1: from npz file]  

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

# (2) gas particles

# (2a) load all gas data within a box of (d_100)^3

if (flag_gas == 0): # from snapshot
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
else: # from npz file
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

# (2b) select gas particles within spheres of d50, d80, d100

gas_density_d20 = gas_density[((gas_x-center[0])**2.0 + (gas_y-center[1])**2.0 + (gas_z-center[2])**2.0)**0.5 < d_20]
gas_density_d50 = gas_density[((gas_x-center[0])**2.0 + (gas_y-center[1])**2.0 + (gas_z-center[2])**2.0)**0.5 < d_50]
gas_density_d80 = gas_density[((gas_x-center[0])**2.0 + (gas_y-center[1])**2.0 + (gas_z-center[2])**2.0)**0.5 < d_80]
gas_density_d100 = gas_density[((gas_x-center[0])**2.0 + (gas_y-center[1])**2.0 + (gas_z-center[2])**2.0)**0.5 < d_100]
del gas_density

gas_mag_d20 = gas_mag[((gas_x-center[0])**2.0 + (gas_y-center[1])**2.0 + (gas_z-center[2])**2.0)**0.5 < d_20]
gas_mag_d50 = gas_mag[((gas_x-center[0])**2.0 + (gas_y-center[1])**2.0 + (gas_z-center[2])**2.0)**0.5 < d_50]
gas_mag_d80 = gas_mag[((gas_x-center[0])**2.0 + (gas_y-center[1])**2.0 + (gas_z-center[2])**2.0)**0.5 < d_80]
gas_mag_d100 = gas_mag[((gas_x-center[0])**2.0 + (gas_y-center[1])**2.0 + (gas_z-center[2])**2.0)**0.5 < d_100]
del gas_mag

gas_sfr_d20 = gas_sfr[((gas_x-center[0])**2.0 + (gas_y-center[1])**2.0 + (gas_z-center[2])**2.0)**0.5 < d_20]
gas_sfr_d50 = gas_sfr[((gas_x-center[0])**2.0 + (gas_y-center[1])**2.0 + (gas_z-center[2])**2.0)**0.5 < d_50]
gas_sfr_d80 = gas_sfr[((gas_x-center[0])**2.0 + (gas_y-center[1])**2.0 + (gas_z-center[2])**2.0)**0.5 < d_80]
gas_sfr_d100 = gas_sfr[((gas_x-center[0])**2.0 + (gas_y-center[1])**2.0 + (gas_z-center[2])**2.0)**0.5 < d_100]
del gas_sfr

gas_u_d20 = gas_u[((gas_x-center[0])**2.0 + (gas_y-center[1])**2.0 + (gas_z-center[2])**2.0)**0.5 < d_20]
gas_u_d50 = gas_u[((gas_x-center[0])**2.0 + (gas_y-center[1])**2.0 + (gas_z-center[2])**2.0)**0.5 < d_50]
gas_u_d80 = gas_u[((gas_x-center[0])**2.0 + (gas_y-center[1])**2.0 + (gas_z-center[2])**2.0)**0.5 < d_80]
gas_u_d100 = gas_u[((gas_x-center[0])**2.0 + (gas_y-center[1])**2.0 + (gas_z-center[2])**2.0)**0.5 < d_100]
del gas_u

gas_mass_d20 = gas_mass[((gas_x-center[0])**2.0 + (gas_y-center[1])**2.0 + (gas_z-center[2])**2.0)**0.5 < d_20]
gas_mass_d50 = gas_mass[((gas_x-center[0])**2.0 + (gas_y-center[1])**2.0 + (gas_z-center[2])**2.0)**0.5 < d_50]
gas_mass_d80 = gas_mass[((gas_x-center[0])**2.0 + (gas_y-center[1])**2.0 + (gas_z-center[2])**2.0)**0.5 < d_80]
gas_mass_d100 = gas_mass[((gas_x-center[0])**2.0 + (gas_y-center[1])**2.0 + (gas_z-center[2])**2.0)**0.5 < d_100]
del gas_mass

gas_x_d20 = gas_x[((gas_x-center[0])**2.0 + (gas_y-center[1])**2.0 + (gas_z-center[2])**2.0)**0.5 < d_20]
gas_x_d50 = gas_x[((gas_x-center[0])**2.0 + (gas_y-center[1])**2.0 + (gas_z-center[2])**2.0)**0.5 < d_50]
gas_x_d80 = gas_x[((gas_x-center[0])**2.0 + (gas_y-center[1])**2.0 + (gas_z-center[2])**2.0)**0.5 < d_80]
gas_x_d100 = gas_x[((gas_x-center[0])**2.0 + (gas_y-center[1])**2.0 + (gas_z-center[2])**2.0)**0.5 < d_100]

gas_y_d20 = gas_y[((gas_x-center[0])**2.0 + (gas_y-center[1])**2.0 + (gas_z-center[2])**2.0)**0.5 < d_20]
gas_y_d50 = gas_y[((gas_x-center[0])**2.0 + (gas_y-center[1])**2.0 + (gas_z-center[2])**2.0)**0.5 < d_50]
gas_y_d80 = gas_y[((gas_x-center[0])**2.0 + (gas_y-center[1])**2.0 + (gas_z-center[2])**2.0)**0.5 < d_80]
gas_y_d100 = gas_y[((gas_x-center[0])**2.0 + (gas_y-center[1])**2.0 + (gas_z-center[2])**2.0)**0.5 < d_100]

gas_z_d20 = gas_z[((gas_x-center[0])**2.0 + (gas_y-center[1])**2.0 + (gas_z-center[2])**2.0)**0.5 < d_20]
gas_z_d50 = gas_z[((gas_x-center[0])**2.0 + (gas_y-center[1])**2.0 + (gas_z-center[2])**2.0)**0.5 < d_50]
gas_z_d80 = gas_z[((gas_x-center[0])**2.0 + (gas_y-center[1])**2.0 + (gas_z-center[2])**2.0)**0.5 < d_80]
gas_z_d100 = gas_z[((gas_x-center[0])**2.0 + (gas_y-center[1])**2.0 + (gas_z-center[2])**2.0)**0.5 < d_100]

del gas_x
del gas_y
del gas_z

# (3) save data

np.savez("npzdata_protocluster_gas_snap"+str(snap)+"_z0id"+str(z0id)+"_"+str(calclabel),
         gas_density_d20 = gas_density_d20,
         gas_density_d50 = gas_density_d50,
         gas_density_d80 = gas_density_d80,
         gas_density_d100 = gas_density_d100,
         gas_mag_d20 = gas_mag_d20,
         gas_mag_d50 = gas_mag_d50,
         gas_mag_d80 = gas_mag_d80,
         gas_mag_d100 = gas_mag_d100,
         gas_sfr_d20 = gas_sfr_d20,
         gas_sfr_d50 = gas_sfr_d50,
         gas_sfr_d80 = gas_sfr_d80,
         gas_sfr_d100 = gas_sfr_d100,
         gas_u_d20 = gas_u_d20, # u : internal energy
         gas_u_d50 = gas_u_d50,
         gas_u_d80 = gas_u_d80,
         gas_u_d100 = gas_u_d100,
         gas_mass_d20 = gas_mass_d20,
         gas_mass_d50 = gas_mass_d50,
         gas_mass_d80 = gas_mass_d80,
         gas_mass_d100 = gas_mass_d100)
