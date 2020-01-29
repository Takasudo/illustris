# radial_info.py
# + need to run create_npz_gas.py and create_npz_bh.py in advance 

from math import *
import numpy as np

# basic data

snap = 33
redshift = 2.
z0id = 0
calclabel = 'test'

basePath="/Users/takasudo/illustris/TNG300-3/output"

flag_gal = 0 # [0: from snapshot, 1: from npz file]

# parameter / unit conversion

hubble = 0.7
sigma_pp = 3.e-26 # [cm^2]
c = 3.e10         # [cm] 

num_gas = (1e10 / hubble) * 1.989e33 / 1.66e-24       # from mass [10^10 M_sun/h] to gas num [Particle Number]
mag_gauss = hubble * (1.0 + redshift)**2.0 * 2.60e-6  # from [(h/a^2)(Unit Mass/Unit Length)^{1/2}/(Unit Time)] to [cgs Gauss]
length_kpc = 1./(1.+redshift)/hubble                  # from [ckpc/h] to [kpc]
km_to_cm = 1.e5
kpc_to_cm = 3.086e21
density_cgs = (1e10 / hubble) * 1.989e33 / (length_kpc*kpc_to_cm)**3.0 # from [(10^10 M_sun/h)/(ckpc/h)^3] to [g/cm^3]   

# (1) load npz files

npz_subhalo = np.load("npzdata_protocluster_subhalo_snap"+str(snap)+"_z0id"+str(z0id)+"_"+str(calclabel)+".npz")
npz_gas = np.load("gas_snap"+str(snap)+"_z0id"+str(z0id)+"_"+str(calclabel)+".npz")
npz_bh = np.load("bh_snap"+str(snap)+"_z0id"+str(z0id)+"_"+str(calclabel)+".npz")

# (2) load subhalo data (center and radius)

center = npz_subhalo['CoM']
d_20 = npz_subhalo['d_20']
d_50 = npz_subhalo['d_50']
d_80 = npz_subhalo['d_80']
d_100 = npz_subhalo['d_100']

# (3) load gas data (location, properties)

gas_x = npz_gas['gas_x']
gas_y = npz_gas['gas_y']
gas_z = npz_gas['gas_z']
gas_distance_from_center = (np.array(gas_x - center[0])**2.0 + np.array(gas_y - center[1])**2.0 + np.array(gas_z - center[2])**2.0)**0.5
gas_distance_min = min(gas_distance_from_center) # [ckpc/h]
gas_distance_max = max(gas_distance_from_center) # [ckpc/h] 

gas_density = npz_gas['gas_density']
gas_mass = npz_gas['gas_mass']
gas_mag = npz_gas['gas_mag'] # three-dimensional (w/ direction)
gas_sfr = npz_gas['gas_sfr']
gas_u = npz_gas['gas_u']

gas_mag_strength = (np.sum(gas_mag**2.0,axis=1))**0.5 * mag_gauss # sqrt(Bx^2 + By^2 + Bz^2)
gas_num = npz_gas['gas_mass'] * num_gas # (gas mass) / (proton mass)

# (4) calculate gas properties 

# (4-1) turbulent magnetic field from internal energy

epsilon_B = 0.1
gas_turb_B = (epsilon_B * 4.0 * pi * gas_u * (km_to_cm)**2.0 * gas_density * density_cgs)**0.5

# (5) load BH data (location, properties)

bh_x = npz_bh['bh_x']
bh_y = npz_bh['bh_y']
bh_z = npz_bh['bh_z']

bh_acc = npz_bh['bh_Mdot']

### Below : not checked yet

# (6) properties within radius R

R_bins = 20
radius_list = np.logspace(log10(gas_distance_min),log10(gas_distance_max),R_bins)
volume_within_R = 4.0 * pi * (radius_list*length_kpc*kpc_to_cm)**3.0 / 3.0 # [cm^3]

sum_of_gas_num_within_R = np.zeros(R_bins)
sum_of_gas_mag_strength_within_R = np.zeros(R_bins)
sum_of_gas_turb_B_within_R = np.zeros(R_bins)
sum_of_gas_sfr_within_R = np.zeros(R_bins)
num_of_sim_particle_within_R = np.zeros(R_bins)

for r in range(R_bins):
    R = radius_list[r]
    # take all particles within R
    gas_num_within_R = gas_num[gas_distance_from_center<R]
    gas_mag_strength_within_R = gas_mag_strength[gas_distance_from_center<R]
    gas_turb_B_within_R = gas_turb_B[gas_distance_from_center<R]
    gas_sfr_within_R = gas_sfr[gas_distance_from_center<R]
    # take sum within R : X(<R)
    sum_of_gas_num_within_R[r] = sum(gas_num_within_R)
    sum_of_gas_mag_strength_within_R[r] = sum(gas_mag_strength_within_R)
    sum_of_gas_turb_B_within_R[r] = sum(gas_turb_B_within_R)
    sum_of_gas_sfr_within_R[r] = sum(gas_sfr_within_R)
    # number of simulated gas particles : N(<R) 
    num_of_sim_particle_within_R[r] = len(gas_num_within_R)

# (7) take diff 

mean_radius = (radius_list[1:] + radius_list[0:-1])/2.0
shell_width = np.diff(radius_list)

total_in_shell_volume = np.diff(volume_within_R)                            # V : volume of spherical shell
total_in_shell_sim_particle_num = np.diff(num_of_sim_particle_within_R)     # N : gas particle num within the shell

sum_in_shell_gas_num = np.diff(sum_of_gas_num_within_R)
sum_in_shell_gas_mag_strength = np.diff(sum_of_gas_mag_strength_within_R)
sum_in_shell_gas_turb_B = np.diff(sum_of_gas_turb_B_within_R)
sum_in_shell_gas_sfr = np.diff(sum_of_gas_sfr_within_R)

# (8) avg value within the shell

avg_gas_num_density_within_shell = sum_in_shell_gas_num / total_in_shell_volume                      # X / V
avg_gas_mag_strength_within_shell = sum_in_shell_gas_mag_strength / total_in_shell_sim_particle_num  # X / N 
avg_gas_turb_B_within_shell = sum_in_shell_gas_turb_B / total_in_shell_sim_particle_num              # X / N

# (9) physical quantities within the shell

# interaction time

T_pp_within_shell = (avg_gas_num_density_within_shell * sigma_pp * c)**(-1.0) / 3.15e7  # [yr]

# escape

eta_g = 1 

Dbohm_PeV_within_shell_from_mag_strength = eta_g * 108. * 3.086e18 * (1.e15/1.e17) * c / 3.0 / (avg_gas_mag_strength_within_shell * 1.e6)
Dbohm_PeV_within_shell_from_turb_B = eta_g * 108. * 3.086e18 * (1.e15/1.e17) * c / 3.0 / (avg_gas_turb_B_within_shell * 1.e6)
