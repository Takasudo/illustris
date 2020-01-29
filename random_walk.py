# random walk on illustris snapshot

from math import *
import numpy as np
import illustris_python as il

# basic data

snap = 33
redshift = 2.
z0id = 0
calclabel = 'test'

runtmp = 'test'

Epev = [1.0]
haloID = [0]

run_num = "run_"+runtmp+"_snap"+str(snap)

####

mag_flag = 0
# 0 : for B
# 1 : for turb B

if mag_flag==0:
    mag_label="B"
elif mag_flag==1:
    mag_label="turbB"

filename_log = ("run_"+runtmp+"/tmp_"+str(mag_label)+"_"+run_num+"_log.dat")
file_log = open(filename_log,'w')

print("")
print("### run number : ",(args[1]))
print("")
file_log.write("### run number : "+(args[1])+"\n")

# -------------------- #
# [I] prepare for calc #
# -------------------- #

# (0) parameters

N_cr = 300

# (1) load npz files

npz_subhalo = np.load("npzdata_protocluster_subhalo_snap"+str(snap)+"_z0id"+str(z0id)+"_"+str(calclabel)+".npz")
npz_gas = np.load("gas_snap"+str(snap)+"_z0id"+str(z0id)+"_"+str(calclabel)+".npz") # 300-3
npz_bh = np.load("npzdata_protocluster_bh_snap"+str(snap)+"_z0id"+str(z0id)+"_"+str(calclabel)+".npz")
npz_all_subhalo = np.load('subhalo_snap'+str(snap)+"_z0id"+str(z0id)+"_"+str(calclabel)+".npz")

# (2) prepare for unit conversion

num_gas = (1e10 / hubble) * 1.989e33 / 1.66e-24       # from mass [10^10 M_sun/h] to gas num [Particle Number]
mag_gauss = hubble * (1.0 + redshift)**2.0 * 2.60e-6  # from [(h/a^2)(Unit Mass/Unit Length)^{1/2}/(Unit Time)] to [cgs Gauss]
length_kpc = 1./(1.+redshift)/hubble                  # from [ckpc/h] to [kpc]
km_to_cm = 1.e5
kpc_to_cm = 3.086e21
density_cgs = (1e10 / hubble) * 1.989e33 / (length_kpc*kpc_to_cm)**3.0 # from [(10^10 M_sun/h)/(ckpc/h)^3] to [g/cm^3]
lightspeed = 3.e10 # [cm/s]

# (3) subhalo : center and size of protocluster

center = npz_subhalo['CoM']
d_20 = npz_subhalo['d_20']
d_50 = npz_subhalo['d_50']
d_80 = npz_subhalo['d_80']
d_100 = npz_subhalo['d_100']

subhalo_pos = npz_all_subhalo['subhalo_pos']
subhalo_sfr = npz_all_subhalo['subhalo_sfr']
subhalo_mass = npz_all_subhalo['subhalo_mass']
subhalo_size = npz_all_subhalo['subhalo_size']
subhalo_bh = npz_all_subhalo['subhalo_bh']

# (4) gas : location and properties

gas_x = npz_gas['gas_x']
gas_y = npz_gas['gas_y']
gas_z = npz_gas['gas_z']

gas_distance_from_center = (np.array(gas_x - center[0])**2.0 + np.array(gas_y - center[1])**2.0 + np.array(gas_z - center[2])**2.0)**0.5
gas_distance_min = min(gas_distance_from_center) # [ckpc/h]
gas_distance_max = max(gas_distance_from_center) # [ckpc/h] 

gas_density = npz_gas['gas_density']
gas_mass = npz_gas['gas_mass']
gas_mag = npz_gas['gas_mag'] * mag_gauss # three-dimensional (w/ direction)
gas_sfr = npz_gas['gas_sfr']
gas_u = npz_gas['gas_u']     # 300-3
gas_id = npz_gas['gas_id']   # 300-3

gas_mag_strength = (np.sum(gas_mag**2.0,axis=1))**0.5
gas_Reff = (3.0 * gas_mass/ (4.0*pi) / gas_density)**(1./3.) # effective size

epsilon_B = 0.01
gas_turb_B = (epsilon_B * 4.0 * pi * gas_u * (km_to_cm)**2.0 * gas_density * density_cgs)**0.5  # 300-3

gas_Lcoh = gas_Reff * length_kpc # coherent length of the magnetic field [kpc]

# (5) gas : characterize diffusion

D_crit = gas_Lcoh * (lightspeed / kpc_to_cm) / 3.0 # [kpc^2/s]
Tdiff_crit = (gas_Reff * length_kpc)**2.0 / (D_crit) / 3.15e7                   # [yr] R^2/D
Tcross = (2.0 * gas_Reff * length_kpc * kpc_to_cm) / (lightspeed) / 3.15e7      # [yr] 2R/D
Tpp = 1./(0.5 * gas_density * density_cgs / 1.66e-24 * 3.e10 * 3.e-26 * 3.15e7) # [yr]
fpp_crit = Tdiff_crit / Tpp
fpp_cross = Tcross / Tpp

# --------------------------------------- #
# [II] random walk through gas particles  #
# --------------------------------------- #

from random import uniform, choice,choices

# parameters

E_cr = [Epev] * N_cr # [PeV]

if mag_flag==0:
    max_time = 1.e8 # [yr]
elif mag_flag==1:
    max_time = 1.e9

# list to save data

time_list = [[] for n in range(N_cr)]
gas_id_list = [[] for n in range(N_cr)]
gas_fpp_list = [[] for n in range(N_cr)]
gas_mag_list = [[] for n in range(N_cr)]
gas_n_list = [[] for n in range(N_cr)] # trace n_gas along CR trajectory 
gas_R_list = [[] for n in range(N_cr)] # trace R_eff along CR trajectory
gas_x_list = [[] for n in range(N_cr)]
gas_y_list = [[] for n in range(N_cr)]
gas_z_list = [[] for n in range(N_cr)]
gas_lc_list = [[] for n in range(N_cr)] # List of Light Crossing Time (for code check)
cr_x_list = [[] for n in range(N_cr)]
cr_y_list = [[] for n in range(N_cr)]
cr_z_list = [[] for n in range(N_cr)]
d_from_init_pos = [[] for n in range(N_cr)]

time_cumsum = [[] for n in range(N_cr)]
fpp_cumsum = [[] for n in range(N_cr)]

fpp_Myr = []
fpp_3Myr = []
fpp_10Myr = []
fpp_30Myr = [] 
fpp_50Myr = []
fpp_90Myr = []

print(" Halo ID = ",haloID," Halo SFR = ",subhalo_sfr[haloID])
file_log.write("Halo ID = "+str(haloID)+" Halo SFR ="+str(subhalo_sfr[haloID])+"\n")

bndsize = 100./length_kpc  # [ckpc/h]

for n in range(N_cr):

    # print

    print("")
    print("START n= ", n)
    print("Energy E_cr [PeV] = ",E_cr[n])
    file_log.write("\n")
    file_log.write("START n= "+str(n)+"\n")
    file_log.write("Energy E_cr [PeV] = "+str(E_cr[n])+"\n")

    # (1) set initial condition

    random_subhalo_pos = subhalo_pos[haloID]
    #tmp_x = random_subhalo_pos[0]
    #tmp_y = random_subhalo_pos[1]
    #tmp_z = random_subhalo_pos[2]
    subhalo_radius = subhalo_size[haloID] 
    tmp_x = random_subhalo_pos[0] + choice([-bndsize,bndsize]) # 1029 bnd tmp
    tmp_y = random_subhalo_pos[1] + choice([-bndsize,bndsize]) # 1029 bnd tmp 
    tmp_z = random_subhalo_pos[2] + choice([-bndsize,bndsize]) # 1029 bnd tmp 
    subhalo_distance_from_gas_particles = (gas_x - tmp_x)**2.0 + (gas_y - tmp_y)**2.0 + (gas_z - tmp_z)**2.0
    d_min = min(subhalo_distance_from_gas_particles)
    tmp = np.where(subhalo_distance_from_gas_particles == d_min)[0][0]

    # initial gas properties that contains cr particle (cr_gas_...)

    cr_gas_id = gas_id[tmp] # 300-3 
    cr_gas_x = gas_x[tmp]
    cr_gas_y = gas_y[tmp]
    cr_gas_z = gas_z[tmp]
    cr_gas_Reff = gas_Reff[tmp]
    cr_gas_Lcoh = gas_Lcoh[tmp]
    cr_gas_ngas = gas_density[tmp]
    cr_gas_Tdiff_crit = Tdiff_crit[tmp]
    cr_gas_fpp_crit = fpp_crit[tmp]
    cr_gas_Tcross = Tcross[tmp]
    cr_gas_fpp_cross = fpp_cross[tmp]

    if mag_flag==0:
        cr_gas_b_strength = gas_mag_strength[tmp]
        cr_gas_b_direction = gas_mag[tmp] / cr_gas_b_strength
    elif mag_flag==1:
        cr_gas_b_strength = gas_turb_B[tmp] # 300-3
        cr_gas_b_direction = [1./sqrt(3),1./sqrt(3),1./sqrt(3)] # temp !!

    # initial cr positions / velocity (cr_...)

    cr_x = cr_gas_x
    cr_y = cr_gas_y
    cr_z = cr_gas_z

    theta1 = uniform(0,2.0*pi)
    theta2 = uniform(0,2.0*pi)

    cr_v = np.array([cos(theta1)*cos(theta2),cos(theta1)*sin(theta2),sin(theta1)])

    # initial cr Larmor radius

    b_perp = cr_gas_b_strength * sqrt(1 - np.dot(cr_v, cr_gas_b_direction))
    R_g = 1.08e-3 * E_cr[n] / (b_perp * 1.e6) / length_kpc  # Larmor radius [ckpc/h]

    # initial cr timestep (= diffusion time in gas particle)

    if cr_gas_Lcoh > R_g:
        cr_timestep = cr_gas_Tdiff_crit * (R_g/cr_gas_Lcoh)**(-1./3.)
        cr_gas_fpp = cr_gas_fpp_crit * (R_g/cr_gas_Lcoh)**(-1./3.)
    else:
        cr_timestep = cr_gas_Tcross
        cr_gas_fpp = cr_gas_fpp_cross

    # print 

    print("Initial properties: R_L/L_coh = ",int(1.e5*R_g/cr_gas_Lcoh)/1.e5, "time_step/max_time = ",int(1.e4*cr_timestep/max_time)/1.e4)
    file_log.write("Initial properties: R_L/L_coh = ")
    file_log.write(str(int(1.e5*R_g/cr_gas_Lcoh)/1.e5))
    file_log.write("time_step/max_time = ")
    file_log.write(str(int(1.e4*cr_timestep/max_time)/1.e4))
    file_log.write("\n")

    # save initial id and position

    init_id = cr_gas_id # 300-3 
    init_x = cr_x
    init_y = cr_y
    init_z = cr_z

    # (2) random walk

    total_time = 0. # [yr]
    i = 0
    j = 0

    while total_time < max_time:

        # [1] increase timestep

        total_time += cr_timestep

        # [2] add values in list

        gas_id_list[n].append(cr_gas_id) # 300-3 
        gas_x_list[n].append(cr_gas_x)
        gas_y_list[n].append(cr_gas_y)
        gas_z_list[n].append(cr_gas_z)
        gas_mag_list[n].append(cr_gas_b_strength)
        gas_n_list[n].append(cr_gas_ngas)
        gas_R_list[n].append(cr_gas_Reff)
        gas_lc_list[n].append(cr_gas_Tcross)
        cr_x_list[n].append(cr_x)
        cr_y_list[n].append(cr_y)
        cr_z_list[n].append(cr_z)
        d_from_init_pos[n].append(((cr_x-init_x)**2.0 + (cr_y-init_y)**2.0 + (cr_z-init_z)**2.0)**0.5)
        time_list[n].append(cr_timestep)

        # [3] move CR particle

        if cr_gas_Lcoh > R_g:
        # Small Larmor Radius : random motion
            theta1 = uniform(-1,1) * pi
            theta2 = uniform(-1,1) * pi
        else:
        # Large Larmor Radius : small angle scattering
            deflect_angle = cr_gas_Lcoh/R_g
            theta1 += deflect_angle * choice([-1,1])
            theta2 += deflect_angle * choice([-1,1])

        # [4] new gas position and velocity

        cr_v = np.array([cos(theta1)*cos(theta2),cos(theta1)*sin(theta2),sin(theta1)])

        cr_x += cr_gas_Reff * cr_v[0]
        cr_y += cr_gas_Reff * cr_v[1]
        cr_z += cr_gas_Reff * cr_v[2]

        # [5] new gas particle ID and properties

        cr_distance_from_gas_particles = (cr_x - gas_x)**2.0 + (cr_y - gas_y)**2.0 + (cr_z - gas_z)**2.0
        d_min = min(cr_distance_from_gas_particles)
        index_min = np.where(cr_distance_from_gas_particles == d_min)[0][0]

        # [6] gas properties for ID == index_min

        cr_gas_id = gas_id[index_min] # 300-3 
        cr_gas_x = gas_x[index_min]
        cr_gas_y = gas_y[index_min]
        cr_gas_z = gas_z[index_min]
        cr_gas_Reff = gas_Reff[index_min]
        cr_gas_Lcoh = gas_Lcoh[index_min]
        cr_gas_ngas = gas_density[index_min]
        cr_gas_Tdiff_crit = Tdiff_crit[index_min]
        cr_gas_fpp_crit = fpp_crit[index_min]
        cr_gas_Tcross = Tcross[index_min]
        cr_gas_fpp_cross = fpp_cross[index_min]

        if mag_flag == 0:
            cr_gas_b_strength = gas_mag_strength[index_min]
            cr_gas_b_direction = gas_mag[index_min] / cr_gas_b_strength
        elif mag_flag == 1:
            cr_gas_b_strength = gas_turb_B[index_min] # 300-3 
            cr_gas_b_direction = [1. / sqrt(3), 1. / sqrt(3), 1. / sqrt(3)]  # temp !!

        # [7] new cr Larmor radius and diffusion time for ID == index_min

        b_perp = cr_gas_b_strength * sqrt(1 - np.dot(cr_v, cr_gas_b_direction))
        R_g = 1.08e-3 * E_cr[n] / (b_perp * 1.e6) / length_kpc  # Larmor radius [ckpc/h]

        if cr_gas_Lcoh > R_g:
            # resonant scattering
            cr_timestep = cr_gas_Tdiff_crit * (R_g / cr_gas_Lcoh) ** (-1. / 3.)
            cr_gas_fpp = cr_gas_fpp_crit * (R_g / cr_gas_Lcoh) ** (-1. / 3.)
        else:
            # diffusion not applicable
            cr_timestep = cr_gas_Tcross
            cr_gas_fpp = cr_gas_fpp_cross

        # [8] back to step 1

        i += 1
        gas_fpp_list[n].append(cr_gas_fpp)
        fpp_cum = np.sum(gas_fpp_list[n])
        fpp_cum_prevstep = np.sum(gas_fpp_list[n]) - cr_gas_fpp
        fpp_cum_avg_two_steps = (fpp_cum+fpp_cum_prevstep)/2.0

        # [8-1] : see f_pp at 1, 3, 10, 30, 50, 90 Myr

        if j==0 and total_time>1.e6:
            j+=1
            fpp_Myr.append(min(1,fpp_cum_avg_two_steps))
        if j==1 and total_time>3.e6:
            j+=1
            fpp_3Myr.append(min(1,fpp_cum_avg_two_steps))
        if j==2 and total_time>10.e6:
            j+=1
            fpp_10Myr.append(min(1,fpp_cum_avg_two_steps))
        if j==3 and total_time>30.e6:
            j+=1
            fpp_30Myr.append(min(1,fpp_cum_avg_two_steps))
        if j==4 and total_time>50.e6:
            j+=1
            fpp_50Myr.append(min(1,fpp_cum_avg_two_steps))
        if j==5 and total_time>90.e6:
            j+=1
            fpp_90Myr.append(min(1,fpp_cum_avg_two_steps))

        # [8-2] : break if fpp > 1

        if fpp_cum > 1:
            print("calc finished because f_pp > 1 ")
            file_log.write("calc finished because f_pp > 1 ")
            if j==0:
                fpp_Myr.append(1)
                fpp_3Myr.append(1)
                fpp_10Myr.append(1)
                fpp_30Myr.append(1)
                fpp_50Myr.append(1)
                fpp_90Myr.append(1)
            elif j==1:
                fpp_3Myr.append(1)
                fpp_10Myr.append(1)
                fpp_30Myr.append(1)
                fpp_50Myr.append(1)
                fpp_90Myr.append(1)
            elif j==2:
                fpp_10Myr.append(1)
                fpp_30Myr.append(1)
                fpp_50Myr.append(1)
                fpp_90Myr.append(1)
            elif j==3:
                fpp_30Myr.append(1)
                fpp_50Myr.append(1)
                fpp_90Myr.append(1)
            elif j==4:
                fpp_50Myr.append(1)
                fpp_90Myr.append(1)
            elif j==5:
                fpp_90Myr.append(1)  
            break

    time_cumsum[n] = np.cumsum(time_list[n])/1.e6 # [Myr]
    fpp_cumsum[n] = np.cumsum(gas_fpp_list[n])

    t2 = time.time()
    print("END n= ",n," Time [s] = ",t2-t1," Step Num = ",i)
    file_log.write("END n= "+str(n)+" Time [s] = "+str(t2-t1)+" Step Num = "+str(i)+"\n")

print("")
file_log.close()

# --------------------- #
# [III] output results  #
# --------------------- #

np.savez(
"run_"+runtmp+"/crDiff_"+str(mag_label)+"_"+run_num+"_halo"+str(haloID)+"_Epev"+str(Epev),
subhaloID = haloID,
subhaloSFR = subhalo_sfr[haloID],
subhaloLOC = subhalo_pos[haloID],
crEpev = Epev,
crnum = N_cr,
trace_time = np.array(time_cumsum),
trace_fpp = np.array(fpp_cumsum), 
trace_mag = np.array(gas_mag_list),
trace_Reff = np.array(gas_R_list),
trace_ngas = np.array(gas_n_list), 
trace_timestep = np.array(time_list),  # yr
trace_lc_time = np.array(gas_lc_list), # yr
trace_cr_x = np.array(cr_x_list),
trace_cr_y = np.array(cr_y_list),
trace_cr_z = np.array(cr_z_list),
list_fpp_Myr = np.array(fpp_Myr),
list_fpp_3Myr = np.array(fpp_3Myr),
list_fpp_10Myr = np.array(fpp_10Myr),
list_fpp_30Myr = np.array(fpp_30Myr),
list_fpp_50Myr = np.array(fpp_50Myr),
list_fpp_90Myr = np.array(fpp_90Myr)
)
