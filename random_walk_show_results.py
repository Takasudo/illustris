# visualize results from random_walk.py

from math import *
import numpy as np

# basic data

snap = 33
redshift = 2.
z0id = 0
calclabel = 'test'

runtmp = 'test'

Epev = [1.0]
haloID = [0]

run_num = "run_"+runtmp+"_snap"+str(snap)

mag_flag = 0

if mag_flag==0:
    mag_label="B"
elif mag_flag==1:
    mag_label="turbB"
    
# -------------------- #
# [I] prepare for calc #
# -------------------- #

# (1) load npz files

npz_subhalo = np.load("npzdata_protocluster_subhalo_snap"+str(snap)+"_z0id"+str(z0id)+"_"+str(calclabel)+".npz")
npz_gas = np.load("gas_snap"+str(snap)+"_z0id"+str(z0id)+"_"+str(calclabel)+".npz")
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

subhalo_idlist = [i for i in range(len(subhalo_sfr))]

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
gas_u = npz_gas['gas_u']       # 300-3
gas_id = npz_gas['gas_id']     # 300-3

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

# ----------------------------------------------------------------------- #
# [II] import matplotlib and prepare function to show gas and subhalo map #
# ----------------------------------------------------------------------- #

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as patches

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes

def showmap(ax1,ax2,ax3):

    # gas 

    Hist = ax1.hist2d(gas_x,gas_y,weights=gas_density,bins=170,norm=mpl.colors.LogNorm(),cmap=mpl.cm.hot)
    fig.colorbar(Hist[3],ax=ax1)
    ax1.set_xlabel('x [ckpc/h]')
    ax1.set_ylabel('y [ckpc/h]') 

    Hist = ax2.hist2d(gas_y,gas_z,weights=gas_density,bins=170,norm=mpl.colors.LogNorm(),cmap=mpl.cm.hot)
    fig.colorbar(Hist[3],ax=ax2)
    ax2.set_xlabel('y [ckpc/h]')
    ax2.set_ylabel('z [ckpc/h]')

    Hist = ax3.hist2d(gas_x,gas_z,weights=gas_density,bins=170,norm=mpl.colors.LogNorm(),cmap=mpl.cm.hot)
    fig.colorbar(Hist[3],ax=ax3)
    ax3.set_xlabel('x [ckpc/h]')
    ax3.set_ylabel('z [ckpc/h]')

    # subhalo

    pos_1 = subhalo_pos.T[0]
    pos_2 = subhalo_pos.T[1]
    pos_3 = subhalo_pos.T[2]
    
    for m in range (len(subhalo_pos)):
        c1 = patches.Circle(xy=(pos_1[m], pos_2[m]), fill=False, radius=subhalo_size[m], ec='grey',linewidth=1.2)
        c2 = patches.Circle(xy=(pos_2[m], pos_3[m]), fill=False, radius=subhalo_size[m], ec='grey',linewidth=1.2)
        c3 = patches.Circle(xy=(pos_1[m], pos_3[m]), fill=False, radius=subhalo_size[m], ec='grey',linewidth=1.2)
        ax1.add_patch(c1)
        ax2.add_patch(c2)
        ax3.add_patch(c3)
        for i in range(len(HaloID)):
            if subhalo_idlist[m]==HaloID[i]:
                ax1.text(pos_1[m],pos_2[m],str(subhalo_idlist[m]),color='blue',fontsize=8)
                ax2.text(pos_2[m],pos_3[m],str(subhalo_idlist[m]),color='blue',fontsize=8)
                ax3.text(pos_1[m],pos_3[m],str(subhalo_idlist[m]),color='blue',fontsize=8)
                c1 = patches.Circle(xy=(pos_1[m], pos_2[m]), fill=True, radius=subhalo_size[m], ec='blue',linewidth=1.5)
                c2 = patches.Circle(xy=(pos_2[m], pos_3[m]), fill=True, radius=subhalo_size[m], ec='blue',linewidth=1.5)
                c3 = patches.Circle(xy=(pos_1[m], pos_3[m]), fill=True, radius=subhalo_size[m], ec='blue',linewidth=1.5)
                ax1.add_patch(c1)
                ax2.add_patch(c2)
                ax3.add_patch(c3)

# -------------------------------------------------------------------------------- #
# [III] load results from crDiff npz file & Show results for each set of halo/Epev #
# -------------------------------------------------------------------------------- # 

# save quantities for later

subhalo_sfr_list = [0 for h in range(len(HaloID))]
avg_fpp_Myr = [[0 for e in range (len(Epev))] for h in range(len(HaloID))]
avg_fpp_10Myr = [[0 for e in range (len(Epev))] for h in range(len(HaloID))]
avg_fpp_90Myr = [[0 for e in range (len(Epev))] for h in range(len(HaloID))]

rgb=0

for e in range(len(Epev)):
    for h in range(len(HaloID)):

        # load data

        filename = "run_"+runtmp+"/crDiff_"+str(mag_label)+"_"+run_num+"_halo"+str(HaloID[h])+"_Epev"+str(Epev[e])+".npz"
        npz_cr = np.load(filename,allow_pickle=True)
        haloSFR = npz_cr['subhaloSFR']
        cr_E = npz_cr['crEpev']
        num = npz_cr['crnum']
        cr_time = npz_cr['trace_time'] 
        fpp = npz_cr['trace_fpp']
        mag = npz_cr['trace_mag'] 
        trace_ngas = npz_cr['trace_ngas']
        trace_R = npz_cr['trace_Reff']
        timestep = npz_cr['trace_timestep'] 
        lc_time = npz_cr['trace_lc_time'] 
        cr_x = npz_cr['trace_cr_x'] 
        cr_y = npz_cr['trace_cr_y'] 
        cr_z = npz_cr['trace_cr_z'] 
        fpp_Myr = npz_cr['list_fpp_Myr'] 
        fpp_3Myr = npz_cr['list_fpp_3Myr'] 
        fpp_10Myr = npz_cr['list_fpp_10Myr'] 
        fpp_30Myr = npz_cr['list_fpp_30Myr'] 
        fpp_50Myr = npz_cr['list_fpp_50Myr'] 
        fpp_90Myr = npz_cr['list_fpp_90Myr']

        # save some quantities

        subhalo_sfr_list[h] = haloSFR
        avg_fpp_Myr[h][e] = np.average(fpp_Myr)
        avg_fpp_10Myr[h][e] = np.average(fpp_10Myr)
        avg_fpp_90Myr[h][e] = np.average(fpp_90Myr)

        # prepare color

        for_color = 0
        color_tmp = [[0,0,0] for n in range(num)]

        for n in range(num):
            for_color +=1
            if for_color % 3 == 0:
                rgb = 0
            elif for_color % 3 == 1:
                rgb = 1
            elif for_color %3 == 2:
                rgb = 2
            color_tmp[n][rgb] += n/float(num)

        # (1) fig for f_pp

        figname = "run_"+runtmp+"/fig_fpp_"+str(mag_label)+"_"+run_num+"_halo"+str(HaloID[h])+"_Epev"+str(Epev[e])+".pdf"
        fig = plt.figure(figsize=(9,18))

        # (1a) evolution of f_pp with time

        ax1 = fig.add_subplot(211)
        for n in range(num):
            ax1.loglog(cr_time[n],fpp[n],color='cyan',alpha=0.7,linewidth=0.7)
        ax1.set_xlabel("Time [Myr]",fontsize=15)
        ax1.set_ylabel("$f_{pp}$",fontsize=15)
        ax1.loglog([1.],np.median(fpp_Myr),marker='*',ms=14,color='blue')
        ax1.loglog([1.],np.average(fpp_Myr),marker='p',ms=11,color='red')
        ax1.loglog([3.],np.median(fpp_3Myr),marker='*',ms=14,color='blue')
        ax1.loglog([3.],np.average(fpp_3Myr),marker='p',ms=11,color='red')
        ax1.loglog([10.],np.median(fpp_10Myr),marker='*',ms=14,color='blue')
        ax1.loglog([10.],np.average(fpp_10Myr),marker='p',ms=11,color='red')
        ax1.loglog([30.],np.median(fpp_30Myr),marker='*',ms=14,color='blue')
        ax1.loglog([30.],np.average(fpp_30Myr),marker='p',ms=11,color='red')
        ax1.loglog([90.],np.median(fpp_90Myr),marker='*',ms=14,color='blue',label="median")
        ax1.loglog([90.],np.average(fpp_90Myr),marker='p',ms=11,color='red',label="avg")
        ax1.text(min(cr_time[n]),min(1.0,np.max(fpp_90Myr))*0.95,"Subhalo SFR = "+str(int(haloSFR))+" [$M_\odot$/yr]",fontsize=13)
        ax1.text(min(cr_time[n]),min(1.0,np.max(fpp_90Myr))*0.75,"CR Energy = "+str(cr_E)+" [PeV]",fontsize=13)
        plt.legend(loc='lower right')

        # (1b) Histogram of f_pp 

        ax2 = fig.add_subplot(212)
        ax2.set_xlabel("log$_{10}$($f_{pp}$)",fontsize=16)
        ax2.set_ylabel("N",fontsize=16)
        fppbins = np.linspace(-6,0,49)
        ax2.hist(np.log10(fpp_Myr),bins=fppbins,color='black',linewidth=4,histtype='step',label="1 Myr") 
        ax2.hist(np.log10(fpp_10Myr),bins=fppbins,color='blue',linewidth=4,histtype='step',label="10 Myr")
        ax2.hist(np.log10(fpp_90Myr),bins=fppbins,color='red',linewidth=4,histtype='step',label="90 Myr")
        ax2.hist(np.log10(fpp_Myr),bins=fppbins,color='black',alpha=0.3)
        ax2.hist(np.log10(fpp_10Myr),bins=fppbins,color='blue',alpha=0.4)
        ax2.hist(np.log10(fpp_90Myr),bins=fppbins,color='red',alpha=0.5)
        plt.legend(loc='upper left')

        plt.savefig(figname)
        plt.close()

        # (2) map

        mapname = "run_"+runtmp+"/map_"+str(mag_label)+"_"+run_num+"_halo"+str(HaloID[h])+"_Epev"+str(Epev[e])+".pdf"
        fig = plt.figure(figsize=(18,16))
        ax1 = fig.add_subplot(221)
        ax2 = fig.add_subplot(224)
        ax3 = fig.add_subplot(223)
        showmap(ax1,ax2,ax3)
        for n in range(num):
            ax1.plot(cr_x[n],cr_y[n],color="cyan",alpha=0.7,linewidth=0.7)
            ax2.plot(cr_y[n],cr_z[n],color="cyan",alpha=0.7,linewidth=0.7)
            ax3.plot(cr_x[n],cr_z[n],color="cyan",alpha=0.7,linewidth=0.7) 
        
        ax2.set_title("Trajectory for Halo Id ="+str(HaloID[h])+" (SFR = "+str(int(haloSFR))+") and E_pev = "+str(cr_E))
        
        # (2a) zoom

        ax1_zoom = zoomed_inset_axes(ax1, 4, loc=2)
        plt.xticks(visible=False)
        plt.yticks(visible=False)
        ax2_zoom = zoomed_inset_axes(ax2, 4, loc=2)
        plt.xticks(visible=False)
        plt.yticks(visible=False)
        ax3_zoom = zoomed_inset_axes(ax3, 4, loc=2)
        plt.xticks(visible=False)
        plt.yticks(visible=False)
        for n in range(num):
            ax1_zoom.plot(cr_x[n],cr_y[n],color=color_tmp[n],alpha=0.7,linewidth=0.5)
            ax2_zoom.plot(cr_y[n],cr_z[n],color=color_tmp[n],alpha=0.7,linewidth=0.5)
            ax3_zoom.plot(cr_x[n],cr_z[n],color=color_tmp[n],alpha=0.7,linewidth=0.5)
        ax1_zoom.set_xlim(cr_x[n][0]-500/length_kpc,cr_x[n][0]+500/length_kpc)
        ax1_zoom.set_ylim(cr_y[n][0]-500/length_kpc,cr_y[n][0]+500/length_kpc)
        ax2_zoom.set_xlim(cr_y[n][0]-500/length_kpc,cr_y[n][0]+500/length_kpc)
        ax2_zoom.set_ylim(cr_z[n][0]-500/length_kpc,cr_z[n][0]+500/length_kpc)
        ax3_zoom.set_xlim(cr_x[n][0]-500/length_kpc,cr_x[n][0]+500/length_kpc)
        ax3_zoom.set_ylim(cr_z[n][0]-500/length_kpc,cr_z[n][0]+500/length_kpc)

        p_1 = subhalo_pos.T[0][HaloID[h]]
        p_2 = subhalo_pos.T[1][HaloID[h]]
        p_3 = subhalo_pos.T[2][HaloID[h]]
        rad = subhalo_size[HaloID[h]]
        print(rad*length_kpc)
        c1 = patches.Circle(xy=(p_1, p_2), fill=False, radius=rad, ec='red',linewidth=1.)
        c2 = patches.Circle(xy=(p_2, p_3), fill=False, radius=rad, ec='red',linewidth=1.)
        c3 = patches.Circle(xy=(p_1, p_3), fill=False, radius=rad, ec='red',linewidth=1.)        
        ax1_zoom.add_patch(c1)
        ax2_zoom.add_patch(c2)
        ax3_zoom.add_patch(c3)

        # (2b) distance travelled

        ax4 = fig.add_subplot(222)
        for n in range(num):
            distance_from_start = (np.array(cr_x[n]-cr_x[n][0])**2.0 + np.array(cr_y[n]-cr_y[n][0])**2.0 + np.array(cr_z[n]-cr_z[n][0])**2.0)**0.5
            ax4.loglog(fpp[n],distance_from_start*length_kpc,color=color_tmp[n],alpha=0.7,linewidth=0.7)
            ax4.set_xlabel("$f_{pp}$")
            ax4.set_ylabel("Distance [kpc]")

        plt.savefig(mapname,bbox_inches='tight')
        plt.close()

        # (3) trace CR trajectory
        
        figname = "run_"+runtmp+"/fig_track_"+str(mag_label)+"_"+run_num+"_halo"+str(HaloID[h])+"_Epev"+str(Epev[e])+".pdf"
        fig = plt.figure(figsize=(9,18))

        ax1 = fig.add_subplot(311)
        for n in range(num):
            ax1.loglog(cr_time[n],fpp[n],color=color_tmp[n],alpha=0.7,linewidth=0.7)
        ax1.set_xlabel("Time [Myr]",fontsize=15)
        ax1.set_ylabel("$f_{pp}$",fontsize=15)

        ax2 = fig.add_subplot(312)
        for n in range(num):
            ax2.loglog(cr_time[n],trace_ngas[n],color=color_tmp[n],alpha=0.7,linewidth=0.7)
        ax2.set_xlabel("Time [Myr]",fontsize=15)
        ax2.set_ylabel("$n_{gas}$",fontsize=15)

        ax3 = fig.add_subplot(313)
        for n in range(num):
            ax3.loglog(cr_time[n],mag[n],color=color_tmp[n],alpha=0.7,linewidth=0.7)
        ax3.set_xlabel("Time [Myr]",fontsize=15)
        ax3.set_ylabel("$B$",fontsize=15)

        plt.savefig(figname,bbox_inches='tight')
        plt.close()
        
        # note

        t2 = time.time()
        print("Epev : ",str(e+1),"/",str(len(Epev))," Halo : ",str(h+1),"/",str(len(HaloID))," time : ",t2-t1)
