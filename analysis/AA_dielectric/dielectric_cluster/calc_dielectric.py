# Script to calculate the dielectric of slabs within a simulation.
# Reads in an xtc file, centered on the largest cluster of protein atoms, and a tpr file for the topography
# calculates the dielectric across each of "nslab" slabs.
# Saves 1) the total dielectric of the system as a function of time ("eps_tot.txt"), and 2) the dielectric constant of each slab as a function of time ("eps_slab.txt").
import sys
import os
import math
import numpy as np
import MDAnalysis as mda

# Inputs. Will need to be modified to be used by others.
xtc_file='../prod_fixed_cluster.xtc'
tpr_file='../prod.tpr'
nslab=str(20)
GMX_path='/home/xclin/bin/gmx2019c/bin/'

# calculates dipole moment as a function of time for L frames, with dipoles given by dipole_tot and volumes given by V_list
def calc_dielectic(L,dipole_tot,V_list):
    eps = []
    for i in range(0, L):
        # split data to only include up to the current timestep
        dipole_tot2 = np.asarray(dipole_tot[0:i+1])
        V_list2 = np.asarray(V_list[0:i+1])
        #print(len(dipole_tot2))

        # subtract out mean, < (M - <M>) * (M - <M>) >
        ave = dipole_tot2.mean(0)
        subtracted = dipole_tot2 - ave
        #print(subtracted)
        dipole_variance = np.mean((subtracted * subtracted).sum(-1), 0)
        # include constants and fix units
        # unit conversions
        nm = 10 ** (-9)
        e = 1.602176634 * (10 ** (-19))
        kB = 1.380649 * (10 ** (-23))
        T = 313.15
        eps0 = 8.8541878128 * (10 ** (-12))
        dipole_variance = dipole_variance * (nm ** 2) * (e ** 2)
        V = np.mean(V_list2, 0)
        V = V * (nm ** 3)
        eps_timestep = 1 + (dipole_variance) / (3 * eps0 * V * kB * T)
        eps.append(eps_timestep)
    eps = np.asarray(eps)
    return eps


# runs gromacs to generate dipole moment
def run_gmx(t,XTC,TPR,nslab,GMX_path):
    slab_file=' slab_'+t+'.xvg '
    dipole_file='Mtot_'+t+'.xvg'
    os.system('rm \#*')
    os.system('echo 0 | '+GMX_path+'gmx dipoles -f '+XTC+' -s '+TPR+' -temp 313.5 -slab'+slab_file+' -axis Z -sl '+nslab+' -b '+t+' -e '+t+' -o '+dipole_file)
    return

# Determine basics of the system such as the number of atoms, trajectory length, time between timesteps, and box size.
def analyze_traj(XTC):
    u1 = mda.Universe(XTC)
    N = len(u1.atoms)
    print('Number of atoms: '+str(N))
    L = len(u1.trajectory)
    print('Length of trajectory: '+str(L))
    dt=u1.trajectory.dt
    print('Time between timesteps: '+str(dt))
    time0 = u1.trajectory.time
    print('Time trajectory begins: ' + str(time0))

    box_list=[]
    for ts in u1.trajectory:
        box=ts.dimensions
        box_list.append(box)
    box_list = np.asarray(box_list)*0.1 #convert to nm

    return N,L,dt,box_list,time0

# Read in dipoles and convert units
def read_dipole(filename):
    dipoles=np.loadtxt(filename,skiprows=27)
    return dipoles[1:4]*0.02083 # convert from D to nm*e

# read in lab file and convert units
def read_slab(filename):
    temp=np.loadtxt(filename,skiprows=26)
    slab=np.zeros((len(temp),4))
    slab[:,0]=temp[:,0]
    slab[:,1:4]=temp[:,1:4]*0.02083 # convert from D to nm*e
    return slab

# Main -------------------------------------------------------------------------------------------------------------
Natoms,Ntimesteps,dt,box,t0=analyze_traj(xtc_file)

# Determine the dipole moment for the 1) the entire system, and 2) each slab in the system.
dipoles=[]
slabs=[]
tot_volume=[]
for i in range(0,Ntimesteps):
    V=box[i,0]*box[i,1]*box[i,2]
    tot_volume.append(V)
    time=str(i*dt+t0)
    # calculate dipole moment at each timestep with GROMACS
    run_gmx(time, xtc_file, tpr_file, nslab,GMX_path)
    # Mtot is the dipole moment of the entire simulation box
    dipole_temp=read_dipole('Mtot_'+time+'.xvg')
    dipoles.append(dipole_temp)
    # slab is the dipole moment of each slab seperately
    slab_temp = read_slab('slab_'+time+'.xvg')
    slabs.append(slab_temp)

# calculate overall dipole of the system
dipoles = np.asarray(dipoles)
E=calc_dielectic(Ntimesteps,dipoles,tot_volume)
np.savetxt('eps_tot.txt',E)
print(E)

# calculate the dipole of each slab
E_slab=[]
for i in range(0,int(nslab)):
    dipole_slab=[]
    slab_V=[]
    for j in range(0,Ntimesteps):
        dipole_slab.append([slabs[j][i][1],slabs[j][i][2],slabs[j][i][3]])
        Z_axis=slabs[j][0][0]*2 # slab size in Z-axis
        V = box[j, 0] * box[j, 1] * Z_axis
        slab_V.append(V)
    dipole_slab = np.asarray(dipole_slab)
    slab_V = np.asarray(slab_V)
    eps_temp=calc_dielectic(Ntimesteps, dipole_slab, slab_V)
    E_slab.append(eps_temp)
E_slab = np.asarray(E_slab)
np.savetxt('eps_slab.txt',E_slab)
print(E_slab)







