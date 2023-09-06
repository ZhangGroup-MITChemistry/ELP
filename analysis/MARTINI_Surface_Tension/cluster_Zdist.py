# Script to calculate the size of clusters from slab simulations.
# Defined as proteins that are members whose center of mass is within "cut"
# distance of another protein in the given cluster
# Script runs off of a tpr file ("TPR"), which has information on the system, and
# a fixed xtc file ("XTC"), which ensures chains are not broken across pbc boundaries and excludes equilibration time for the initial run.
# the user must also specify the number of protein chains in the simulation box ("nchain")
# Written by Andrew Latham
# Note the inputs at top. These need to be customized for the system

# Output files:
# cluster_size.txt - Size of the 4 largest clusters, determined by the depth first search algorithm as a function of time
# cluster_list.txt - trajectory of whether each protein is in the largest cluster (1) or is not in the largest cluster (0)
# protein_solvent_hist.txt - density of protein / solvent from Z_pos, determined by relative Z-distance to the center of mass of the largest cluster

# Inputs
# tpr file used to run the simulation
TPR='prod_old.tpr'
# xtc file used to run the simulation
XTC='prod_fixed2.xtc'
# number of protein chains in the simulation box
nchain=40
# number of equilibration frames
eq=-1
# distance cuttoff
cut=40


import sys
import os
import math
import numpy
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
from MDAnalysis import transformations

# executes depth first search algorithm
def dfs_iterative(graph, start):
    stack, path = [start], []

    while stack:
        vertex = stack.pop()
        if vertex in path:
            continue
        path.append(vertex)
        if vertex in graph:
            for neighbor in graph[vertex]:
                stack.append(neighbor)

    return path

# converts adjacency matrix into dictionary
def mat_to_dict(mat):
    dict1={}
    for i in range(0,len(mat)):
        for j in range(0,len(mat[i])):
            if i in dict1 and mat[i][j]==1:
                dict1[i].append(j)
            elif i not in dict1 and mat[i][j]==1:
                dict1[i] = [j]
    return dict1

# function to handle pbc
def wrap_coord(pos,side):
    while pos<0:
        pos=pos+side
    while pos>side/2:
        pos=pos-side
    return pos


# main function of the code. Returns 1) the size of the largest cluster, 2) which proteins are in the largest cluster,
# and 3) the density of protein, water, and ions around the largest cluster.
def contact_mat(tpr_file,xtc_file,nchain,start,cutoff):
    u1 = mda.Universe(tpr_file, xtc_file)

    # caclulate number of atoms per protein chain and number of timesteps
    protien=u1.select_atoms("protein")
    N1=int(len(protien)/nchain)
    print('Number of protein atoms per chain: '+str(N1))
    timesteps=len(u1.trajectory)-(start+1)
    print('Number of timesteps: '+str(timesteps))
    N = len(u1.atoms)

    # counter for number of timesteps, initial empty pas
    count=0
    path_tot_timestep = []
    Z_mat2 = numpy.zeros((timesteps, 4))
    Z_pos = numpy.zeros((N,timesteps))


    # loop over all timesteps. only analyze frames after equilibration time
    for ts in u1.trajectory:
        print('frame: ' + str(ts.frame))
        if ts.frame < start:
            pass
        else:
            # Determine the adjacency matrix of all proteins
            com = numpy.zeros((nchain, 3))
            box = ts.dimensions
            mass = numpy.zeros((nchain, 1))
            # calculate com of each protein chain
            for i in range(0, nchain):
                index1 = i * N1
                index2 = ((i + 1) * N1)
                atoms = u1.atoms[index1:index2]
                com[i, :] = atoms.center_of_mass()
                mass[i] = atoms.total_mass()
            # convert to adjacency matrix
            mat = numpy.zeros((nchain, nchain))
            for i in range(0, nchain):
                for j in range(i + 1, nchain):
                    dist = distance_array(com[i, :], com[j, :], box)
                    if dist < cutoff:
                        mat[i, j] = mat[i, j] + 1
                        mat[j, i] = mat[i, j]

            # Calculate dfs at each timestep from matrix
            graph2 = mat_to_dict(mat)
            visited = []
            path_tot = []
            for i in range(0, len(mat)):
                flag = 0
                for k in range(0, len(visited)):
                    if visited[k] == i:
                        flag = 1

                if flag == 0:
                    path1 = []
                    path1 = dfs_iterative(graph2, i)
                    for j in range(0, len(path1)):
                        visited.append(path1[j])
                    path_tot.append(path1)
                path_tot_timestep.append(path_tot)
                #print(path_tot_timestep)

            # find atoms in largest cluster
            n = len(path_tot)
            max = 0
            l2 = -1
            index = -1
            #print(path_tot)
            l_list=[]
            for i in range(0, n):
                l = len(path_tot[i])
                l_list.append(l)
                if l > l2:
                    # reset old variables
                    l2 = l
                    index = i
            atoms_in_cluser = path_tot[index]
            mass_tot = 0
            for i in range(0, len(atoms_in_cluser)):
                index = atoms_in_cluser[i]
                mass_tot = mass_tot + mass[index]
            com2 = (mass[:, 0] / mass_tot) * numpy.cos((com[:, 2] / box[2]) * 2 * numpy.pi)
            com3 = (mass[:, 0] / mass_tot) * numpy.sin((com[:, 2] / box[2]) * 2 * numpy.pi)

            # Calculate COM of cluster. Use angles to avoid issues with periodicity
            Z1 = 0
            Z2 = 0
            # l_list: list of cluster sizes
            l_list=numpy.asarray(l_list)
            l_list[::-1].sort()
            l_list.resize((4))
            Z_mat2[count,:]=l_list
            print(Z_mat2[count,:])
            for i in range(0, l2):
                index = atoms_in_cluser[i]
                Z1 = Z1 + com2[index]
                Z2 = Z2 + com3[index]
            Z1 = Z1 / l2
            Z2 = Z2 / l2
            theta = numpy.arctan2(-1 * Z2, -1 * Z1) + numpy.pi
            Z_com = (box[2] / (2 * numpy.pi)) * theta


            #calculate wrapped list of atomistic distances to the largest cluster
            for i in range(0,N):
                atom = u1.atoms[i]
                Z = atom.position[2]
                Z_diff=Z-Z_com
                Z_final = wrap_coord(Z_diff, box[2])
                Z_pos[i,count]=abs(Z_final) / box[2]



            count = count + 1
            # steps to cut short for debugging
            #if count>10:
            #    break
    #Z_pos=Z_pos[:,0:count]

    # calculate density of protein / solvent from Z_pos
    # presets for making histogram
    nbins=50
    maxZ=0.5
    dZ=maxZ/nbins
    hist=numpy.zeros((nbins,5))
    # set X-axis
    for i in range(0,nbins):
        hist[i,0]=dZ*i+dZ/2
    hist[0,0]=hist[0,0]-dZ/2
    hist[nbins-1,0] = hist[nbins-1,0] + dZ / 2
    # number of water and ions
    nW=len(u1.select_atoms('name W'))
    nNA = len(u1.select_atoms('name NA'))
    nCL = len(u1.select_atoms('name CL'))
    #calculate histogram of atomic masses
    for i in range(0,N):
        for j in range(0, count):
            bin=int(Z_pos[i,j]/dZ)
            # atoms in protein
            if i<N1*nchain:
                hist[bin,1]=hist[bin,1]+u1.atoms.masses[i]
            # atoms in water
            elif i <N1*nchain+nW:
                hist[bin,2]=hist[bin,2]+u1.atoms.masses[i]
            # NA salt
            elif i <N1*nchain+nW+nNA:
                hist[bin,3]=hist[bin,3]+u1.atoms.masses[i]
            # CL salt
            else:
                hist[bin,4]=hist[bin,4]+u1.atoms.masses[i]
    # Header for output file
    key_string = '\tr\t\t\tProtein\t\t\tW\t\t\tNA\t\t\tCL'

    # calculate density of protein / solvent from Z_pos
    # Calculate number of residue types in protein
    res_types={}
    for i in range(0,N1):
        found=0
        res=u1.atoms.resnames[i]
        if res not in res_types:
            res_types[res]=len(res_types)
    #print(res_types)
    nres=len(res_types)
    # initialize histogram
    hist2 = numpy.zeros((nbins, nres+1))
    # set X-axis
    for i in range(0, nbins):
        hist2[i, 0] = dZ * i + dZ / 2
    hist2[0, 0] = hist2[0, 0] - dZ / 2
    hist2[nbins - 1, 0] = hist2[nbins - 1, 0] + dZ / 2
    # calculate histogram of atomic masses
    for i in range(0, N1*nchain):
        for j in range(0, count):
            bin = int(Z_pos[i, j] / dZ)
            resname=u1.atoms.resnames[i]
            hist_index=res_types[resname]+1
            hist2[bin, hist_index] = hist2[bin, hist_index] + u1.atoms.masses[i]
    # Header for output file
    key_list = list(res_types.keys())
    key_string2='\tr\t\t\t'
    for i in range(0,len(key_list)):
        key_string2=key_string2+key_list[i]+'\t\t\t'

    return Z_mat2,hist,hist2,key_string,key_string2


# run main function
cluster_size,hist,hist2,key_string,key_string2=contact_mat(TPR,XTC,nchain,eq,cut)

# save outputs
numpy.savetxt('cluster_size.txt',cluster_size)
numpy.savetxt('protein_solvent_hist.txt',hist,header=key_string)
numpy.savetxt('protein_residue_hist.txt',hist2,header=key_string2)

