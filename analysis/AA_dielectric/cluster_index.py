# Script to write indicies of the proteins in the largest cluster.
# Defined as proteins that are members whose center of mass is within "cut"
# distance of another protein in the largest cluster for "cluster_cut" fraction of the simulation.
# Script runs off of a tpr file ("TPR"), which has information on the system, and
# a fixed xtc file ("XTC"), which ensures chains are not broken across pbc boundaries and excludes equilibration time for the initial run.
# the user must also specify the number of protein chains in the simulation box ("nchain")
# The script also outputs the density from slab simulations, which was used to verify our choice of cutoff.
# Written by Andrew P. Latham
# Note the inputs at top. These need to be customized for the system.

# Output files:
# index_cluster.ndx - index file with 2 groups. 1. System, is all atoms in the system. 2. Cluster is all the atoms in proteins in the largest cluster
# cluster_size.txt - Size of the 4 largest clusters, determined by the depth first search algorithm as a function of time
# cluster_list.txt - trajectory of whether each protein is in the largest cluster (1) or is not in the largest cluster (0)
# protein_solvent_hist.txt - density of protein / solvent from Z_pos, determined by relative Z-distance to the center of mass of the largest cluster


# Inputs
# tpr file used to run the simulation
TPR='prod.tpr'
# xtc file used to run the simulation
XTC='prod_fixed.xtc'
# number of protein chains in the simulation box
nchain=40
# number of equilibration frames
eq=-1
# distance cuttoff
cut=40
# fraction of simulation in the largest cluster to include as part of the cluster
cluster_cut=0.99


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


# main function of the code. Determines which atoms are in the largest cluster, saves the indicies of all atoms for
# proteins in the largest cluster, and returns 1) the size of the largest cluster, 2) which proteins are in the largest cluster,
# and 3) the density of protein, water, and ions around the largest cluster.
def contact_mat(tpr_file,xtc_file,nchain,start,cutoff,cluster_cut):
    u1 = mda.Universe(tpr_file, xtc_file)

    # caclulate number of atoms per protein chain and number of timesteps
    protien=u1.select_atoms("protein")
    N1=int(len(protien)/nchain)
    print('Number of protein atoms per chain: '+str(N1))
    timesteps=len(u1.trajectory)-(start+1)
    print('Number of timesteps: '+str(timesteps))
    N = len(u1.atoms)

    # counter for number of timesteps, initialize output variables
    count=0
    path_tot_timestep = []
    Z_mat2 = numpy.zeros((timesteps, 4))
    Z_pos = numpy.zeros((N,timesteps))
    cluster_list=numpy.zeros((nchain,timesteps))

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
            atoms_in_cluster = path_tot[index]
            mass_tot = 0
            for i in range(0, len(atoms_in_cluster)):
                index = atoms_in_cluster[i]
                mass_tot = mass_tot + mass[index]
            com2 = (mass[:, 0] / mass_tot) * numpy.cos((com[:, 2] / box[2]) * 2 * numpy.pi)
            com3 = (mass[:, 0] / mass_tot) * numpy.sin((com[:, 2] / box[2]) * 2 * numpy.pi)

            # Loop through atoms_in_cluster. Assign 0 to not in cluster, 1 to in cluster
            for i in range(0,len(atoms_in_cluster)):
                index=atoms_in_cluster[i]
                cluster_list[index,count]=1
            print(cluster_list[:,count])


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
                index = atoms_in_cluster[i]
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
            #if count>1:
                #break
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
    # generate histrogram for protein positions
    # 1. get all protein atoms. 2. Loop over all protein atoms and all timesteps.
    # 3. Get the Z_position and mass of the protein atoms
    indexP=u1.select_atoms('protein').indices
    for i in range(0,len(indexP)):
        for j in range(0, count):
            index=indexP[i]
            bin = int(Z_pos[index, j] / dZ)
            hist[bin, 1] = hist[bin, 1] + u1.atoms.masses[index]
    # repeat for water
    indexW = u1.select_atoms('name OW').indices
    for i in range(0, len(indexW)):
        for j in range(0, count):
            index = indexW[i]
            bin = int(Z_pos[index, j] / dZ)
            hist[bin, 2] = hist[bin, 2] + u1.atoms.masses[index]
    # repeat for NA
    indexNA = u1.select_atoms('name SOD').indices
    for i in range(0, len(indexNA)):
        for j in range(0, count):
            index = indexNA[i]
            bin = int(Z_pos[index, j] / dZ)
            hist[bin, 3] = hist[bin, 3] + u1.atoms.masses[index]
    # repeat for CL
    indexCL = u1.select_atoms('name CLA').indices
    for i in range(0, len(indexCL)):
        for j in range(0, count):
            index = indexCL[i]
            bin = int(Z_pos[index, j] / dZ)
            hist[bin, 4] = hist[bin, 4] + u1.atoms.masses[index]
    # Header for output file
    key_string = '\tr\t\t\tProtein\t\t\tOW\t\t\tSOD\t\t\tCLA'

    # determine the fraction of simulation for which proteins are in the largest cluster
    proteins_in_cluster=numpy.sum(cluster_list,1)/count
    atoms_list=[]
    # loop over all proteins, if they are in the largest cluster for at least "cluster_cut" frequency,
    # select their atoms and add them to atoms_list (the list of all indecies of atoms in the largest cluster)
    for i in range(0,len(proteins_in_cluster)):
        if proteins_in_cluster[i]>cluster_cut:
            start_chain=(N1)*(i)
            end_chain = (N1) * (i+1)
            for j in range(start_chain,end_chain):
                atoms_list.append(j+1)
    #print(atoms_list)
    print(len(atoms_list))
    # write index file
    # output the entire system
    new=open('index_cluster.ndx','w')
    new.write('[ System ]\n')
    count_atoms=0
    # split each line at 10 atoms for readability
    for i in range(1,N+1):
        if count_atoms%10==9:
            new.write(str(i)+'\n')
        else:
            new.write(str(i)+' ')
        count_atoms = count_atoms+1
    # Add the atoms of the proteins in the largest cluster to the output
    new.write('\n\n[ Cluster ]\n')
    count_atoms = 0
    for i in range(0,len(atoms_list)):
        if count_atoms % 10 == 9:
            new.write(str(atoms_list[i]) + '\n')
        else:
            new.write(str(atoms_list[i]) + ' ')
        count_atoms = count_atoms + 1
    new.write('\n')
    new.close()

    return Z_mat2,hist,key_string,cluster_list


# run main function
cluster_size,hist,key_string,cluster_list=contact_mat(TPR,XTC,nchain,eq,cut,cluster_cut)

# save outputs
numpy.savetxt('cluster_size.txt',cluster_size)
numpy.savetxt('cluster_list.txt',cluster_list)
numpy.savetxt('protein_solvent_hist.txt',hist,header=key_string)

