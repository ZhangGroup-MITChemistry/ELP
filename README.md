# Code to calculate the surface tension and micropolarity from multiscale simulations of ELP condensates
By Andrew P. Latham and Bin Zhang
If you use this code, please cite:
Songtao Ye, Andrew P. Latham, Chia-Heng Hsiung, Yuqi Tang, Junlin Chen, Feng Luo, Yu Liu, Bin Zhang, Xin Zhang, Micropolarity governs the structural organization of biomolecular condensates, Nature Chemical Biology, accepted manuscript (https://www.biorxiv.org/content/10.1101/2023.03.30.534881v1).

Here, we walk through an example of the simulations and analysis used in our manuscript.
All steps for performing simulations and analysis are demonstrated using the V10 system.
For steps involving GROMACS simulations, we demonstrate how to prepare tpr files, which can then be input to mdrun.

Finally, for those interested in performing AA simulations of condensates, we recommend exploring our updated procedure for
MOFF configurations to all-atom (AA) structures (see OpenABC, https://github.com/ZhangGroup-MITChemistry/OpenABC).

# MOFF Simulations ----------------------------------------------------------------------------------------------
// MOFF simulations are used to set up initial configuration of the cluster
// Found in: simulation/MOFF
// Run energy minimization / equilibration at NPT of polyV
gmx_mpi grompp -f min.mdp -c start.gro -p polyV.top -o min.tpr
gmx_mpi grompp -f eq.mdp -c min.gro -p polyV.top -o eq.tpr  -maxwarn 1
# Run energy minimization / equilibration at NVT
# Note the switch to add sequence specificity at this step 
gmx_mpi grompp -f min.mdp -c eq.gro -p ELP.top -o min_eq.tpr
gmx_mpi grompp -f eq3.mdp -c min_eq.gro -p ELP.top -o eq3.tpr  -maxwarn 1
# Center the output configuration, and convert to AA.
echo 0 0 | gmx_mpi trjconv -f eq3.gro -s eq3.tpr -pbc mol -o eq3.pdb -center
# Convert to all atom. This configuration is used as input to MARTINI. Note to ensure the box dimensions specified in the leap file align with those from the pdb file
tleap -f build.leap

# MARTINI Simulations -------------------------------------------------------------------------------------------
# MARTINI simulations are used to calculate the surface tension
# Found in: simulation/MARTINI
# Run MARTINI to generate initial configuration
martinize2 -f V10_AA1.pdb -nt -o PROA.top -ff martini3001 -x V10_MARTINI1.pdb -maxwarn 25
# Run MARTINI on a single chain to ensure the interpretation of the topology is correct
martinize2 -f V10_AA3.pdb -nt -o PROA.top -ff martini3001 -x V10_MARTINI3.pdb
# edit box (XXX- pbc dimensions of the PDB)
gmx_mpi editconf -f V10_MARTINI1.pdb -c -o V10_MARTINI.pdb -box XXX XXX 40.000
# solvate
gmx_mpi solvate -cp V10_MARTINI.pdb -cs water.gro -o solvated.gro -p PROA.top -radius 0.21
# add ions (specify number of ions to neutralize protein if charged)
gmx_mpi grompp -f ions.mdp -c solvated.gro -p PROA.top -o ions.tpr
gmx_mpi genion -s ions.tpr -o start.gro -p PROA.top -pname NA -nname CL -conc 1 
# run minimization
gmx_mpi grompp -f step4.1_minimization.mdp -c start.gro -p PROA.top -o min.tpr
# run nvt
gmx_mpi grompp -f step4.2_equilibration.mdp -c min.gro -p PROA.top -o nvt.tpr
# run npt
gmx_mpi grompp -f step4.3_equilibration.mdp -c nvt.gro -t nvt.cpt -p PROA.top -o npt.tpr
# run production- half 1
gmx_mpi grompp -f step5.0_production.mdp -c npt.gro -t npt.cpt -p PROA.top -o prod.tpr
# run production- half 2
gmx_mpi grompp -f step5.0_production.mdp -c prod.gro -t prod.cpt -p PROA.top -o prod2.tpr
# combine the two runs into a single xtc file
gmx_mpi trjcat -f prod.xtc prod2.xtc -settime -o trj_tot.xtc
# center on the largest cluster and ensure molecules are not split across pbc
echo 1 1 0 | gmx_mpi trjconv -f trj_tot.xtc -b 25000000.000 -center -pbc cluster -o prod_cluster.xtc -s prod.tpr
echo 0 | gmx_mpi trjconv -f prod_cluster.xtc -s prod.tpr -pbc mol -o prod_fixed2.xtc
# extract the desired frame
echo 0 | gmx_mpi trjconv -f prod_fixed2.xtc -s prod.tpr -o prod_100.gro -b 100000000 -e 100000000
# Run the backward script on this frame to get a starting AA configuration
./initram-v5.sh -f V10_100.gro -o V10_100_CHARMM.gro -to charmm36 -p topol.top -em 1000 -nb 5000 -np 16
# If necessary, write old tpr file that can be read in by MDAnalysis
gmx_mpi grompp -f step5.0_production.mdp -c npt.gro -p PROA.top -o prod_old.tpr

# AA Simulations ------------------------------------------------------------------------------------------------
# All-atom simulations are used to determine the micropolarity
# Found in: simulation/AA
# Run simulations. Topology of 1 chain was previously generated via CHARMM-GUI
# run energy minimization
gmx_mpi grompp -f step4.0_minimization.mdp -c V10_CHARMM.gro -r V10_CHARMM.gro -p topol.top -o min.tpr
# run nvt equilibration
gmx_mpi grompp -f step4.1_equilibration.mdp -c min.gro -r min.gro -p topol.top -o nvt.tpr
# run first npt equilibration
gmx_mpi grompp -f step4.2_equilibration.mdp -c nvt.gro -t nvt.cpt -p topol.top -o npt1.tpr
# run second npt equilibration
gmx_mpi grompp -f step4.3_equilibration.mdp -c npt1.gro -t npt1.cpt -p topol.top -o npt2.tpr
# run production
gmx_mpi grompp -f step5_production.mdp -c npt2.gro -t npt2.cpt -p topol.top -o prod.tpr
# ensure molecules are not split across the pbc to begin analysis, and exclude equilibration time
gmx_mpi trjconv -f prod.xtc -s prod.tpr -pbc mol -center -o prod_fixed.xtc -b 50000.000

# Analyze surface tension (MARTINI) -----------------------------------------------------------------------------
# Code to analyze surface tension
# Found in: analysis/MARTINI_Surface_Tension
# Surface tension * nsurf can be extracted directly from the edr file (not included due to file size constraints). This is done in 2 parts:
# Note: The trajectory example stored in this folder has been greatly shortened due to storage considerations. However, the results shown are for the full trajectory
gmx_mpi energy -f prod.edr -o surf_tens.xvg
gmx_mpi energy -f prod2.edr -o surf_tens2.xvg
# The number of surfaces for each configuration is calculated by cluster_Zdist.py. 
# The number of surfaces is defined as the number of clusters with 10 or more proteins (which represents 25% of the total system size).
python cluster_Zdist.py
# At each time point, the surface tension of the system is the surface_tension*nsurf from GROMACS divided by the number of clusters in the system

# Analyze Dielectric (AA) ---------------------------------------------------------------------------------------
# Code to analyze dielectric constant
# Found in: analysis/AA_dielectric
# Determine which indicies correspond to the largest protein cluster
# Note: The trajectory example stored in this folder has been greatly shortened due to storage considerations. However, the results shown are for the full trajectory
python cluster_index.py
# Create xtc file that centers the system on the largest cluster
echo 1 0 | gmx_mpi trjconv -f prod.xtc -s prod.tpr -pbc mol -center -o prod_fixed_cluster.xtc -b 50000.000 -n index_cluster.ndx 
# Run dielectric calculation
python calc_dielectric.py
