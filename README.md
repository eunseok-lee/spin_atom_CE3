##########################################################
### Spin-Atom Cluster Expansion for quinary alloys
#### written by Eunseok Lee
#### v1: May 23, 2018
#### v2: June 5, 2018 (major update)
##########################################################

This computational package is to perform cluster expansion (CE) for quinary alloys such as Li_pNi_qMn_rCo_sVa_{2-p-q-r-s}O2 (NMC cathode for Li ion batteries). The cluster functions are formulated from the coupled configuration of atomic species and magnetic moment at each lattice site. For further information, refer to Physical Review B 95, 085134 (2017) or http://atom.uah.edu.

The package consists of five sub-programs, listed in the following list. Each program can run independently if all required parameters are provided correctly.

1. clusterlist: formulate the clusters and cluster functions based on the geometrical information of lattice site.
2. data_to_corr_mat3: convert the coupled configuration of atomic species and magnetic moment to the correlation matrix of CE.
3. findcluster3: select the most representative cluster functions (expansion basis) and the corresponding effective cluster interactions (expansion coefficient).
4. predictstructure_ce3: predict the lowest energy structure using the result of 3.
5. link_to_ann: construct and optimize an artificial neural network to fit correlation matrix to formation energy.

The main function of the package is realized by sub-programs 3 and 4. You may use sub-program 5 instead of sub-program 3, but its performance will be inferior to the one of sub-program 3 if the nunmber of data sets is small.

The every code was written in C. The most time consuming part of the package comes from the calculation of the correlation matrix. Hence, sub-programs 2, 3, and 4 were developed in parallel version, using MPI. Sub-program 1 is light to run and a serial program. 

* * *
### Installation 
Requirement: mpicc, GSL, qhull

GNU Scientific Library (GSL) was used in several parts of the package and hence the installation of GSL is pre-requisite for compilation (refer to https://www.gnu.org/software/gsl/ for further information).

The library from Qhull was also used to construct the convex hull and to calculate the formation energy above the convex hull. Hence, the user needs to install Qhull on their local machine and add the location of library to the PATH. The user can do this by downloading and extracting qhull package, 'make', 'make install', and 'export LD_LIBRARY_PATH=$PWD/lib:$LD_LIBRARY_PATH' (refer to http://www.qhull.org for detailed information). 

For comparison, spin_atom_CE package (not spin_atom_CE3) uses its own code to construct the convex hull, not the library from Qhull.

Although the sub-programs 2~4 will be compiled by mpicc, they can run on single-cpu.

To compile each sub-program, move in the src_(sub-program) and follow the direction in README there.

This version targets R-3m space group materials, such as the layered Li_[p]Ni_[q]Mn_[r]Co_[s]Va_[2-p-q-r-s]O_2. Each cationic lattice site is be occupied by any of Li, Ni, Mn, Co, or Va and can have any of up-, down-, and zero-magnetic moment. Hence, the total degrees of freedom (DOF) is 5x3=15.

Although they developed for the layered Li_[p]Ni_[q]Mn_[r]Co_[s]Va_[2-p-q-r-s]O_2, the sub-programs 2)~4) can also be applied for any material system with 5x3 DOF at each lattice site if the classification of cluster functions and mapping from atomic clusters to cluster function are provided.


* * *
### Updates
- k-NN algorithm for anomaly detection was implemented into data_to_corr_mat (sub-program 2).
- Backpropagation in artificial neural network architectures was implemended into link_to_ann (sub-program 5) in version 2.
