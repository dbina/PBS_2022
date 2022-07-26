###### PBS_2022
###### PBS modelling code
###### readme file version 2022-07-21 by D.B.

This set of .m files performs analysis of energy transfer in
phycobilisome, based on structure saved as .pdb file within the framework of Foerster theory
the theory behind the computations can be found in the accompanying publication
[REFERENCE TO THE PAPER HERE, I suppose]

A. workflow
1. .pdb is translated to coordinate file (.xyz) of coordinates of c
arbon atoms of bilins (ascii) and cysteines (these are used to assign the bilin binding)

2. using sequence data (can be extracted from the same pdb), bilins are assigned to 
spectral classes (= given a numerical code) using the <dictionary_v01.txt> file

3. spectra for pigment classes and their overlap integrals are generated based on bilin spectra from literature. The Boltzmann 
equilibrium condition for backward (uphill) energy transfer is applied in the form of modification of the overlap integral
overlaps are stored as a matrix (J) whose i-th,j-th element corresponds to the overlap integral between 
donor pigment of class i and acceptor pigment of class j
this can also output the pigment spectra for simulation of the transient spectroscopy data

4. coordinate data (.xyz) are read and converted to molecular transition dipoles
the resulting file format is in the script codes denoted as mRe variable. This stands for: mju-R-e, i.e.:
(transition) dipole size ([debye])-position vector (nanometers)-unit vector in the direction of mju, that is mRe is a matrix:
wih columns [dipole in debye,x,y,z,x,y,z]
the dipole data are extracted from bilin coordinates by performing singular value decomposition on the coordinates
to find the geometrical axis of the conjugated system;
the same approach is used to extract dipole for the carotenoids; which atoms to use was determined manually by inspection of the pigment 
.pdb structure, therefore two different scripts exist, one for bilins, the other for carotenoids (..._car.m)
the mRe are saved as <*_vecdata.txt> files that are then loaded and processed to compute dipole couplings etc.
also generated is a file <index_....txt> which contains numerical codes for spectral classes of pigments
so that a line in the mRe matrix corresponds to a line in the index_ file

5. mRe are loaded from the files and used to compute dipole-dipole couplings using the dipole corrdinates and sizes, combined with 
overlaps these give the transfer rates. Resulting matrix of pairwise rate constants is fed into the ordinary differential equation 
solver to yield kinetics of excited states. Initial conditionas can be chosen to simulate  excitation into different pigments.
kinetics are saved as text files

6. kinetics can be combined with the spectra (step 3) to simulate transient spectroscopy (fluorescence) data, that can then be 
processed to yield decay-associated spectra (DAS). Alternatively, DAS can be computed directly using rate constant matrix. However, fitting the simulated data like one would a normal dataset generates more experiment-like DAS directly.  

B. the scripts (numbers correspond to workflow stages given in A)
1. pdb2xyz.m
2. identify_hetresidues.m uses <dictionary_v0.txt> to assign classes to pigments
pdb_2_seq.m is provided as a tool to convert pdb file to fasta file of protein sequence, so that only pdb is needed as input for the 
work
3. generate_sample_spectra_01_boltzmann.m
4. convert_xyz_to_mRe_vecfile(_car).m
5. PBS_w_overlaps_and_ode_gill.m 
(the _gill part stands for the option to run a Gillespie algorithm simulation instead of solving ode; this e.g. allows to incorporate
structural dynamics into the system for stochastic simulation; not used in the original publication)
this script contains the actual simulation code; the steps 1-4 have to be run only once to generate the vector files,spectra and indexes
it outputs the excited state concentrations as a function of time and, as an option the matrix of the pairwise rate constants for the whole pigment system

6. spectra_from_ode.m - outputs the spectra and DAS

notes:
i) most scripts require that their folder contains a text file <current_path.txt> with a path to the folder where the structure data are stored
ii) follow comments inside the script code for more explanations
iii) the scrips were developed as a tool to help its author analyze a particular dataset, hence with minimal consideration for user interfacing, consequently parameters such as file names for output, values of constants, simulation settings etc. are to be input directly into the code
and optional functionality is achieved by commenting/uncomenting section of code

