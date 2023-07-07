One-body Potential Pins -- version 0.1
Written by Trung

===== DESCRIPTION

This routine diagonalize within a given conformal Hilbert space a one-body Hamiltonian. The Hamiltonian is a set of N direct delta Delta functions located at x_1, x_2, ..., x_N:

H = delta^2(x_1) + delta^2(x_2) + ... + delta^2(x_N)

The conformal Hilbert space is specified by specifying a list of Jack polynomials. 

This routine works on the disk geometry.

===== SOURCE CODE
    PotentialPins.jl
    Potentials.jl

Also requires FQH_states_v2.jl, Misc.jl, and Density.jl from _qhe-library

Other accompanying routine:
    admissibility_full.jl
    jack_get.py
    jack


===== INSTRUCTION
 0) Before running the routine, make sure to prepare two things: (i) a single file containing all the root configurations for the Jacks that span the desired Hilbert space, and (ii) All the corresponding jack polynomials stored in separate files in a subfolder called "jacks". Each jack must be named as "J_<root configuration>"

 	        For example, to get a one-quasihole Hilbert space of Laughlin state with two electrons, there has to be a file with the following content:
 	        	3
 	        	01001
 	        	0.0
 	        	10001
 	        	0.0
 	        	10010
 	        	0.0
                
 	        And inside the folder "jacks" there must be three files named "J_01001", "J_10001", and "J_10010"

    If the above requirements are not yet fulfilled, they may be prepared by first running

julia admissibility_full.jl

    This will generate a list of admissibile root given a number of electrons and a number of orbitals. And the list of jacks can be generated with.

python3 jack_get.py

    (See NOTE 1 below for more details)

 1) After the above is ready, the main routine can be run with

julia PotentialPins.jl

The program will prompt the user to input the followings:

    Input root file name:
> Input the name of the file that contains all root configurations

    How many potential pins?
> Input the number of potential traps

    Input location of x y of each potential pins
> Input each pair of (x,y) coordinate for each potential pins. For example, if one wants to place a pin at x=2, y=-1.4, then input "2 1.4" (without quotation marks). Press Enter/Return to register each pair.


2) Suppose <filename> is the name of the file that contains all the root configurations and <N> is the number of pins, the eigenvalues will be stored in a file named <filename>_gs_<N>pins_eigen. Two lowest eigenstates will also be saved as <filename>_gs_<N>pins_0 and <filename>_gs_<N>pins_1.

===== NOTES
1) For larger system sizes, it may take very long to generate all the jacks. Before running python3 jack_get.py, check the folder /media/sdd2/Jack_polynomials to see whether the jacks needed are already there.

2) For the mathematical derivation behind this routine, see https://hqtr.notion.site/One-body-potential-0068c6263da84e08b5c14cfbbbd959fc