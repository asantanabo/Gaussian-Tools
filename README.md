# Gaussian tools

This is a collection of scripts to automatize or analyze output produced by the Quantum Chemistry code Gaussian. 
The first code (frag.py) split a xyz file into molecular fragments for further calculations such as the BSSE. 

The second code analyze the transition dipole moment and the magnetic dipole moment produced in a simulation of
ECD spectrum. AS a result the code produces a .vmd file which can be easily visulized with vmd -e command.

The third code is a refinement of the Projective method used to compute the mobility within the marcus theory formalism. 
The code reads 3 Gaussian output files: Molecule A, Molecule B and Molecule AB. As a result the program produces the values 
of the transfer integrals for selected energy levels.
