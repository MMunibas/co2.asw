CO2-CCSD(T)-F12 PES


Potential energy surfaces for the ground state of CO2:

The grid used for scan is in Jacobi coordinates where r= distance between the diatomic(CO_A)[0.8, 2.1] Angstrom, R= distance between center of mass of diatomic and atom(O_B)[0.9, 5.9] Angstrom and theta= angle between R and r [0, 180] degree. 


**Requirements**
(1) RKHS toolkit : Download from https://github.com/MeuwlyGroup/RKHS

**Running the executable**

Compilation of the program: gfortran RKHS.f90  pesCO2.f90  test.f90


Before running the executable make sure that the asymp.coeff, .csv and/or .kernel files for that PES present in the current directory (or change the file path in the fortran program).
