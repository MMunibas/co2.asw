CO2-CCSD(T)-F12 PES


Potential energy surfaces for the ground state of CO2:

TODO: add details on how the scan was performed.

TODO: Description of the files. 


**Requirements**
(1) RKHS toolkit : Download from https://github.com/MeuwlyGroup/RKHS

**Running the executable**

Compilation of the program: gfortran RKHS.f90  pesCO2.f90  test.f90


Before running the executable make sure that the asymp.coeff, .csv and/or .kernel files for that PES present in the current directory (or change the file path in the fortran program).
