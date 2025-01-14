# Channel Flow with Control

MPI-parallelized high-performance code for channel flow in Fortran 90, building upon an older formulation by Flores & Jimenez (2006), with the opposition control strategy. The code solves coupled equations for vorticity and wall-normal velocity using pseudo-spectral direct numerical simulations, with Fourier-Chebyshev spatial discretization, Runge-Kutta timestepping scheme and constant mass flux. The opposition control strategy is implemented as a proportionality between the boundary conditions on velocity at the wall and velocity at a given location above it. See Guseva & Jiménez (2022) for more details.


Contributors: A. Guseva, O. Flores and J. Jiménez     

REFERENCES 
- Guseva, A. and Jiménez, J. (2022) Linear instability and resonance effects in large-scale opposition flow control. J. Fluid Mech., 935, A35.   
- Kim,J., Moin,P. and Moser,R. (1987) Turbulence Statistics in fully developped channel flow at low Reynolds numbers, J. Fluid Mech., 177, 133-166

### Short description of files

- main.f90 - principal file reading initial conditions and calling the time loop  
- ctes3D - constants defining numerical resolution and locations of spectra
- in_hre - namelists with input parameters (Reynolds number, domain size, control parameters, input and output)
- cross.v2.f90 - timestepper; linear and nonlinear terms
- io.f90 - input and output
- cfdiff8.v7.f90 - compact finite difference subroutines
- modules.v7.f90 - definitions of common variables for the code parts
- fou3D.f, rftsingle.f, cftsingle.f - Fourier transform subroutines
- laps.v7.f - solving Poisson problem v'' - rK v = phi 
- makefile - gathers all the code files for compilation

### To use the code
- install MPI / Fortran MPI compilers
- set up compiler parameters in the makefile (e.g. mpiifort or mpif90)
- choose resolution in ctes3D
- compile the code with make CTRL
- set up flow parameters in in_hre
- run the code using your parallel architecture (e.g. mpirun -np 48 ./CTRL < in_hre)
