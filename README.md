# Channel Flow with Control

MPI-parallelized high-performance code for channel flow in Fortran 90, building upon an older formulation by Flores & Jimenez (2006), with the opposition control strategy. The code solves coupled equations for vorticity and wall-normal velocity using pseudo-spectral direct numerical simulations, with Fourier-Chebyshev spatial discretization, Runge-Kutta timestepping scheme and constant mass flux. The opposition control strategy is implemented as a proportionality between the boundary conditions on velocity at the wall and velocity at a given location above it. See Guseva & Jiménez (2022) for more details.


Contributors: A. Guseva, O. Flores and J. Jiménez     

REFERENCES 
- Guseva, A. and Jiménez, J. (2022) Linear instability and resonance effects in large-scale opposition flow control. J. Fluid Mech., 935, A35.   
- Kim,J., Moin,P. and Moser,R. (1987) Turbulence Statistics in fully developped channel flow at low Reynolds numbers, J. Fluid Mech., 177, 133-166 
