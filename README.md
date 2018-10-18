# SpectralCodes

A collection of codes (in MATLAB & Fortran 77), and examples, for solving reaction-diffusion equations in one and two space dimensions. In areas of the mathematical community integrating factors together with spectral methods are used to remove the stiffness associated with the diffusive terms in a reaction-diffusion model allowing explicit high order timestepping to be used. This is particularly valuable for two (and higher) space dimension problems.

For a description of the codes please refer to the technical report: 
R. V. Craster and R. Sassi, *Spectral algorithms for reaction-diffusion equations*, Technical Report, Note del Polo No. 99, Università degli Studi di Milano, 2006. Handle:[`2434/24276`](http://hdl.handle.net/2434/24276). [arXiv:1810.07431 [math.NA]](https://arxiv.org/abs/1810.07431).
A PDF copy of the paper can be downloaded from [here](https://air.unimi.it/retrieve/handle/2434/24276/12468/NotePolo99.pdf).

## Getting Started

The collection of codes are placed into a directory with the following structure:

```
|-- Fortran
|   |-- OneD
|   |   |-- CK45
|   |   |   |-- auto_CK45.f
|   |   |   |-- epidemic_CK45.f
|   |   |   |-- fisher_CK45.f
|   |   |   `-- gray1D_CK45.f
|   |   |-- ETDRK4_B
|   |   |   |-- auto_ETDRK4_B.f
|   |   |   |-- epidemic_ETDRK4_B.f
|   |   |   |-- fisher_ETDRK4_B.f
|   |   |   `-- gray1D_ETDRK4_B.f
|   |   `-- RK4
|   |       |-- auto_RK4.f
|   |       |-- epidemic_RK4.f
|   |       |-- fisher_RK4.f
|   |       `-- gray1D_RK4.f
|   `-- TwoD
|       |-- CK45
|       |   |-- gray2D_CK45.f
|       |   `-- labyrinthe2D_CK45.f
|       |-- ETDRK4_B
|       |   |-- gray2D_ETDRK4_B.f
|       |   `-- labyrinthe2D_ETDRK4_B.f
|       `-- RK4
|           |-- gray2D_RK4.f
|           `-- labyrinthe2D_RK4.f
|-- Matlab
|   |-- OneD
|   |   |-- CK45
|   |   |   |-- auto_CK45.m
|   |   |   |-- epidemic_CK45.m
|   |   |   |-- fisher1D_CK45.m
|   |   |   `-- gray1D_CK45.m
|   |   |-- ETDRK4
|   |   |   `-- gray1D_ETDRK4.m
|   |   |-- ETDRK4_B
|   |   |   |-- auto_ETDRK4_B.m
|   |   |   |-- epidemic_ETDRK4_B.m
|   |   |   |-- fisher1D_ETDRK4_B.m
|   |   |   `-- gray1D_ETDRK4_B.m
|   |   `-- RK4
|   |       |-- auto_RK4.m
|   |       |-- epidemic_RK4.m
|   |       |-- fisher1D_RK4.m
|   |       `-- gray1D_RK4.m
|   `-- TwoD
|       |-- CK45
|       |   |-- fisher2D_CK45.m
|       |   |-- gray2D_CK45.m
|       |   `-- labyrinthe2D_CK45.m
|       |-- ETDRK4
|       |   `-- labyrinthe2D_ETDRK4.m
|       |-- ETDRK4_B
|       |   |-- fisher2D_ETDRK4_B.m
|       |   |-- gray2D_ETDRK4_B.m
|       |   `-- labyrinthe2D_ETDRK4_B.m
|       `-- RK4
|           |-- fisher2D_RK4.m
|           |-- gray2D_RK4.m
|           `-- labyrinthe2D_RK4.m
`-- Useful
    |-- adifisher.m
    |-- fourierupsample.m
    |-- fourierupsample2D.m
    |-- plot_fisher2D.m
    |-- plot_gray2D.m
    `-- plot_labyrinthe2D.m
```

The Fortran 77 and Matlab codes are separated and further sub-divided into One- and Two-dimensional codes. Each illustrative
example of the text is in each of the Cash-Karp adaptive (`CK45`), constant time step Exponential time differencing (`EDTRK4-B` and for
a few cases `EDTRK4`), and constant time step Runge-Kutta (`RK4`) folders.

Each code is provided as a completely self-contained unit and is internally commented. The Fortran 77 codes have been compiled
using several different compilers (g77, Intel, Salford, Compaq) and compile cleanly with no flags or options required.

There is additionally a directory, `Useful`, that contains the Matlab function routines to interpolate the data using FFTs, 
plotting functions and contains the ADI code from appendix C.

## Authors

The code were initially developed in 2003-2004 by:
* **Richard V. Craster** - [Imperial College, UK](https://www.imperial.ac.uk/people/r.craster)
* **Roberto Sassi** - [University of Milan, Italy](https://homes.di.unimi.it/sassi/)

See also the list of [contributors](https://github.com/roesassi/SpectralCodes/contributors) who participated in this project.

## Cited by

To the best of our knowledge, the following is an (incomplete) list of the scientific works which employed the codes:

* K.M. Owolabi and K.C. Patidar, *Effect of spatial configuration of an extended nonlinear Kierstead–Slobodkin reaction-transport model with adaptive numerical scheme*, (2016) SpringerPlus, 5 (1), art. no. 303. [DOI: 10.1186/s40064-016-1941-y](http://dx.doi.org/10.1186/s40064-016-1941-y)

* K.M. Owolabi and K.C. Patidar, *Solution of pattern waves for diffusive fisher-like non-linear equations with adaptive methods*, (2016) International Journal of Nonlinear Sciences and Numerical Simulation, 17 (6), pp. 291-304. [DOI: 10.1515/ijnsns-2015-0173](http://dx.doi.org/10.1515/ijnsns-2015-0173)

* K.M. Owolabi and K.C. Patidar, *Higher-order time-stepping methods for time-dependent reaction-diffusion equations arising in biology*, (2014) Applied Mathematics and Computation, 240, pp. 30-50. [DOI: 10.1016/j.amc.2014.04.055](http://dx.doi.org/10.1016/j.amc.2014.04.055)

* K.M. Owolabi and K.C. Patidar, *Numerical solution of singular patterns in one-dimensional Gray-Scott-like models*, (2014) International Journal of Nonlinear Sciences and Numerical Simulation, 15 (7-8), pp. 437-462. [DOI: 10.1515/ijnsns-2013-0124](http://dx.doi.org/10.1515/ijnsns-2013-0124)

* M.A. Dahlem and T.M. Isele, *Transient localized wave patterns and their application to migraine*, (2013) Journal of Mathematical Neuroscience, 3 (1), pp. 1-28. [DOI: 10.1186/2190-8567-3-7](http://dx.doi.org/10.1186/2190-8567-3-7)

* R.M. Bradley and P.D. Shipman, *A surface layer of altered composition can play a key role in nanoscale pattern formation induced by ion bombardment*, (2012) Applied Surface Science, 258 (9), pp. 4161-4170. [DOI: 10.1016/j.apsusc.2011.07.003](http://dx.doi.org/10.1016/j.apsusc.2011.07.003)

* R. Carlson, *Harmonic analysis for star graphs and the spherical coordinate trapezoidal rule*, (2011) Journal of Computational and Applied Mathematics, 235 (8), pp. 2077-2089. [DOI: 10.1016/j.cam.2010.10.006](http://dx.doi.org/10.1016/j.cam.2010.10.006)

* N. Apreutesei and G. Dimitriu, *On a prey-predator reaction-diffusion system with Holling type III functional response*, (2010) Journal of Computational and Applied Mathematics, 235 (2), pp. 366-379. [DOI: 10.1016/j.cam.2010.05.040](http://dx.doi.org/10.1016/j.cam.2010.05.040)



