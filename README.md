# SpectralCodes

A collection of codes (in MATLAB & Fortran 77), and examples, for solving reaction-diffusion equations in one and two space dimensions. In areas of the mathematical community integrating factors together with spectral methods are used to remove the stiffness associated with the diffusive terms in a reaction-diffusion model allowing explicit high order timestepping to be used. This is particularly valuable for two (and higher) space dimension problems.

For a description of the codes please refer to the technical report: 
R. V. Craster and R. Sassi, *Spectral algorithms for reaction-diffusion equations*, Technical Report, Note del Polo – Ricerca N.99, Università di Milano, 2006. Handle:[`2434/24276`](http://hdl.handle.net/2434/24276). [arXiv:1810.07431 [math.NA]](https://arxiv.org/abs/1810.07431).
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

