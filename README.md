# Uncertainty quantification in Bayesian earthquake source inversions (HPC ready)
Routines for evaluating approximate covariance matrices of Green's functions, designed for uncertainty quantification in Bayesian inversions of earthquake sources (ready for high-performance computing).
***************************************

This repository provides a Fortran90 and MATLAB routines for evaluation of the 
covariance (and cross-covariance) matrices of Green's functions. It is 
specifically designed as a high-performance toolset to efficiently handle 
model uncertainties in Bayesian waveform-based earthquake source inversions. 
It features a dual-language implementation to balance computational power 
(Fortran90) with ease of testing and visualization (MATLAB). List of four 
types of methods implemented in this package:
*   Approximate Covariance Function (ACF)
*   Approximate Cross-covariance Function (AXCF)
*   Stationarized Approximate Covariance Function (SACF)
*   Stationarized Approximate Cross-covariance Function (SAXCF)

1 METHODOLOGY
===================

  Hallo, M., Gallovič, F. (2016). Fast and cheap approximation of Green
functions uncertainty for waveform-based earthquake source inversions,
Geophys. J. Int., 207, 1012-1029. [https://doi.org/10.1093/gji/ggw320](https://doi.org/10.1093/gji/ggw320)

2 TECHNICAL IMPLEMENTATION
===================

Fortran90, High Performance Computing (HPC), Fourier Transform, Cross-Platform (Windows, Linux)

*   **High Performance:** Core routines implemented in Fortran90 for efficient evaluation of large covariance and cross-covariance matrices.
*   **Dual-Language Implementation:** Fortran90 functions are supplemented with respective MATLAB functions, using identical mathematical logic for seamless testing and prototyping.
*   **Efficient Uncertainty Quantification:** Specifically optimized for handling model uncertainties in Bayesian source inversion frameworks.
*   **Full Matrix Output:** The toolset returns complete covariance structures ready for integration into inverse solvers.

The official software version is archived on Zenodo:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19343279.svg)](https://doi.org/10.5281/zenodo.19343279)

3 PACKAGE CONTENT
===================

  1. `approxc.f90` - Fortran90 subroutines for determining covariance matrices (ACF, AXCF, SACF, SAXCF)
  2. `axcf.m` - MATLAB function for determining covariance matrices (ACF, AXCF)
  3. `saxcf.m` - MATLAB function for determining stationarized covariance matrices (SACF, SAXCF)
  4. `example.m` - Example of MATLAB code using these functions

4 REQUIREMENTS
===================

  FORTRAN: Codes fulfill Fortran 90 Standard (Modern Fortran).

  MATLAB: Version R2025b
  
  MATLAB Toolboxes: Matlab Curve Fitting Toolbox (`smooth.m`), Matlab Signal Processing Toolbox (`filtfilt.m`)

5 COPYRIGHT
===================

Copyright (C) 2016-2018  Miroslav Hallo and František Gallovič

This program is published under the GNU General Public License (GNU GPL).

This program is free software: you can modify it and/or redistribute it
or any derivative version under the terms of the GNU General Public
License as published by the Free Software Foundation, either version 3
of the License, or (at your option) any later version.

This code is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY. We would like to kindly ask you to acknowledge the authors
and don't remove their names from the code.

You should have received copy of the GNU General Public License along
with this program. If not, see <http://www.gnu.org/licenses/>.
