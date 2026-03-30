# Uncertainty quantification in earthquake source inversions
Approximate covariance matrix of Green's functions for uncertainty quantification in earthquake source inversions
***************************************

  Open-source functions for determining the (cross-)covariance matrix
of Green's functions by approximate covariance function (ACF, AXCF) and
stationarized approximate covariance function (SACF, SAXCF). Functions
are distributed with an intuitive example, and return the full
(cross-)covariance matrices.

1 METHODOLOGY
===================

  Hallo, M., Gallovič, F. (2016). Fast and cheap approximation of Green
functions uncertainty for waveform-based earthquake source inversions,
Geophys. J. Int., 207, 1012-1029. [https://doi.org/10.1093/gji/ggw320](https://doi.org/10.1093/gji/ggw320)

2 TECHNICAL IMPLEMENTATION
===================

Fourier Transform (complex numbers), Cross-Platform (Windows, Linux)

3 PACKAGE CONTENT
===================

  a) "approxc.f90" - Fortran subroutines for determining the
  (cross-)covariance matrix by (stationarized) approximate covariance
  function (AXCF, SAXCF)
  
  b) "axcf.m" - Matlab function for the (cross-)covariance matrix by
  approximate covariance function (ACF, AXCF)
  
  c) "saxcf.m" - Matlab function for the (cross-)covariance matrix by
  stationarized covariance function (SACF, SAXCF)
  
  d) "example.m" - Example of Matlab code using these functions

4 REQUIREMENTS
===================

  FORTRAN: Codes should fulfill Fortran 90 Standard. The codes were
successfully compiled by GFortran (GCC) and ifort (Intel) compilers on
Ubuntu 14 operation system.

  MATLAB: Version R2016b, "smooth" function from Matlab Curve Fitting Toolbox, and
"filtfilt" function from Matlab Signal Processing Toolbox.

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
