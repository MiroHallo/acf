% RUN_EXAMPLE Plots covariance matrices from Fortran program run_example
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This plots covariance matrices produced by Fortran program run_example
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Miroslav HALLO
% Charles University in Prague, Faculty of Mathematics and Physics
% E-mail: hallo@karel.troja.mff.cuni.cz
% Revision 2026/04: First version
% Tested in Matlab R2025b
% Method:
% Hallo, M., Gallovic, F. (2016): Fast and cheap approximation of Green 
%      functions uncertainty for waveform-based earthquake source inversions,
%      Geophys. J. Int., 207 1012-1029. https://doi.org/10.1093/gji/ggw320
%
% Copyright (C) 2026 Miroslav Hallo
%
% This program is published under the GNU General Public License (GNU GPL).
%
% This program is free software: you can modify it and/or redistribute it
% or any derivative version under the terms of the GNU General Public
% License as published by the Free Software Foundation, either version 3
% of the License, or (at your option) any later version.
%
% This code is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY. We would like to kindly ask you to acknowledge the authors
% and don't remove their names from the code.
%
% You should have received copy of the GNU General Public License along
% with this program. If not, see <http://www.gnu.org/licenses/>.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INIT:
close all;
clearvars;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:

N = 200;      % Length of data [samples]
dt = 0.1;     % Sample interval [s]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prepare time
t = (0 : N-1) * dt;

% Read covariance matrices from Fortran program run_example
Caa = readmatrix('example_ACF.txt');
Cxa = readmatrix('example_AXCF.txt');
Cas = readmatrix('example_SACF.txt');
Cxs = readmatrix('example_SAXCF.txt');

% Plot covariance matrix by ACF
fig = figure('Color','white');
subplot(2,2,1)
imagesc(t,t,Caa);
axis image;
colorbar;
xlabel('Time (s)'); ylabel('Time (s)')
title('Cov. matrix by ACF', 'FontWeight', 'normal')

% Plot cross-covariance matrix by AXCF
subplot(2,2,2)
imagesc(t,t,Cxa);
axis image;
colorbar;
xlabel('Time (s)'); ylabel('Time (s)')
title('Cross-cov. matrix by AXCF', 'FontWeight', 'normal')

% Plot covariance matrix by SACF
subplot(2,2,3)
imagesc(t,t,Cas);
axis image;
colorbar;
xlabel('Time (s)'); ylabel('Time (s)')
title('Cov. matrix by SACF', 'FontWeight', 'normal')

% Plot covariance matrix by SAXCF
subplot(2,2,4)
imagesc(t,t,Cxs);
axis image;
colorbar;
xlabel('Time (s)'); ylabel('Time (s)')
title('Cross-cov. matrix by SAXCF', 'FontWeight', 'normal')

% Save figure
exportgraphics(fig, 'example_Fortran.png', 'Resolution', 300)

