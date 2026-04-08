% EXAMPLE Example using functions axcf.m (ACF, AXCF) and saxcf.m (SACF, SAXCF).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code creates two signals and then compute their covariance and
% cross-covariance matrices (ACF, AXCF, SACF, SAXCF). -> Routines for evaluating
% approximate covariance matrices of Green's functions, designed for uncertainty
% quantification in Bayesian earthquake source inversions.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Miroslav HALLO, František GALLOVIC
% Charles University in Prague, Faculty of Mathematics and Physics
% E-mail: hallo@karel.troja.mff.cuni.cz
% Revision 2016/01: First version
% Revision 2018/01: Corrected for real signals
% Revision 2026/04: Updated for new Matlab version
% Tested in Matlab R2025b
% Method:
% Hallo, M., Gallovic, F. (2016): Fast and cheap approximation of Green 
%      functions uncertainty for waveform-based earthquake source inversions,
%      Geophys. J. Int., 207 1012-1029. https://doi.org/10.1093/gji/ggw320
%
% Copyright (C) 2016-2018 Miroslav Hallo and František Gallovič
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
% REQUIREMENTS:
%
% smooth (MatLab Curve Fitting Toolbox)
% filtfilt (MatLab Signal Processing Toolbox)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INIT:
close all;
clearvars;
assert(~isempty(which('smooth')), 'Curve Fitting Toolbox is required.')
assert(~isempty(which('filtfilt')), 'Signal Processing Toolbox is required.')
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:

N = 200;      % Length of synthetic data [samples]
dt = 0.1;     % Sample interval [s]
freq1 = 0.3;  % Frequency of the synthetic (cosine) signal [Hz]
T = 10;       % Duration of the synthetic signal [s]
gShift = 2.0; % Time shift of the generated signal g relatively to f [s]
L1 = 2.0;     % Width of joint (f and g) uniform distribution of time-shifts [s]
L12 = 1.0;    % Width of relative (f vs g) uniform distribution of time-shifts [s]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prepare synthetic signals f and g
t = (0 : N-1) * dt;
tau = t - t(round(N/2));
f = cos(tau*freq1*2*pi) .* exp(-tau.^2/T);
g = circshift(f',round(gShift/dt))';

% Covariance matrix by ACF
Caa = axcf(f,L1,dt);

% Cross-covariance matrix by AXCF
Cxa = axcf(f,g,L1,L12,dt);

% Covariance matrix by SACF
Cas = saxcf(f,L1,dt,T);

% Cross-covariance matrix by SAXCF
Cxs = saxcf(f,g,L1,L12,dt,T);

% Plot f signal
fig = figure('Color','white','Units','normalized','Position',[0.3 0.1 0.4 0.8]);
subplot(3,2,1)
plot(t,f,'LineWidth',1.2)
xlabel('Time (s)')
ylabel('Amplitude')
title('Signal f', 'FontWeight', 'normal')

% Plot f and g signals
subplot(3,2,2)
plot(t,f,t,g,'LineWidth',1.2)
xlabel('Time (s)')
ylabel('Amplitude')
legend('f', 'g')
title('Signals f and g', 'FontWeight', 'normal')

% Plot covariance matrix by ACF
subplot(3,2,3)
imagesc(t,t,Caa);
axis image;
colorbar;
xlabel('Time (s)'); ylabel('Time (s)')
title('Cov. matrix by ACF', 'FontWeight', 'normal')

% Plot cross-covariance matrix by AXCF
subplot(3,2,4)
imagesc(t,t,Cxa);
axis image;
colorbar;
xlabel('Time (s)'); ylabel('Time (s)')
title('Cross-cov. matrix by AXCF', 'FontWeight', 'normal')

% Plot covariance matrix by SACF
subplot(3,2,5)
imagesc(t,t,Cas);
axis image;
colorbar;
xlabel('Time (s)'); ylabel('Time (s)')
title('Cov. matrix by SACF', 'FontWeight', 'normal')

% Plot covariance matrix by SAXCF
subplot(3,2,6)
imagesc(t,t,Cxs);
axis image;
colorbar;
xlabel('Time (s)'); ylabel('Time (s)')
title('Cross-cov. matrix by SAXCF', 'FontWeight', 'normal')

% Save figure
exportgraphics(gcf, 'example.png', 'Resolution', 300)

