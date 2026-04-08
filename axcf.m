function C = axcf(varargin)
% AXCF Returns covariance matrix by Approximate Covariance Functions (ACF or AXCF)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Returns Auto-Covariance matrix by Approximate Covariance Function (ACF) or  
% Cross-covariance matrix by Approximate Cross-covariance Function (AXCF)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Miroslav HALLO, František GALLOVIC
% Charles University in Prague, Faculty of Mathematics and Physics
% E-mail: hallo@karel.troja.mff.cuni.cz
% Revision 2016/01: First version
% Revision 2018/01: Corrected by inserting zeros before and after the signal
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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:
% f - First input signal (vector)
% g - Second input signal (vector)
% L1 - Width of the joint uniform distribution of time-shifts [s] (scalar)
% L12 - Width of the relative uniform distribution of time-shifts [s] (scalar)
% dt - Signal sampling [s] (scalar)
% 
% OUTPUT:
% C - (Cross-)covariance matrix (2D matrix)
%
% EXAMPLES:
% axcf(f, L1, dt) - Compute auto-covariance matrix from f
% axcf(f, g, L1, L12, dt) - Compute cross-covariance matrix from f and g
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Manage input arguments
if nargin == 5   % cross-covariance
    f = varargin{1};
    g = varargin{2};
    L1 = varargin{3};
    L12 = varargin{4};
    dt = varargin{5};
else  
    if nargin == 3   % auto-covariance
        f = varargin{1};
        g = f;
        L1 = varargin{2};
        L12 = 0;
        dt = varargin{3};
    else
        error('axcf: Incorrect number of input arguments');
    end
end

% The number of samples in f
nsampl = length(f);

% Check the size of the final covariance matrix
if nsampl > 22360
    disp('axcf: Covariance matrix require more than 4GB of memory');
elseif nsampl > 11180
    disp('axcf: Covariance matrix require more than 1GB of memory');
elseif nsampl > 5590
    disp('axcf: Covariance matrix require more than 250MB of memory');
end

% Check the length and preallocate output C matrix
if nsampl ~= length(g)
    error('axcf: Signals f ang g must have the same length');
else
    C = zeros(nsampl);
end

% Prepare smoothing windovs
L1s = ceil(L1/dt);
if L1s < 3
    L1s = 3;
    disp(['axcf: L1 changed to the minimum allowed value: ',num2str(L1s*dt),' (s)'])
end
L12s = max(ceil(L12/dt),1);

% Put zeros (values) before and after the signal
Lzeros = max(L1s,L12s);
f = [zeros(Lzeros,1)+f(1); f(:); zeros(Lzeros,1)+f(end)];
g = [zeros(Lzeros,1)+g(1); g(:); zeros(Lzeros,1)+g(end)];

% New number of samples
nsamplN = nsampl + 2*Lzeros;

% Smooth f and g signals
fSmooth = smooth(f, L1s, 'moving');
gSmooth = smooth(g, L12s, 'moving');

% Compute ACF for time-lags
ACF = zeros(nsamplN,nsamplN*2);
for tshift = 1 : nsamplN*2 % loop for time-lags (nsamplN+1 is zero time-lag)
    gShift = circshift(gSmooth,tshift-(nsamplN+1));
    ACF(:,tshift) = smooth( f.*gShift, L1s, 'moving' ) ...
        - (fSmooth.*smooth( gShift, L1s, 'moving' )  );
end

% Fill the covariance matrix by ACF
for i=1:nsampl % loop for rows
    for j=1:nsampl % loop for columns
        C(i,j) = ACF(i+Lzeros,nsamplN+1-(j-i));
    end
end

% Check if covariance matrix is symmetric (in case any numerical issues)
if isequal(f,g) % only auto-covariance
    issym = @(x) all(all(x==x.'));
    if ~issym(C)
        % Symmetrize it
        Cm = triu(ones(nsampl));
        C = C.*Cm;
        C = C + tril(C.',-1);
        % disp('axcf: The covariance matrix was symmetrized')
    end
end

end
