function pos = posrnd(hgt0, lat0)
% Set random position = [latitude; logititude; height], where
% lat within [-pi/2, pi/2], lon within [-pi, pi], hgt within [0, hgt].
%
% Prototype: pos = posrnd(hgt0, lat0)
% Inputs: hgt0 - random height range
%         lat0 - random latitude range
% Output: pos - random position
% 
% See also  attrnd, avpset, llh, avpset, poserrset.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 12/05/2025
    if nargin<2, lat0=83*pi/180; end
    if nargin<1, hgt0=0; end
    pos = [rand*2*lat0-lat0; rand*2*pi-pi; rand*hgt0];