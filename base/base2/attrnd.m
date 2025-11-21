function att = attrnd(pr0)
% Set random attitude = [pitch; roll; yaw], where
% pitch within [-pr0, pr0], roll within [-pr0, pr0], yaw within [-pi, pi].
%
% Prototype: att = attrnd(pr0)
% Input: pr0 - random pitch/roll range
% Output: att - random attitude
% 
% See also  posrnd, a2qua.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 12/05/2025
    if nargin<1, pr0=10*pi/180; end
    att = [rand(2,1)*2*pr0-pr0; rand*2*pi-pi];