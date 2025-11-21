function Cew = pos2cew(pos, afa)
% Convert geographic pos = [lat; lon; *] & wander azimuth to wander frame
%   transformation matrix Cew ( from Earth-frame to wander-frame ).
%
% Prototype: Cew = pos2cew(pos, afa)
% Inputs: pos - geographic position
%         afa - wander azimuth
% Output: Cew - transformation matrix from Earth-frame to wander-frame
%
% See also  pos2cen, pos2ceg, cew2pos.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 25/03/2025
    if nargin<2, afa=0; end
    Cew = pos2cen(pos)*rxyz(-afa,'z');
