function [CeG, slat, clat, slon, clon, zeta] = pos2ceg(pos)
% Convert geographic pos = [lat; lon; *] to transformation matrix CeG (
% from Earth-frame to grid-frame ).
%
% Prototype: [CeG, slat, clat, slon, clon, zeta] = pos2ceg(pos)
% Input: pos - geographic position [lat,lon,hgt]
% Output: CeG - transformation matrix from Earth-frame to grid-frame
%         slat, clat, slon, clon, zeta - see the code
%
% See also  pos2cng, pos2cen, pos2cew, cen2pos, blh2xyz, xyz2blh, nv2cen.

% Copyright(c) 2009-2024, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 28/11/2024
    slat = sin(pos(1)); clat = cos(pos(1));  slon = sin(pos(2)); clon = cos(pos(2));
    gs = clat*slon; gc = clat*clon;  zeta = 1/sqrt(1-gs*gs);
    CeG = [ -zeta*gs*gc,  -zeta*slat, gc
             1/zeta,       0,         gs
            -zeta*gs*slat, zeta*gc,   slat ];
