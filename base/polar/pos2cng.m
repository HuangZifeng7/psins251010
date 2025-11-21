function [CnG, CeG, Cen] = pos2cng(pos)
% Convert geographic pos = [lat; lon; *] to transformation matrix CnG (
% from ENU-frame to grid-frame £©.
%
% Prototype: CnG = pos2cng(pos)
% Input: pos - geographic position [lat,lon,hgt]
% Outputs: CnG - transformation matrix from ENU-frame to grid-frame
%          CeG - transformation matrix from Earth-frame to grid-frame
%          Cen - transformation matrix from Earth-frame to ENU-frame
%
% See also  pos2ceg, pos2cen, cen2pos, blh2xyz, xyz2blh, nv2cen.

% Copyright(c) 2009-2024, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 28/11/2024
    slat = sin(pos(1)); clat = cos(pos(1));  slon = sin(pos(2)); clon = cos(pos(2));
    gs = clat*slon; gc = clat*clon;  zeta = 1/sqrt(1-gs*gs);
    ssigma = slon*slat*zeta; csigma = clon*zeta;
    CnG = [  csigma,  ssigma,  0
            -ssigma,  csigma,  0
             0,       0,       1 ];
    if nargout>1
        CeG = [ -zeta*gs*gc,  -zeta*slat, gc
                 1/zeta,       0,         gs
                -zeta*gs*slat, zeta*gc,   slat ];        
    end
    if nargout>2
        Cen = [ -slon,  -slat*clon,  clat*clon
                 clon,  -slat*slon,  clat*slon
                 0,      clat,       slat      ];        
    end    