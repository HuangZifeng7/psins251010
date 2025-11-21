function [tpos, Ctn] = pos2tpos(pos)
% Convert pos=[lat;lon;hgt] to transverse tpos=[T_lat;T_lon;hgt].
%  NOTE: XYZ-ECEF -> YZX-transverse ECEF
%
% Prototype: tpos = pos2tpos(pos)
% Input: pos - =[lat;lon;hgt] in ENU frame
% Outputs: tpos - =[lat;lon;hgt] in transverse ENU frame
%          Ctn - rotation from transverse ENU frame to ENU frame
%
% Example
%    pos0=[0;pi;1];
%    [tpos, Ctn]=pos2tpos(pos0); [pos1, Cnt]=tpos2pos(tpos);
%    pos1-pos0, Ctn-Cnt'
%
% See also  tpos2pos, cen2pos, pos2ceg, blh2xyz, xyz2blh, a2mat, pp2vn.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 05/05/2025
    slat = sin(pos(1)); clat = cos(pos(1)); 
    slon = sin(pos(2)); clon = cos(pos(2));
    cs = clat*slon;  sqcs = sqrt(1-cs*cs);
    tpos = [atan(cs/sqcs); atan2(clat*clon,slat); pos(3)];
    if nargout>1
        cafa = -slat*slon/sqcs;  safa = clon/sqcs;
        Ctn = [cafa,-safa,0; safa,cafa,0; 0,0,1];
    end