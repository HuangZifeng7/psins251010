function [pos, Cnt] = tpos2pos(tpos)
% Convert transverse tpos=[T_lat;T_lon;hgt] to pos=[lat;lon;hgt].
%
% Prototype: tpos = pos2tpos(pos)
% Input: tpos - =[lat;lon;hgt] in transverse ENU frame
% Outputs: pos - =[lat;lon;hgt] in ENU frame
%          Cnt - rotation from ENU frame to transverse ENU frame
%
% See also  pos2tpos, cen2pos, pos2ceg, blh2xyz, xyz2blh, a2mat, pp2vn.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 05/05/2025
    slat = sin(tpos(1)); clat = cos(tpos(1)); 
    slon = sin(tpos(2)); clon = cos(tpos(2));
    cc = clat*clon;  sqcc = sqrt(1-cc^2);
    pos = [atan(cc/sqcc); atan2(slat,clat*slon); tpos(3)];
    if nargout>1
        cafa = -slat*clon/sqcc;  safa = slon/sqcc;
        Cnt = [cafa,safa,0; -safa,cafa,0; 0,0,1];
    end