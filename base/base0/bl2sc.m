function [slat, clat, slon, clon, ss, sc, cs, cc] = bl2sc(bl)
% Calculate sin/cos for bl=[lat;lon;*]. Note: BLH for Breite/Laenge/Hoehe?
%
% Prototype: [slat, clat, slon, clon, ss, sc, cs, cc] = bl2sc(bl)
% Input: bl - geographic coordinate bl=[lat;lon;*],
%             where lat & lon in radians
% Outputs: slat, clat, slon, clon - sin/cos
%          ss, sc, cs, cc - sin*sin,sin*cos...
%
% See also  blh2xyz, xyz2blh, Dblh2Dxyz, pos2cen, blh2xyzBatch.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 06/02/2025
    if size(bl,2)>1  % batch process
        slat = sin(bl(:,1)); clat = cos(bl(:,1));
        slon = sin(bl(:,2)); clon = cos(bl(:,2));
        if nargout>4
            ss = slat.*slon; sc = slat.*clon;
            cs = clat.*slon; cc = clat.*clon;
        end
    else
        slat = sin(bl(1)); clat = cos(bl(1));
        slon = sin(bl(2)); clon = cos(bl(2));
        if nargout>4
            ss = slat*slon; sc = slat*clon;
            cs = clat*slon; cc = clat*clon;
        end
    end