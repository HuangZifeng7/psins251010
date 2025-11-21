function [Cen, xyz] = tpos2cen(pos)
% Convert transverse geographic pos = [lat; lon; *] to transformation
%  matrix Cen ( from Earth-frame to nav-frame ).
%
% Prototype: Cen = tpos2cen(pos)
% Input: pos - geographic position
% Outputs: Cen - transformation matrix from Earth-frame to transverse nav-frame
%          xyz - xyz in ECEF frame
%
% See also  tpos2pos, pos2cen, cen2pos, blh2xyz, xyz2blh, a2mat, pp2vn.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 05/02/2025
global glv
    if size(pos,1)>3  % batch process
        [slat, clat, slon, clon, ss, sc, cs, cc] = bl2sc(pos);
        Cen = [ clon,  -ss,   cs, ...
                ss*0, clat, slat, ...
               -slon,  -sc,   cc  ];
        if nargout>1
            RN = glv.Re./sqrt(1-glv.e2*cc.^2); RNh=RN+pos(:,3);
            xyz = [RNh.*clat.*slon, RNh.*slat, (RN*(1-glv.e2)+pos(:,3)).*cc];
        end
        return;
    end
    if length(pos)<2, pos = [pos(1);0;0]; end  % tpos2cen(lat)
    [slat, clat, slon, clon, ss, sc, cs, cc] = bl2sc(pos);
    Cen = [ clon,  -ss,   cs
               0, clat, slat
           -slon,  -sc,   cc ];
    if nargout>1
        RN = glv.Re/sqrt(1-glv.e2*cc^2);
        xyz = [(RN+pos(3))*[clat*slon;slat]; (RN*(1-glv.e2)+pos(3))*cc];
    end

% syms slon clon slat clat
% Cne=[-slon,clon,0; -slat*clon,-slat*slon,clat; clat*clon,clat*slon,slat]*[0,0,1; 1,0,0; 0,1,0];
% Cen=Cne.'