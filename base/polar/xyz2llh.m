function [llh, RN, sl] = xyz2llh(xyz)
% Convert ECEF Cartesian coordinate [x;y;z] to geographic
% coordinate [lat;lon;height].
%
% Prototype: [llh, RN, sl] = xyz2llh(xyz)
% Input: xyz - ECEF Cartesian coordinate vector, in meters 
% Outputs: llh - geographic coordinate blh=[lat;lon;height],
%               where lat & lon in radians and hegtht in meter
%          RN, sl - radius of curvature in prime vertical, sin(lat)
%
%  Example
%    xyz = blh2xyz([-89.9*glv.deg; 10*glv.deg; 12345]);
%    llh = xyz2llh(xyz); llh1 = xyz2blh(xyz); llh1-llh
%
% See also xyz2blh, blh2xyz, pos2ceg

% Copyright(c) 2009-2024, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 28/10/2024
global glv
    s = sqrt(xyz(1)^2+xyz(2)^2);
    if s>abs(xyz(3))  % low latitude
        tanL = xyz(3)/s;
        for k=1:10
            tanL0 = tanL;
            tanL = (xyz(3)+glv.Re*glv.e2*tanL0/sqrt(1+(1-glv.e2)*tanL0^2))/s;
            if abs(tanL-tanL0)<1e-16, break; end
        end
        lat = atan(tanL);  sl = sin(lat);
        RN = glv.Re/sqrt(1-glv.e2*sl^2);
        h = s/cos(lat)-RN;
    else  % high latitude
        cotL = s/xyz(3);
        for k=1:10
            cotL0 = cotL;
            cotL = s/(xyz(3)+sign(cotL0)*glv.Re*glv.e2/sqrt(cotL0^2+1-glv.e2));
            if abs(cotL-cotL0)<1e-16, break; end
        end
        lat = acot(cotL);  sl = sin(lat);
        RN = glv.Re/sqrt(1-glv.e2*sl^2);
        h = xyz(3)/sl-RN*(1-glv.e2);
    end
% 	if s<1e-6, lon=0; else, lon=atan2(xyz(2),xyz(1)); end  % no need, atan2(0,0)=0
	lon = atan2(xyz(2),xyz(1));
    llh = [lat; lon; h];  % [llh(1:2)/glv.deg; h]',
    