function [llh, RN] = xyz2llhBatch(xyz)
% Convert ECEF Cartesian coordinate [x;y;z] to geographic
% coordinate [lat;lon;height].
%
% Prototype: [llh, RN] = xyz2llhBatch(xyz)
% Input: xyz - ECEF Cartesian coordinate vector, in meters 
% Outputs: llh - geographic coordinate blh=[lat;lon;height],
%               where lat & lon in radians and hegtht in meter
%          RN - radius of curvature in prime vertical
%
% Example
%   llh = [(0:90)'*glv.deg, randn(91,1), randn(91,1)*1000];
%   xyz = blh2xyzBatch(llh);   llh1 = xyz2llhBatch(xyz);
%   myfig, plot(llh(:,1)/glv.deg, llh(:,1:2)-llh1(:,1:2));
%   myfig, plot(llh(:,1)/glv.deg, llh(:,3)-llh1(:,3));
%
% See also xyz2llh, xyz2blh, blh2xyzBatch, pos2ceg.

% Copyright(c) 2009-2024, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 28/10/2024
global glv
    s = sqrt(xyz(:,1).^2+xyz(:,2).^2);
    idx = s>abs(xyz(:,3));  llh=xyz; RN=llh(:,1);
    xyz=llh(idx,:);
    if ~isempty(xyz)  % low latitude
        tanL = xyz(:,3)./s(idx);
        for k=1:10
            tanL0 = tanL;
            tanL = (xyz(:,3)+glv.Re*glv.e2*tanL0./sqrt(1+(1-glv.e2).*tanL0.^2))./s(idx);
            if max(abs(tanL-tanL0))<1e-16, break; end
        end
        lat = atan(tanL);  sl = sin(lat);
        RN(idx,1) = glv.Re./sqrt(1-glv.e2*sl.^2);
        h = s(idx)./cos(lat)-RN(idx,1);
        lon = atan2(xyz(:,2),xyz(:,1));
        llh(idx,1:3) = [lat,lon,h];
    end
    idx = ~idx;
    xyz=llh(idx,:);
    if ~isempty(xyz)  % high latitude
        cotL = s(idx)./xyz(:,3);
        for k=1:10
            cotL0 = cotL;
            cotL = s(idx)./(xyz(:,3)+sign(cotL0).*glv.Re*glv.e2./sqrt(cotL0.^2+1-glv.e2));
            if max(abs(cotL-cotL0))<1e-16, break; end
        end
        lat = acot(cotL);  sl = sin(lat);
        RN(idx,1) = glv.Re./sqrt(1-glv.e2*sl.^2);
        h = xyz(:,3)./sl-RN(idx,1)*(1-glv.e2);
        lon = atan2(xyz(:,2),xyz(:,1));
        llh(idx,1:3) = [lat,lon,h];
    end
    