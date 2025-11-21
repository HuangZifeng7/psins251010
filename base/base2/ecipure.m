function [avpi, avpn] = ecipure(imu, ap0, t0)
% ECI(Earth Centered Inertial)-frame pure INS Alogrithm.
%
% Prototype: [avpi, avpn] = ecipure(imu, ap0, t0)
% Inputs: imu - [x/y/z gyro, x/y/z acc, t]
%         ap0 - initial attitude&position [pitch;roll;yaw; lat;lon;hgt]
% Output: avpi - ECI-frame AVP =[pitch/roll/yaw in i-frame, vx,vy,vz, x,y,z, t]
%         avpn - ENU-frame AVP =[pitch,roll,yaw, vE,vN,vU, lat,lon,hgt, t];
%
% Examples:
%   avp0 = avpset([0;0;0]*glv.deg, 0, [34;116;0]*glv.deg, 0);
%   imu = imustatic(avp0, 0.01, 300);
%   [avpi, avpn] = ecipure(imu, avp0([1:3,7:9]));
%
%   load trj10ms.mat;
%   [avpi, avpn] = ecipure(trj.imu, trj.avp0([1:3,7:9]), 12*3600);
%   avpcmpplot(trj.avp, avpn); 
%   avp1=inspure(trj.imu,trj.avp0,'f'); avpcmpplot(trj.avp,avp1);
%
% See also  lcipure, lcefpure, inspure.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 07/05/2025
global glv
    if nargin<3, t0=0; end
    [nn,ts,nts,nts2] = nnts(2,diff(imu(1:2,end)));
    t0 = t0 + imu(1,end)-ts;
    Cie = rv2m([0;0;glv.wie*t0]);   [re, Cen] = blh2xyz(ap0(4:6));  Cnb = a2mat(ap0(1:3));
    qib = m2qua(Cie*Cen*Cnb);  vi = Cie*cross([0;0;glv.wie],re);  ri = Cie*re;
    len = length(imu);    avpi = zeros(fix(len/nn), 10);  avpn = avpi;
    ki = timebar(nn, len, 'ECI navigation processing.');
    for k=1:nn:len-nn+1
        k1 = k+nn-1;  wvm = imu(k:k1,1:6);  t = imu(k1,end);
        [phim, dvbm] = cnscl(wvm);
        Cie = rv2m([0;0;glv.wie*(t0+t-nts2)]);
        fe = gravj4(Cie'*(ri+vi*nts2),6);  % gravity calculate
        vi1 = vi + qmulv(qib,dvbm) + Cie*fe*nts;  % vel update
        ri = ri + (vi+vi1)*nts2;  vi=vi1;  % pos update
        qib = qupdt(qib, phim);  % att update
        avpi(ki,:) = [q2att(qib); vi; ri; t];  % ECI frame
        % output
        Cie = rv2m([0;0;glv.wie*(t0+t)]);
        re = Cie'*ri;  [pos,Cen] = xyz2blh(re);  Cni = (Cie*Cen)';
        avpn(ki,:) = [m2att(Cni*q2mat(qib)); Cni*vi-Cen'*cross([0;0;glv.wie],re); pos; t];
        ki = timebar;
    end
    insplot(avpn,'avp');


