function [avp, avpn] = lcefpure(imu, ap0, A0)
% LCEF-frame (launch centered Earth-fixed) pure INS for launch vehicle.
%
% Prototype: avp = lcefpure(imu, ap0, A0)
% Inputs: imu - [x/y/z gyro, x/y/z acc, t], x-forward, y-upward, z-rightward
%         ap0 - initial launch attitude&position [pitch;roll;yaw; lat;lon;hgt]
%                          where pitch;roll;yaw in pitch-yaw-roll sequence
%         A0 - firing/shooting angle, NOTE: counter-clockwise to be +
% Output: avp - LCEF-frame AVP =[pitch,roll,yaw, vx,vy,vz, x,y,z, t]
%                          FUR pitch-yaw-roll for att, 
%                          x-forward, y-upward, z-rightward for vel&pos
%
% Example
%    avp0 = [0;0;0; 0;0;0; glv.pos0]; imu = imulci(imustatic(avp0, 0.01, 100));
%    [avp, avpn] = lcefpure(imu, avp0, 0);  % insplot(avpn);
%
%    load trj10ms.mat;  avp0=[trj.avp0(1:2);0;trj.avp0(7:9)];  inspure(trj.imu, trj.avp0);
%    [avp, avpn] = lcefpure(imulci(trj.imu), avp0, trj.avp0(3));  % insplot(avpn);
%
% See also  lcipure, ecipure, inspure, atttrans.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 01/05/2025
global glv
    if size(ap0,1)>7, ap0=ap0([1:3,7:9]); end
    [nn,ts,nts,nts2] = nnts(2,diff(imu(1:2,end)));
    Cnp0a = a2mat([0;A0;0]); % [F;U;R]  C^np0_a
    [re0, Cenp0] = blh2xyz(ap0(4:6));  Cenp0 = CenNUE(Cenp0);  Cea = Cenp0*Cnp0a;
    waie = Cea'*[0;0;1]*glv.wie;
    qab = a2qua1(ap0(1:3));  % ^a for LCEF frame, _b for body frame
    va = zeros(3,1);  pa = zeros(3,1);
    len = length(imu);   avp = zeros(fix(len/nn), 10);  avpn = avp;
    ki = timebar(nn, len, 'Launch vehicle (LCEF) navigation processing.');
    for k=1:nn:len-nn+1
        k1 = k+nn-1;  wvm = imu(k:k1,1:6);  t = imu(k1,end);
        [phim, dvbm] = cnscl(wvm);
        re = re0 + Cea*(pa+va*nts2);
        [fe, ge] = gravj4(re);  % gravity calculate
        va1 = va + rotv(-waie*nts2, qmulv(qab,dvbm)) + (-2*cros(waie,va) + Cea'*ge)*nts;  % vel update
        pa = pa + (va+va1)*nts2;  va=va1;  % pos update
        qab = qupdt2(qab, phim, waie*nts);  % att update
        avp(ki,:) = [q2att1(qab); va; pa; t];  % RFU-LCEF launch frame
        [pos, Cen] = xyz2blh(re0+Cea*pa);  Cen=CenNUE(Cen); Cna = Cen'*Cea;  % trans to ENU AVP
        Cnb = Cna*q2mat(qab);  vn = Cna*va;
        avpn(ki,:) = [m2att(Cnb([3,1,2],[3,1,2])); vn([3,1,2]); pos; t];
        ki = timebar;
    end
    insplot(avp,'avpi'); subplot(222), title(sprintf('\\psi_0 = %.3f \\circ',A0/glv.deg));


