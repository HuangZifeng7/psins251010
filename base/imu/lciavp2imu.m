function [imu, avp] = lciavp2imu(avp, pos0, A0)
% Transform LCI AVP to FUR-IMU for launch vehicle use.
%
% Prototype: [imu, avp] = lciavp2imu(avp, pos0, A0)
% Inputs: avp - LCI-frame AVP =[pitch,roll,yaw, vx,vy,vz, x,y,z, t],
%                               FUR pitch-yaw-roll for att,
%                               x-forward, y-upward, z-rightward for vel&pos
%         pos0 - initial launch position [lat;lon;hgt]
%         A0 - firing/shooting angle, NOTE: counter-clockwise to be +
% Output: imu - IMU out =[x/y/z_gyro, x/y/z_acc, t]
%                         x-forward, y-upward, z-rightward for vel&pos
%
% See also  lciap2avp, avp2imu, lcipure, atttrans, q2att1, lciavp2avp.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 01/05/2025
global glv
    if avp(1,end)>0, avp=[[2*avp(1,1:9)-avp(2,1:9),0]; avp];  end  % make sure from t=0
    Cn0a = a2mat([0;A0;0]);
    [re0,Ce0n0] = blh2xyz(pos0); Ce0n0 = CenNUE(Ce0n0); Ce0a = Ce0n0*Cn0a;  ra0 = Ce0a'*re0;
    qab0 = a2qua1(avp(1,1:3)');  va0 = avp(1,4:6)';  pa0 = avp(1,7:9)';  t0 = avp(1,end);
    len = length(avp);  imu = zeros(len-1,7);
    ki = timebar(1, len, 'Translate LCI AVP to IMU.');
    for k=2:len
        qab = a2qua1(avp(k,1:3)');  va = avp(k,4:6)';  pa = avp(k,7:9)';  t = avp(k,end);
        ts = t-t0;  ts2 = ts/2;  t0 = t;
        Ceka = rv2m([0;0;-glv.wie*(t-ts2)])*Ce0a;  % C^ek_e0*C^e0_a
        re = Ceka*(ra0+(pa0+pa)/2);  pa0 = pa;
        fe = gravj4(re);  % gravity calculate
        dvbm = qmulv(qconj(qab0), va-va0-Ceka'*fe*ts);  va0 = va;
        phim = qq2rv(qab, qab0);  qab0 = qab;
        if k==2,  wm_1=phim; vm_1=dvbm;  end
        wm = (glv.I33+askew(wm_1/12))^-1*phim;
        vm = (glv.I33+askew(wm/2+wm_1/12))^-1*(dvbm-cros(vm_1,wm)/12);
        imu(ki,:) = [wm; vm; t]';  wm_1 = wm;  vm_1 = vm;
        ki = timebar;
    end
