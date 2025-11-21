function [imu, avp0] = avp2imu_ecef(avp)
% Simulate SIMU sensor outputs from attitude, velocity & position profile.
%
% Prototype: [imu, avp0] = avp2imu_ecef(avp)
% Input: avp = [attECEF,ve,pe,t]
% Outputs: imu = [wm,vm,t]
%          avp0 = init [att,vn,pos] in ENU-frame
%
% See also  ap2avp_ecef, imupolar_ecef, ap2avp, avp2imu, trajsimu, insupdate.

% Copyright(c) 2009-2024, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 28/11/2024
global glv
    if size(avp,2)<9, avp = ap2avp(avp, diff(avp(1:2,end))); end
    len = size(avp,1);  ts = avp(2,10)-avp(1,10);  ts2 = ts/2;
    Cbe_1 = a2mat(avp(1,1:3)')';  ve_1 = avp(1,4:6)';  pe_1 = avp(1,7:9)';
    wm_1 = [0;0;0];  vm_1 = [0;0;0];  weie = [0;0;glv.wie];
    imu = zeros(len, 6);
    timebar(1, len, 'Trajectory inversion avp-ecef->imu.');
    for k=2:len  % begin from 2
        Ceb = a2mat(avp(k,1:3));   ve = avp(k,4:6)';   pe = avp(k,7:9)';
        pe01 = (pe_1+pe)/2;
        [llh, RN, sl] = xyz2llh(pe01);
        g = glv.g0*(1+5.2790414e-3*sl^2+2.32718e-5*sl^4)-3.086e-6*llh(3);
        ge01 = -g*[pe01(1:2)/(RN+llh(3)); sl];
        gcc = -2*cross(weie,(ve_1+ve)/2) + ge01;
        phim = m2rv(Cbe_1*rv2m(weie*ts)*Ceb);
        wm = (glv.I33+askew(wm_1/12))^-1*phim; % using previous subsample: phim = wm + 1/12*cross(wm_1,wm)
        dvbm = Cbe_1*qmulv(rv2q(weie*ts2), ve-ve_1-gcc*ts); % sins
        vm = (glv.I33+askew(wm/2+wm_1/12))^-1*(dvbm-cros(vm_1,wm)/12);  % dvbm = vm + 1/2*cross(wm,vm)
        imu(k,:) = [wm; vm]';
        Cbe_1 = Ceb'; ve_1 = ve; pe_1 = pe; wm_1 = wm; vm_1 = vm;
        timebar;
    end
    imu = [imu(2:end,:), avp(2:end,10)];
    avp0 = avptrans(avp(1,1:9)','e2n');

