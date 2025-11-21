function avpi = lcipure(imu, ap0, A0)
% LCI-frame pure INS for launch vehicle.
%
% Prototype: avpi = lcipure(imu, ap0, A0)
% Inputs: imu - [x/y/z gyro, x/y/z acc, t], x-forward, y-upward, z-rightward
%         ap0 - initial launch attitude&position [pitch;roll;yaw; lat;lon;hgt]
%                          where pitch;roll;yaw in pitch-yaw-roll sequence
%         A0 - firing/shooting angle, NOTE: counter-clockwise to be +
% Output: avpi - LCI-frame AVP =[pitch,roll,yaw, vx,vy,vz, x,y,z, t]
%                          FUR pitch-yaw-roll for att, 
%                          x-forward, y-upward, z-rightward for vel&pos
%
% See also  lcefpure, ecipure, lciavp2imu, lciavp2avp, atttrans, q2att1, inspure.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 01/05/2025
global glv
    if size(ap0,1)>7, ap0=ap0([1:3,7:9]); end
    [nn,ts,nts,nts2] = nnts(2,diff(imu(1:2,end)));
    Cn0a = a2mat([0;A0;0]); % [F;U;R]
    [re0,Ce0n0] = blh2xyz(ap0(4:6));  Ce0n0 = CenNUE(Ce0n0);  Ce0a = Ce0n0*Cn0a;  ra0 = Ce0a'*re0;
    qab = a2qua1(ap0(1:3));  % ^a for launch inertial frame, _b for body frame
    va = Ce0a'*cross([0;0;glv.wie],re0);  pa = zeros(3,1);
    len = length(imu);    avpi = zeros(fix(len/nn), 10);
    ki = timebar(nn, len, 'Launch vehicle navigation processing.');
    for k=1:nn:len-nn+1
        k1 = k+nn-1;  wvm = imu(k:k1,1:6);  t = imu(k1,end);
        [phim, dvbm] = cnscl(wvm);
        Ceka = rv2m([0;0;-glv.wie*(t-nts2)])*Ce0a;  % C^ek_e0*C^e0_a
        re = Ceka*(ra0+pa+va*nts2);
        fe = gravj4(re);  % gravity calculate
        va1 = va + qmulv(qab,dvbm) + Ceka'*fe*nts;  % vel update
        pa = pa + (va+va1)*nts2;  va=va1;  % pos update
        qab = qupdt(qab, phim);  % att update
        avpi(ki,:) = [q2att1(qab); va; pa; t];  % RFU-LCI launch frame
        ki = timebar;
    end
    insplot(avpi,'avpi'); subplot(222), title(sprintf('\\psi_0 = %.3f \\circ',A0/glv.deg));


