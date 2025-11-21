% Trajectory generation for launch vehicle simulation use. Please run
% 'test_SINS_launch_vehicle_trj.m' to generate 'trj10ms.mat' beforehand!!!
% See also  test_SINS_launch_vehicle_trj, lcipure, lciavp2imu, test_SINS_trj.
% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 27/04/2025
glvs
load trj10ms.mat;  % insplot(trj.avp);
[nn,ts,nts,nts2] = nnts(2,trj.ts);
A0 = trj.avp0(3);  Cn0a = a2mat([0;0;A0]);
pos0 = trj.avp0(7:9);  [re0,Ce0n0] = blh2xyz(pos0);  Ce0a = Ce0n0*Cn0a;  ra0 = Ce0a'*re0;
qab = a2qua([trj.avp0(1:2);0]);  % ^a for launch inertial frame, _b for body frame
va = Ce0a'*cross([0;0;glv.wie],re0);  pa = zeros(3,1);
len = length(trj.imu);    avpa = zeros(fix(len/nn), 10);
ki = timebar(nn, len, 'Launch vehicle navigation processing.');
for k=1:nn:len-nn+1
    k1 = k+nn-1;  wvm = trj.imu(k:k1,1:6);  t = trj.imu(k1,end);
    [phim, dvbm] = cnscl(wvm);
    Ceka = rv2m([0;0;-glv.wie*(t-nts2)])*Ce0a;  % C^ek_e0*C^e0_a
    re = Ceka*(ra0+pa+va*nts2);
    fe = gravj4(re);  % gravity calculate
    va1 = va + qmulv(qab,dvbm) + Ceka'*fe*nts;  % vel update
    pa = pa + (va+va1)*nts2;  va=va1;  % pos update
    qab = qupdt(qab, phim);  % att update
    avpa(ki,:) = [q2att(qab); va; pa; t];  % RFU-LCI launch frame
    ki = timebar;
end
insplot(avpa,'avpi'); subplot(222), title(sprintf('\\psi_0 = %.3f \\circ',A0/glv.deg));
avpn = lciavp2avp(avpa, pos0, A0, 0);
avpcmpplot(trj.avp, avpn);
%% directly in launch FUR-frame
avpi = lcipure(trj.imu(:,[2,3,1,5,6,4,7]), [trj.avp0(1:2);0;pos0], A0);
avpn = lciavp2avp(avpi, pos0, A0, 1);
avpcmpplot(trj.avp, avpn);
%% PYR-AVP to FUR-IMU
avpi1 = lciap2avp(avpi, 0.01);
[imu,avpr] = lciavp2imu(avpi1,pos0,A0);
imuplot(imu,trj.imu(:,[2,3,1,5,6,4,7]));
avpi1 = lcipure(imu, [trj.avp0(1:2);0;pos0], A0); 
avpcmpplot(avpr, avpi1, 'avpi');
avpn1 = lciavp2avp(avpi1, pos0, A0, 1);
avpcmpplot(trj.avp, avpn1);

% avp1=inspure(trj.imu, trj.avp0, 'f');
% avpcmpplot(trj.avp, avp1);

