% Potter square-root Kalman filter simulation for SINS/GPS, 
% compared with single & doulbe date type.
% Please run 'test_SINS_trj.m' to generate 'trj10ms.mat' beforehand!!!
% See also  test_SINS_trj, test_SINS_GPS_153.
% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 12/04/2025
glvs
psinstypedef(153);
trj = trjfile('trj10ms.mat');  % insplot(trj.avp);
% initial settings
[nn, ts, nts] = nnts(1, trj.ts);
imuerr = imuerrset(0.1, 200, 0.01, 10);
%% 
imu = imuadderr(trj.imu, imuerr);
davp0 = avperrset([1;-1;30], 0.1, [5;5;5]);
ins = insinit(avpadderr(trj.avp0,davp0), ts);
% KF filter
rk = poserrset([5;5;5]);
kf = kfinit(ins, davp0, imuerr, rk);
% kf.Pmin = [avperrset(0.01,1e-4,0.1); gabias(1e-3, [1,10])].^2;  kf.pconstrain=1;
len = length(imu); [avp, xkpk] = prealloc(fix(len/nn), 10, 2*kf.n+1);  sxkpk=xkpk;  pxkpk=xkpk;
ki = 1;
sflag = 0;  % single or double datatype flag selection, =1 for single, =0 for double
if sflag==1, Qk=single(diag(kf.Qk(1:6,1:6))); Rk=single(diag(kf.Rk)); Pk=single(kf.Pxk);  Xk=single(kf.xk); % CommonKF setting
    timebar(nn, len, 'SINS/GPS Potter-SRKF(single datatype) Simulation.'); 
else, Qk=diag(kf.Qk(1:6,1:6)); Rk=diag(kf.Rk); Pk=kf.Pxk;  Xk=kf.xk;
    timebar(nn, len, 'SINS/GPS Potter-SRKF(double datatype) Simulation.');  end
sQk = sqrt(Qk); sRk = sqrt(Rk); Dk = chol(Pk);  pXk = Xk;  % PotterKF setting
for k=1:nn:len-nn+1
    k1 = k+nn-1;  
    wvm = imu(k:k1,1:6);  t = imu(k1,end);
    ins = insupdate(ins, wvm);
    kf.Phikk_1 = kffk(ins);  % double datatype as reference
    kf = kfupdate(kf);
    if sflag==1, Phikk_1=single(kf.Phikk_1); else, Phikk_1=kf.Phikk_1; end
    [Xk, Pk] = CommonKF(Phikk_1, [], Qk, Xk, Pk);  % common KF
    [pXk, Dk] = PotterKF(Phikk_1, [], sQk, pXk, Dk);  % Potter SRKF
    if mod(t,1)==0
        posGPS = trj.avp(k1,7:9)' + davp0(7:9).*randn(3,1);  % GPS pos simulation with some white noise
        Zk = ins.pos-posGPS;
        kf = kfupdate(kf, Zk, 'M');
        if sflag==1, Zk=single(Zk); Hk=single(kf.Hk); else, Hk=kf.Hk; end
        [Xk, Pk] = CommonKF([], [], [], Xk, Pk, Hk, Rk, Zk, 'M');
        [pXk, Dk, pPk] = PotterKF([], [], [], pXk, Dk, Hk, sRk, Zk, 'M');
        sxkpk(ki,:) = [Xk; diag(Pk); t]';
        pxkpk(ki,:) = [pXk; diag(pPk); t]';
        %[kf, ins] = kffeedback(kf, ins, 1, 'avp');
        avp(ki,:) = [ins.avp', t];
        xkpk(ki,:) = [kf.xk; diag(kf.Pxk); t]';  ki = ki+1;
    end
    timebar;
end
avp(ki:end,:) = [];  xkpk(ki:end,:) = [];   sxkpk(ki:end,:) = [];   pxkpk(ki:end,:) = []; 
% insplot(avp);
% avperr = avpcmpplot(trj.avp, avp);
% kfplot(xkpk, avperr, imuerr);
% kfplot(sxkpk, avperr, imuerr);
% kfplot(pxkpk, avperr, imuerr);
% avpcmpplot(xkpk(:,[1:9,end]),sxkpk(:,[1:9,end]),pxkpk(:,[1:9,end]));
avpcmpplot(xkpk(:,[1:9,end]),sxkpk(:,[1:9,end]));
avpcmpplot(xkpk(:,[1:9,end]),pxkpk(:,[1:9,end])); 
% avpcmpplot(xkpk(:,[16:24,end]),sxkpk(:,[16:24,end]),pxkpk(:,[16:24,end]));
avpcmpplot(xkpk(:,[16:24,end]),sxkpk(:,[16:24,end]));
avpcmpplot(xkpk(:,[16:24,end]),pxkpk(:,[16:24,end]));

