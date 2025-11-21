% Particle Kalman filter(PF) simulation with large misalignment angles.
% See also  test_align_some_methods, test_align_ekf, test_align_ukf.
% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 12/02/2025
glvs;
psinstypedef('test_align_pf_def');
[nn, ts, nts] = nnts(4, 0.1);
T = 600;
avp0 = avpset([10;0;-180], [0;0;0], [34;108;380]); qnb0 = a2qua(avp0(1:3)');
imuerr = imuerrset(0.01, 50, 0.001, 5);
[imu, eth] = imustatic(avp0, ts, T, imuerr);   % imu simulation
afa = [-5; 7; 80]/100*glv.deg;  % large misalignment angles
qpb = qaddafa(qnb0,afa);  vn = zeros(3,1);
kf = kfinit(nts, imuerr, 500); kf.s = 1.01; % forgetting factor
kf.particle(1:3,:) = kf.particle(1:3,:)+afa;
% kf.alpha = 1; kf.beta = 0; kf.kappa = 0;
len = length(imu); [res, xkpk] = prealloc(fix(len/nn), 6, 2*kf.n+1);
ki = timebar(nn, len, 'PF align simulation.');
for k=1:nn:len-nn+1
    k1 = k+nn-1;
    wvm = imu(k:k1,1:6);  t = imu(k1,7);
    [phim, dvbm] = cnscl(wvm);
    dvn = qmulv(qpb,dvbm); vn = vn + dvn + eth.gn*nts;
    qpb = qupdt(qpb, phim-qmulv(qconj(qpb),eth.wnie*nts));
    res(ki,:) = [qq2afa(qpb, qnb0); vn]';
    kf.px = [eth.wnie; dvn/nts; nts];
    for kk=1:kf.Npar
        kf.particle(:,kk) = afamodel(kf.particle(:,kk), kf.px) + kf.sQ.*randn(6,1);
        Verr = vn - kf.particle(4:6,kk);
        kf.weight(kk) = exp(-0.5*Verr'*(kf.Rk\Verr));
    end
    kf.weight = kf.weight+eps;  kf.weight = kf.weight/sum(kf.weight);
    kf.particle = kf.particle(:,randsample(kf.Npar,kf.Npar,true,kf.weight));
    kf.xk = mean(kf.particle,2);  kf.Pxk = cov(kf.particle');
    xkpk(ki,:) = [kf.xk; diag(kf.Pxk); t]';
    ki = timebar;
end
kfplot(xkpk, res);
