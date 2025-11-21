%https://www.psins.org.cn/newsinfo/4552096.html

glvs

ts = 1/100;

load stim300_cdgnss_static.mat;

imuplot(imu); gpsplot(gnss);

t0 = 1; t1 = 300;

att0 = alignsb(imu(1:10/ts,:),gnss(:,4:6)); att0(3)=0;

ins = insinit([att0;gnss(1,4:6)'], ts);

avperr = avperrset([60;300], 1, 10);

imuerr = imuerrset(100, [100;100;10000], 0.1, [10;10;100]);

Pmin = [avperrset([0.2,1.0],0.01,0.2); gabias(0.01, [10,10]); [0.01;0.01;0.01]; 0.001].^2;

[avp, xkpk, zkrk, sk, ins1, kf] = sinsgps(imu(t0/ts:t1/ts,:), gnss, ins, avperr, imuerr, rep3(0), 0, vperrset(0.01,0.01), Pmin, 0, 'avped');

avpcmpplot(gnss, avp(:,[4:9,end]), 'vp');