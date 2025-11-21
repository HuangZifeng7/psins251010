glvs;

ts = 1/100;

load mpu6500gsp.mat; % imuplot(imu); gpsplot(gps);

t0 = 2;

avp0 = [[0;0;-100]*glv.deg; 0;0;0; getat(gps,t0)];

ins = insinit(avp0, ts);

avperr = avperrset([10*60;30*60], 10, 100);

imuerr = imuerrset(500, 5000, 5, 500);

Pmin = [avperrset([0.5,3],0.1,0.1); gabias(1.0, [100,100]); [0.01;0.01;0.01]; 0.001].^2;

Rmin = poserrset(0.1).^2;

[avp1, xkpk, zkrk, sk] = sinsgps(imu, gps, ins, avperr, imuerr, rep3(1), 0.1, poserrset(10), Pmin, Rmin, 'avped');

avpcmpplot(gps, adddt(avp1(:,[7:9,end]),0.0), 'p');