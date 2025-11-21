glvs
load tbinstrj_b.mat; ap=adddt(ap,-169.87); % insplot(ap, 'ap');
ts = 0.01;
[imu, avp0, avp] = ap2imu(ap, ts);  imu=imulci(imu); % insplot(avp);  imuplot(imu);
avpi = lcipure(imu,avp0,0);
qvpi = [q42q3(a2qua1Batch(avpi(:,1:3))), avpi(:,4:10)]; % insplot(qvpi,'qvpi1');
% qbb1 = trjattrot(a2qua1([-0;0;0]*glv.deg),'x-30', 10*glv.dps, ts, 1);
% qbb1 = trjattrot(qbb1,'y179y-179', 10*glv.dps, ts, 10, 4); % myfig, plot(qbb1(:,end), q2attBatch(qbb1(:,1:3))/glv.deg);
% qvpi(:,1:3) = qmulBatch(qvpi(:,1:3), qbb1(1:length(qvpi),1:4));
qbb1 = trjattrot(a2qua1([0;0;0]*glv.deg),'y179y-179', 10*glv.dps, ts, 10, 4); % myfig, plot(qbb1(:,end), q2attBatch(qbb1(:,1:3))/glv.deg);
qvpi(:,1:3) = qmulBatch(qvpi(:,1:3)*0, qbb1(1:length(qvpi),1:4));
[imu1, qvp] = lciqvp2imu(qvpi,avp0(7:9),0);  % imuplot(imu1); insplot(qvpi,'qvpi1');
avpi1 = lcipure(imu1,qvp2avp(qvpi(1,:)',1),0);

lcinserrcov(qvpi(:,[1:3,7:9,end]), avp0, 0);
