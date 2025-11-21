%https://www.psins.org.cn/newsinfo/2370168.html
glvs

dd = load('data2.txt');

ts=0.01; tt = dd(:,1);

imu = [[dd(:,2:4)*glv.dps,dd(:,5:7)]*ts,tt]; imuplot(imu)

avp0 = [dd(:,8:10)*glv.deg,dd(:,11:13),dd(:,14:15)*glv.deg,dd(:,16),tt];

avp0(:,3) = yawcvt(avp0(:,3),'c360cc180');

insplot(avp0(65/ts:end,:));

gps = [dd(:,17:19),dd(:,20:21)*glv.deg,dd(:,22),tt]; gps = gps(gps(:,4)>0,:); gpsplot(gps);

gpsYaw = [dd(:,23), tt]; gpsYaw = gpsYaw(gpsYaw(:,1)>0,:);

figure, plot(avp0(:,end), avp0(:,3)/glv.deg, gpsYaw(:,2), gpsYaw(:,1)); legend('Ref Yaw', 'GPS Yaw'); xygo('GPSYaw / \circ');

magplot([dd(:,26:28),tt]);

baroplot([dd(:,29),tt], gps(:,[end-1,end]));