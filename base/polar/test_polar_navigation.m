% function test_polar_navigation()
% test_polar_navigation
glvs
%% static near pole
ts=0.01; avp0 = avpset([0;0;0], 0, [89.995+0.005; 0; 0]);
imu = imustatic(avp0, ts, 500, imuerrset(0,0,0.1,0));
avp = appendt(repmat(avp0',length(imu),1),ts);
avp0 = avpadderr(avp0, avperrset([1;1;10],[1;1;0],[1000;500;0]));
glvdata(10,6);  avpn = inspure_enu(imu, avp0, 'V');  avpcmpplot_polar(avp, avpn, 'nn2e');  glvdata; plotn(glv.gdata);
avpg = inspure_grid(imu, avp0, 'V'); avpcmpplot_polar(avpn, avpg, 'ng2g');
%% enter->leave pole
ts=0.01; avp0 = avpset([0;0;0], 0, [89.999; 0; 0]);
imu = imustatic(avp0, ts, 50, imuerrset([0;0;3600],0,0,0));
T=17; t=(0:ts:T)'; dbx=0.25*glv.g0*sin(2*pi/T*t); lt=length(t);
imu(1:lt,5)=imu(1:lt,5)+dbx*ts;  imu(lt+1:2*lt,4)=imu(lt+1:2*lt,4)+dbx*ts;
avpg = inspure_grid(imu, avp0, 'V');
avpn = inspure_enu(imu, avp0, 'V');  avpcmpplot_polar(avpg, avpn, 'gn2e');
%% constant course line
ts = 0.01;  avp0 = avpset(-45, 10, 89.99);
seg = trjsegment([],  'init',     0);
seg = trjsegment(seg, 'uniform', 300);
trj = trjsimu(avp0, seg.wat, ts, 1);  % imuplot(trj.imu); insplot(trj.avp);
avpe = avptrans(trj.avp(1:800:end,:));  % insplot_polar(avpe,'e'); exclude singular point
avpe = ap2avp_ecef(avpe, ts); % insplot_polar(avpe,'e');
[imu, avp0] = avp2imu_ecef(avpe); % imuplot(imu);
avpn = inspure_enu(imu, avp0);  avpcmpplot_polar(avpe, avpn, 'en2e');
%% cross pole in straight line (grid trajectory)
[imu, avp0, avp] = imupolar_grid(xyz2blh([1500;-.1;glv.Rp+100]), 100, .01, 30, imuerrset(0,0,0,0)); % imuplot(imu);
avp0 = avpadderr(avp0, avperrset([1;1;10],[1;1;0],[1000;500;0])); % insplot_polar(avp,'g');
avpn = inspure_enu(imu, avp0, 'V');  avpcmpplot_polar(avp, avpn, 'gn2e');
avpg = inspure_grid(imu, avp0, 'V'); avpcmpplot_polar(avp, avpg, 'gg2g');
avpe = inspure_ecef(imu, avp0, 'f'); avpcmpplot_polar(avp, avpe, 'ge2g');
avpv = inspure_nv(imu, avp0, 'f');   avpcmpplot_polar(avp, avpv, 'ge2e'); % avpcmpplot_polar(avpe, avpv, 'ee2e');
avpt = inspure_tenu(imu, avp0, 'f'); avpcmpplot_polar(avp, avpt, 'ge2e'); % avpcmpplot_polar(avpe, avpv, 'ee2e');
%% turn round pole (enu-trajectory)
[imu, avp0, avp] = imupolar_enu(posset(89.995,0,0), 10, .01, 600); % insplot_polar(avp,'n');
avp0 = avpadderr(avp0, avperrset([1;1;10],[1;1;0],[1000;500;0]));
avpn = inspure_enu(imu, avp0, 'V');  avpcmpplot_polar(avp, avpn, 'nn2e');
avpg = inspure_grid(imu, avp0, 'V'); avpcmpplot_polar(avp, avpg, 'ng2g');
%% random flight near pole
[imu, avp0, avp] = imupolar_rnd(0.01, 1000.0, 100.0, 200, (89.9999+0.0001)*glv.deg);  % imuplot(imu)
imu = imuadderr(imu, imuerrset(0,0,10,10));
avp00 = avpadderr(avp0, avperrset([1;1;10],[1;1;0],[10;500;0]));
avpn = inspure_enu(imu, avp0);   avpcmpplot_polar(avp, avpn, 'gn2e');
avpt = inspure_tenu(imu, avp0);   avpcmpplot_polar(avp, avpt, 'gt2e');
avpg = inspure_grid(imu, avp0, 'f');  avpcmpplot_polar(avp, avpg, 'gg2e');
avpe = inspure_ecef(imu, avp0, 'f');  avpcmpplot_polar(avp, avpe, 'ge2e');
avpv = inspure_nv(imu, avp0, 'f');   avpcmpplot_polar(avp, avpv, 'ge2e'); % avpcmpplot_polar(avpe, avpv, 'ee2e');
avpv1 = inspure_nv1(imu, avp0, 'f');   avpcmpplot_polar(avp, avpv1, 'ge2e'); % avpcmpplot_polar(avpe, avpv, 'ee2e');
%% trajectory move from low latitude
load ap_flight.mat; ap=appendt([double(ap_a),ap_ll,double(ap_h)],ap_ts);  % insplot(ap,'ap');
[imu, avp0, avp] = imupolar_move(datacut(ap,0,600), 0.01, 89.9818*glv.deg-0/glv.Re);  % imuplot(imu);
avpg = inspure_grid(imu, avp0);  avpcmpplot_polar(avp, avpg, 'gg2e');
avp1 = avp0; avp1(1:3)=avp1(1:3)+[1;1;30]*glv.min;
avpg1 = inspure_grid(imu, avp1); 
avpn1 = inspure_enu(imu, avp1);   avpcmpplot_polar(avpg1, avpn1, 'gn2e');

