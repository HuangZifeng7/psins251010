% Schuler period length simulation.
% See also  test_SINS_static.
% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 04/06/2025
glvs
T = [];
for lat=0:89
    pos = [lat*glv.deg;0;0];
    eth = earth(pos);
    T(lat+1,:) = 2*pi*sqrt([eth.RMh,eth.RNh]/eth.g);
end
myfig
subplot(211), plot(0:89, [T,(T(:,1)+T(:,2))/2,sqrt(T(:,1).*T(:,2))]/60);
xygo('lat', 'Schuler T / min');
legend('T_{S-N}', 'T_{E-W}', 'T_{mean}', 'T_{geometric mean}');
subplot(212), plot(0:89, T(:,2)-T(:,1));
xygo('lat', 'T_{E-W} - T_{S-N} / s');
%%
[nn, ts, nts] = nnts(2, 1);
avp0 = zeros(9,1);  avp0(7)=0*glv.deg;
imu = imustatic(avp0, ts, 84.6*60);   % SIMU simulation
avp0(4:5) = 0.01;
tscalepush('t/m');
avp = inspure(imu, avp0, 'V');  % pure inertial navigation
tscalepop();
