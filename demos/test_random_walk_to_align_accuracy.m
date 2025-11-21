% Angular&Velocity Random Walk effect to initial alignment accuracy simulation
% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 27/09/2025
glvs
[ts,pos0] = initp(1, 34);  T = 200;
ap0 = [zeros(3,1); pos0];  [wnie,g,gn,wN]=wnieg(pos0);
ARW = 0.00;  % in deg/sqrt(h)
VRW = 2;  % in ug/sqrt(Hz)
imuerr = imuerrset(0, 0, ARW, VRW);  imuerrKF = imuerrset(0, 0, 0.000, 1);
N = 20;  res = zeros(N,4);  att=zeros(3,1);
for k=1:N
    imu = imustatic(ap0, ts, T, imuerr);
    att = aligni0fitp(imu, pos0, 0);
    avp = inspure(imu,ap0,'O',0);  vN=avp(:,5);  t=avp(:,end);
    pp = polyfit(t, vN, 2);  dphiU=-2*pp(1)/(g*wN);
    res(k,:) = [att; dphiU]';
end
dy = norm([ARW*glv.dpsh/sqrt(T), 4*VRW*glv.ugpsHz/g/T^1.5]/wN)/glv.min;
myfig, plot([res(:,3)-ap0(3),res(:,4)-ap0(3)]/glv.min,'--o');  xygo('k', 'dyaw'); hline([-dy,dy]);
legend('i0fit', 'polyfit'); ptitle('std',[dy, std(res(:,3:4)/glv.min)]);
return;

% ARW to phiU
myfig; [wnie,g,gn,wN]=wnieg(glv.pos0); ARW = 0.0001*glv.dpsh; a2=[]; a1=[];
for iter=1:50
    phiE = cumsum(randn(1000,1)); vN = cumsum(g*phiE);  kk=1;
    T = 50:10:length(vN);
    for k=T
        pp = polyfit(1:k, vN(1:k), 2);
        a2(kk,:) = pp(1);  a1(kk,:) = pp(2);  kk=kk+1;
    end
    subplot(121); plot(T, ARW*a1/g/glv.sec); hold on
    subplot(122); plot(T, -ARW*2*a2/(g*wN)/glv.min); hold on
end
subplot(121); s = plotstd(3); grid on; xygo('T / s', 'phiE');
plot(T, ARW/glv.sec.*sqrt(T)*1/2*3,'--','linewidth',2);  % 3sigma
subplot(122); s = plotstd(3); grid on; xygo('T / s', 'phiU');
plot(T, ARW/wN/glv.min./sqrt(T)*3,'--','linewidth',2);  % 3sigma

% VRW to phiU
myfig; [wnie,g,gn,wN]=wnieg(glv.pos0); VRW = 1*glv.ugpsHz; a2=[]; a1=[];
for iter=1:50
    vN = cumsum(randn(500,1));  kk=1;
    T = 50:10:length(vN);
    for k=T
        pp = polyfit(1:k, vN(1:k), 2);
        a2(kk,:) = pp(1);  a1(kk,:) = pp(2);  kk=kk+1;
    end
    subplot(121); plot(T, VRW*a1/g/glv.sec); hold on
    subplot(122); plot(T, -VRW*2*a2/(g*wN)/glv.min); hold on
end
subplot(121); s = plotstd(3); grid on; xygo('T / s', 'phiE');
plot(T, VRW/g/glv.sec./(T.^0.5)*2*3,'--','linewidth',2);  % 3sigma
subplot(122); s = plotstd(3); grid on; xygo('T / s', 'phiU');
plot(T, VRW/(g*wN)/glv.min./(T.^1.5)*4*3,'--','linewidth',2);  % 3sigma

