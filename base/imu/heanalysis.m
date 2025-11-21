function mebxy = heanalysis(gyac, yaw, lat, eb, db)
% Heading effect analysis for SIMU multi angular position around vertical axis
% on static base, usually aided by single-axis turntable.
%
% Prototype: heanalysis(gyac, yaw, lat, eb, db)
% Inputs: gyac - gyro (in rad/s) & acc (in m/s^2), multi angular position
%                test on static base with known yaw angles
%         yaw - reference yaw angles in rad
%         lat - test latitude, or position [lat,lon,hgt]
%         eb,db - gyro drift & acc bias compensation
% Output: mebxy - mean of ebx,eby in dph
%          
% See also  magplot, magellipfit, imuplot, dvlplot, baroplot.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 05/09/2025
global glv
    if nargin<5, db=zeros(3,1);  end
    if nargin<4, eb=zeros(3,1);  end
    [M, N] = size(gyac);
    for k=1:3, gyac(:,k)=gyac(:,k)-eb(k); end
    if N>3, for k=1:3, gyac(:,3+k)=gyac(:,3+k)-db(k); end; end
    [wnie,g,gn] = wnieg(lat);
    wbie = zeros(M,3);  attk = zeros(M,3);
    for k=1:M
        if N<4, acc=[0;0;g]; else, acc=gyac(k,4:6)'; end
        [qnb, att, Cnb] = sv2atti([0;0;g], acc, yaw(k));
        wbie(k,:) = wnie'*Cnb;
        [qnb, att, Cnb] = dv2atti([0;0;g], wnie, acc, gyac(k,1:3)');
        attk(k,:) = att';
    end
    err = gyac(:,1:3)-wbie;
    yawerr = attk(:,3)-yaw;  yaw=yawcvt(yaw,'cc180c360');
    idx=yawerr>pi; yawerr(idx)=yawerr(idx)-2*pi;  idx=yawerr<-pi; yawerr(idx)=yawerr(idx)+2*pi;
    myfig,
    if N>3
        subplot(221), plot(yaw/glv.deg, err(:,1:2)/glv.dph, '--o');  xygo('y', 'ebxy');
        subplot(223), plot(yaw/glv.deg, yawerr/glv.min, '--o');  xygo('y', 'dyaw');
        subplot(322), plot(yaw/glv.deg, attk(:,1)/glv.deg, '--o');  xygo('y', 'p');
        subplot(324), plot(yaw/glv.deg, attk(:,2)/glv.deg, '--o');  xygo('y', 'r');
        subplot(326), plot(yaw/glv.deg, err(:,3)/glv.dph, '--o');  xygo('y', 'ebz');
    else
        subplot(221), plot(yaw/glv.deg, err(:,1:2)/glv.dph, '--o');  xygo('y', 'ebxy');
        subplot(223), plot(yaw/glv.deg, yawerr/glv.min, '--o');  xygo('y', 'dyaw');
        subplot(222), plot(yaw/glv.deg, err(:,3)/glv.dph, '--o');  xygo('y', 'ebz');
    end
    subplot(221), title(['mean=',sprintf('%.4f, ', mean(err(:,1:2))'/glv.dph),'(\circ/h)']);
    subplot(223), title(sprintf('\\psi_0=%.4f(\\circ), mean=%.3f,std=%.3f(\\prime)',yaw(1)/glv.deg,mean(yawerr)/glv.min,std(yawerr)/glv.min));
    mebxy = mean(err(:,1:2))'/glv.dph;
