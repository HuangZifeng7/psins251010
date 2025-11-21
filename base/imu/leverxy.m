function [lx, ly] = leverxy(imu, avp0, isfig)
% Inner lever arm estimation for level x/y acc, rotated by vertical fixed axis.
%
% Prototype: [lx, ly] = leverxy(imu, avp0, isfig)
% Inputs: imu - SIMU data array
%         avp0 - initial AVP
%         isfig - figure flag
% Outputs: lx,ly - [lxx;lxy;0], [lyx;lyy;0] with x/y acc to fixed axis
%
% Example:
%    [lx, ly] = leverxy(imu, avp0, 1);
%
% See also  pp2lever, inslever, imulever, imuclbt.

% Copyright(c) 2009-2021, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 21/01/2021
    if nargin<3, isfig=1; end
    ts = diff(imu(1:2,end));  len=length(imu);  % imuplot(imu), insplot(avp);
    CSS = zeros(3,9);
    Hx = zeros(len,9);  Hy = Hx;
    [wnie, g, gn] = wnieg(avp0(end-2:end));  wniets = wnie*ts;  Cnn = rv2m(-wniets/2);
    qnb = a2qua(avp0(1:3)); vn = zeros(3,1); vnk = zeros(len,3);
    for k=2:len
        [phim, dvbm] = cnscl(imu(k,1:6));
        Cnb = q2mat(qnb);
        dvn = Cnn*Cnb*dvbm;
        vn = vn + dvn + gn*ts;
        qnb = qupdt2(qnb, phim, wniets);
        vnk(k,:) = vn';
        wb = imu(k,1:3)'/ts;
        dwb = (imu(k,1:3)-imu(k-1,1:3))'/ts/ts;
        SS = imulvS(wb, dwb); % fL = SS*[clbt.rx;clbt.ry;clbt.rz];
        CSS = CSS+Cnb*SS*ts;
        Hx(k,:) = CSS(1,:);  Hy(k,:) = CSS(2,:);
    end
    lxy = lscov([Hx(:,[1:2,4:5]);Hy(:,[1:2,4:5])], [vnk(:,1);vnk(:,2)]);
    lx = [lxy(1:2);0];  ly = [lxy(3:4);0];
    if isfig==1
        myfig, plot(imu(:,end), [Hx(:,[1:2,4:5])*lxy,Hy(:,[1:2,4:5])*lxy, vnk(:,1:2)]); xygo('dv');
        title([sprintf('%.3f, ',lxy*100),'(cm)']);
    end
    
