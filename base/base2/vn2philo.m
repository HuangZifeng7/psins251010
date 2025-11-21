function [phi0, en, phi] = vn2philo(vn, lti, isfig)
% Calculating misalign angles from long-time (>20min) open-loop pure SINS velocity error.
%
% Prototype: [phi0, eN, eU] = vn2phil(vn, lti, isfig)
% Inputs: vn - long-time open-loop pure SINS velocity error, in most case for static base
%         lti - latitude
%         isfig - figure flag
% Output: phi0 - misalignment between calculating navigation frame and real
%                navigation frame at the start t0
%         en - esternward/northward/upward gyro-bias estimate
%         phi - misalignment angle at the end
%
% Example:
%   ap0 = [[0;0;0]*glv.deg;glv.pos0];
%   imu = imustatic(ap0, 0.1, 3600, imuerrset([0;0.01;0.02], 0, 0.000, .0));
%   avp = inspure(imu, [q2att(qaddphi(a2qua(ap0(1:3)),[.1;.1;3]*glv.min));glv.pos0], 'O');
%   phi = vn2philo(avp(:,[4:6,end]), glv.pos0);
%
% See also  vn2phil, vn2phi, vn2phiu, vn2phistd, phiu2vn, aa2phi, vn2att.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 23/05/2025
global glv
    if nargin<3, isfig = 1; end
    if length(lti)>1, lti=lti(1); end  % lti = pos;
    if size(vn,2)>4, vn=vn(:,[4:6,end]); end
    wiet = glv.wie*vn(:,end);  sL = sin(lti); cL = cos(lti);
    AE = [ones(length(wiet),1), wiet, wiet.^2, wiet-sin(wiet), wiet.^2/2-1+cos(wiet)];
    aE = lscov(AE, vn(:,1));
    AN = [AE(:,1:2), 1-cos(wiet), AE(:,4)];
    aN = lscov(AN, vn(:,2));
    phi0 = [aN(2); -aE(2); -(aE(2)*sL+aN(3))/cL]*glv.wie/glv.g0;
    eN = (2*aE(3)-aN(2)*sL)*glv.wie^2/glv.g0;
    eU = (aN(4)+2*aE(3)*sL+aN(2)*cL^2)*glv.wie^2/glv.g0/cL;
    a = phi0(2)*sL-phi0(3)*cL;  b = -phi0(1)*sL-eN/glv.wie;
    g = phi0(1)*cL-eU/glv.wie;  x = -phi0(1)-eN*sL/glv.wie+eU*cL/glv.wie;
    phit = [ phi0(1)+a*sin(wiet)+x*(1-cos(wiet)), ...
             phi0(2)+b*wiet-a*sL*(1-cos(wiet))-x*sL*(wiet-sin(wiet)), ...
             phi0(3)+g*wiet-a*cL*(1-cos(wiet))+x*cL*(wiet-sin(wiet)) ];
    en = [0;eN;eU];
    phi = phit(end,1:3)';
    if isfig==1
        myfig, 
        subplot(211), plot(vn(:,end), [vn(:,1:2)]); xygo('VEN');
        hold on, plot(vn(:,end), [AE*aE,AN*aN], '--');
        subplot(212), plot(vn(:,end), phit/glv.min); xygo('phi');
        title(sprintf('\\phi_0=[%.3f, %.3f, %.3f]^\\prime;  \\epsilon_{N,U}=[%.4f, %.4f]\\circ/h', ...
            phi0(1)/glv.min,phi0(2)/glv.min,phi0(3)/glv.min, eN/glv.dph,eU/glv.dph));
    end

