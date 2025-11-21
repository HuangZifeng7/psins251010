function avp = inspure_nv1(imu, avp0, href, isfig)
% Process Normal Vector based SINS pure inertial navigation with SIMU log data and
% using initial condition avp0 = [att0,vn0,pos0].
%
% Prototype: avp = inspure_ecef(imu, avp0, href, isfig)
% Inputs: imu - SIMU data array
%         avp0 - initial parameters, avp0 = [att0,vn0,pos0]
%         href -    'v' - velocity fix-damping, =vn0
%                   'V' - vertical velocity fix-damping, =vn0(3)
%                   'p' - position fix-damping, =pos0
%                   'P' - position fix-damping, =pos0 & vertical velocity fix-damping, =vn0(3)
%                   'H' - height fix-damping, =pos0(3)
%                   'f' - height free.
%         isfig - figure on/off flag
% Output: avp - navigation results, avp = [attE,ve,pe,t]
%
% Example:
%   [imu, avp0, avp] = imupolar_ecef(posset(89.9,10,100), 100, .1, 200);
%   avpe = inspure_nv1(imu, avp0, 'V');
%   avpcmpplot_polar(avp, avpe, 'e');
%
% See also  insinit_ecef, insupdate_ecef, inspure, inspure_ecef, inspure_grid.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 05/02/2024
global glv
    if nargin<3, href='f'; end
    vp_fix = href(1);
    [nn, ts, nts] = nnts(2, imu(:,end));
    %ins0 = insinit_nv(avp0, ts);  
    ins = insinit_nv1(avp0, ts);  vn0 = avp0(4:6); pos0 = avp0(7:9);  ve0 = ins.ve; ne0 = ins.ne;
    len = length(imu);    avp = zeros(fix(len/nn), 10);
    ki = timebar(nn, len, 'Normal-vector pure inertial navigation processing.');
    for k=1:nn:len-nn+1
        k1 = k+nn-1;
        wvm = imu(k:k1, 1:6);  t = imu(k1,end);
        %ins0 = insupdate_nv(ins0, wvm);
        ins = insupdate_nv1(ins, wvm);
        if vp_fix=='v',      ins.ve = ve0;
        elseif vp_fix=='V',  ins.ve = ins.ve - (ins.ve-[0;0;ve0(3)])'*ins.ne*ins.ne;
        elseif vp_fix=='p',  ins.ne = ne0;  ins.h = pos0(3);
        elseif vp_fix=='P',  ins.ne = ne0;  ins.h = pos0(3);  ins.ve = ins.ve - ins.ve'*ins.ne*ins.ne;
        elseif vp_fix=='H',  ins.h = pos0(3);
        elseif vp_fix=='f',  ins.ve(3) = ins.ve(3);  % free, no need
        end
        avp(ki,:) = [ins.avp; t]';   %  avp(ki,1:3) = sum(imu(k:k1, 4:6))+ins.eth.gn'*0.01;
        ki = timebar;
    end
    if nargin<4, isfig=1; end
    if isfig==1,
        insplot_polar(avp,'e');
    end
