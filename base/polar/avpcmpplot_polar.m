function err = avpcmpplot_polar(avpr, avp, typ)
% Polar navigation error plot
%
% Prototype: avpcmpplot_polar(avpr, avp, typ)
% Inputs: avpr - reference AVP
%         avp - avp to compare, in the same frame with avpr among ENU/ECEF/gird
%         typ - error display frame type,
%               'n' for ENU-frame, 'e' for ECEF-frame, 'g' for grid-frame
% Output: err - AVP errors
%
% Example:
%   [imu, avp0, avp] = imupolar_grid(posset(89.9,1.0,100), 100, .1, 200, imuerrset(0,0,0,0));
%   avp0(:,4:5) = avp0(:,4:5)+1;
%   avpn = inspure(imu, avp0, 'f');       avpcmpplot_polar(avp, avpn, 'gn2g');
%   avpe = inspure_ecef(imu, avp0, 'f');  avpcmpplot_polar(avp, avpe, 'ge2g');
%   avpg = inspure_grid(imu, avp0, 'f');  avpcmpplot_polar(avp, avpg, 'gg2g');
%   avpv = inspure_nv(imu, avp0, 'f');    avpcmpplot_polar(avp, avpg, 'ge2g');
%
% See also  avptrans, avp2imu_ecef, imustatic.

% Copyright(c) 2009-2024, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 28/11/2024
global glv
    if nargin<3, typ='g'; end
    if length(typ)==4  %  usage like, err=avpcmpplot_polar(avpr, avp, 'eg2n')
        if typ(1)~=typ(4), avpr=avptrans(avpr,typ([1,3:4])); end
        if typ(2)~=typ(4), avp=avptrans(avp,typ(2:4)); end
        err = avpcmpplot_polar(avpr, avp, typ(4));
        return;
    end
    avpi = interp1(avpr(:,end), avpr(:,1:9), avp(:,end));
    idx = isnan(avpi(:,1));  avpi(idx,:)=[]; avp(idx,:)=[];
    err = avp(:,1:9)-avpi;  t=avp(:,end);  err(:,10)=t;
    idx=err(:,3)>pi/2; err(idx,3)=err(idx,3)-2*pi;
    idx=err(:,3)<-pi/2; err(idx,3)=err(idx,3)+2*pi;
    myfig;
    % trajectory(position) shown in ECEF-frame
	subplot(3,2,[1,3,5]);  % worldmap('North Pole')?
    if typ=='n', avp=avptrans(avp,'n2e'); avpi=avptrans(avpi,'n2e'); end % 'g2e' no need
    lon = atan2(avp(:,8),avp(:,7)); rho = sqrt(avp(:,8).^2+avp(:,7).^2)/glv.Re/glv.deg;
    polar(lon, rho, 'r'); view(90,90); hold on;
    polar(lon(1), rho(1), '*m');
    lon1 = atan2(avpi(:,8),avpi(:,7)); rho1 = sqrt(avpi(:,8).^2+avpi(:,7).^2)/glv.Re/glv.deg;
    polar(lon1, rho1, 'b--');  % legend('Ref', 'Calcu', 'Orientation','horizontal');
    [nxy,idx] = min(normv(avp(:,7:8)));  [nxy1,idx] = max(normv(avp(:,7:8)));
    xlabel(['90-\itL\rm / ( \circ ), ',sprintf('(%.3f ~ %.0fm)',nxy,nxy1)]);  ylabel('\it\lambda\rm / ( \circ )');
    % error shown in n/e/g-frame
    if typ=='n'
        idx=err(:,8)>pi/2; err(idx,8)=err(idx,8)-2*pi;
        idx=err(:,8)<-pi/2; err(idx,8)=err(idx,8)+2*pi;
        subplot(322), plot(t, err(:,1:3)/glv.min), xygo('datt');
        subplot(324), plot(t, err(:,4:6)), xygo('dV');
        subplot(326), plot(t, err(:,7:8)/glv.min), xygo('dll');
    elseif typ=='e' || typ=='g'
        subplot(322), plot(t, err(:,1:3)/glv.min), xygo('datt');
        subplot(324), plot(t, err(:,4:6)), xygo('dV');
        subplot(326), plot(t, err(:,7:9)), xygo('dxyz');
    end
        
