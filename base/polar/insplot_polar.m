function insplot_polar(avp, typ)
% Polar navigation result plot
%
% Prototype: insplot_polar(avp, typ)
% Inputs: avp - AVP data
%         typ - 'n' for ENU-frame, 'e' for ECEF-frame, 'g' for grid-frame
%
% Example:
%   [imu, avp0, avp] = imupolar_grid(posset(89.0,10,100), 100, 1, 2000);
%   insplot_polar(avp,'g');
%   insplot_polar(avptrans(avp,'g2n'),'n');
%   insplot_polar(avptrans(avptrans(avp,'g2n'),'n2g'),'g');
%   
% See also  avp2imu_ecef, imustatic.

% Copyright(c) 2009-2024, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 28/11/2024
global glv
    if size(avp,2)<9, avp=[avp(:,1:3),zeros(size(avp,1),3), avp(:,4:end)];  end  % [a,p,t]
    if nargin<2,
        if norm(avp(1,7:9))<glv.Re/2, typ='n'; else, typ='g'; end
    end
    if norm(avp(1,7:9))<glv.Re/2 && typ=='g', avp=avptrans(avp,'n2g'); end
    if norm(avp(1,7:9))<glv.Re/2 && typ=='e', avp=avptrans(avp,'n2e'); end
    if norm(avp(1,7:9))>glv.Re/2 && typ=='n', avp=avptrans(avp,'g2n'); end
    myfig;
    subplot(321), plot(avp(:,end), avp(:,1:2)/glv.deg), xygo('pr');
    subplot(322), plot(avp(:,end), avp(:,3)/glv.deg), xygo('y');
    subplot(323), plot(avp(:,end), [avp(:,4:6),normv(avp(:,4:6))]), xygo('V');
    if typ=='n'
        subplot(325), plot(avp(:,end), avp(:,9)), xygo('hgt');
        subplot(3,2,[4,6]), plot(avp(:,8)/glv.deg, (avp(:,7)-pi/2)/glv.deg), xygo('lon','L / +90(\circ)');
        plot(avp(1,8)/glv.deg, (avp(1,7)-pi/2)/glv.deg, '*m');
    elseif typ=='e'
        subplot(325), plot(avp(:,end), avp(:,9)-glv.Rp), xygo('z / -Rp(m)');
        subplot(3,2,[4,6]), plot(avp(:,8), avp(:,7)), xygo('y / m','x / m');
        plot(avp(1,8), avp(1,7), '*m'); set(gca,'YDir','reverse');
    elseif typ=='g'
        subplot(325), plot(avp(:,end), avp(:,9)-glv.Rp), xygo('z / -Rp(m)');
        subplot(3,2,[4,6]), 
        lon = atan2(avp(:,8),avp(:,7)); rho = sqrt(avp(:,8).^2+avp(:,7).^2)/glv.Re/glv.deg;
        polar(lon, rho, 'r'); view(90,90); hold on;
        polar(lon(1), rho(1), '*m');
        [nxy,idx] = min(normv(avp(:,7:8)));  [nxy1,idx] = max(normv(avp(:,7:8)));
        xlabel(['90-\itL\rm / ( \circ ), ',sprintf('(%.3f ~ %.0fm)',nxy,nxy1)]);  ylabel('\it\lambda\rm / ( \circ )');
    end
