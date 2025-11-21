function ap = aprot(ap0, dlat, dlon, dyaw, posr, isfig)
% att & pos rotation
%
% Prototype: ap = aprot(ap0, dlat, dlon, dyaw, posr)
% Inputs: ap0 - input AP data in n-frame
%         dlat,dlon,dyaw - latitude,longitude,yaw roation angles
%         posr - reference position near ap0 to rotate
% Output: ap - output AV data
%
% Example:
%   load ap_flight.mat; ap=appendt([double(ap_a),ap_ll,double(ap_h)],ap_ts);
%   ap1 = aprot(ap, 55.98*glv.deg, 0.0*glv.deg, 90*glv.deg, [], 1);
%   insplot_polar(avptrans(insertclm(ap1,4:6),'n2g'), 'g');
%   
% See also  apmove, avprot.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 13/01/2025
global glv
    if nargin<6, isfig=0; end
    if nargin<5, posr=ap0(1,4:6)'; end
    if isempty(posr), posr=ap0(1,4:6)'; end
    if nargin<4, dyaw=0; end
    if nargin<3, dlon=0; end
    ap = ap0;
    % set rotation matrix
    Cu = pos2cen(posr);   Cu = rv2m(Cu(:,3)*dyaw);
    Clat = rv2m([cos(posr(2)-pi/2);sin(posr(2)-pi/2);0]*dlat);
    Clon = rv2m([0;0;1]*dlon);
    T = Clon*Clat*Cu;
    % rotate pos
    Cen = pos2cen(ap(:,4:6)); % old Cen
    Cen1 = m3xm3(T,Cen);  % new Cen
    ap(:,4:6) = cen2pos(Cen1,ap(:,6));
    % rotate att
    Ceb = m3xm3(Cen1,a2matBatch(ap(:,1:3)));
    Cen = pos2cen(ap(:,4:6));
    Cnb = m3xm3(Cen(:,[1,4,7,2,5,8,3,6,9]),Ceb);
    ap(:,1:3) = m2attBatch(Cnb);
    if isfig
        myfig
        subplot(221); plot(ap0(:,end),[ap0(:,1:2),ap(:,1:2)]/glv.deg); xygo('pr');
        subplot(223); plot(ap0(:,end),[ap0(:,3),ap(:,3)]/glv.deg); xygo('y');
        subplot(222); plot(ap0(:,5)/glv.deg,ap0(:,4)/glv.deg); xygo('lon','lat'); axis equal
                      plot(ap0(1,5)/glv.deg,ap0(1,4)/glv.deg,'r*');
        subplot(224);
            if isfig==1
                plot(ap(:,5)/glv.deg,ap(:,4)/glv.deg); xygo('lon','lat'); %axis equal
                plot(ap(1,5)/glv.deg,ap(1,4)/glv.deg,'r*');
            elseif isfig==2  % polar plot
                subplot(224);
                avp = avptrans(insertclm(ap,4:6),'n2e');
                lon = atan2(avp(:,8),avp(:,7)); rho = sqrt(avp(:,8).^2+avp(:,7).^2)/glv.Re/glv.deg;
                h=polar(lon, rho); set(h,'linewidth',1.5); view(90,90); hold on;
                polar(lon(1), rho(1), '*m');
                xygo('lat','lon');
            end
    end
