function [dlon, rk, sy] = dllpolar(rk_1, vn, an, ts)
% dlon - delta_longitide
% ps - ps/RN=delta_latitude
% [dlon, ps] = dllpolar(10, [2;0;0], [0;0;0], 1);  dlon/glv.deg
%     xy = vn*ts+an*0.5*ts^2;  ds = norm(xy(1:2));  % horizontal moving distance
%     syaw = atan2(xy(1),xy(2));
%     ps = sqrt(ps0^2+ds^2 - 2*ps0*ds*cos(syaw)); % pole to s distance
%     if max(ps0,ps)<1e-6
%         dlon = 0;
%     else
%         dlon = asin(sin(syaw)/ps*ds);  % sin(syaw)/ps=sin(dlon)/ds
%         if ds^2>ps0^2+ps^2,  % obtuse angle
%             if xy(1)>0, dlon=pi-dlon;
%             elseif xy(1)<0, dlon=-pi-dlon; end
%         end
%     end
    xy = vn*ts+an*0.5*ts^2;  ds = norm(xy(1:2));  % horizontal moving distance
    syaw = atan2(xy(1),xy(2));
    rk = sqrt(rk_1^2+ds^2 - 2*rk_1*ds*cos(syaw));
    r0 = 10.0e-0;
%     if rk_1<=r0 && rk<=r0  % both in
%         dlon = 0;  rk=rk+ds;  % distance cumulate
%     elseif rk_1>r0 && rk<=r0  % go in
%         dlon = pi;
%     elseif rk_1<=r0 && rk>r0  % go out
%         dlon = 0;
%     else  % both out
%         dlon = sign(xy(1))*acos((rk_1^2+rk^2-ds^2)/(2*rk_1*rk));
%     end
    rr = rk_1*rk;  rr0 = r0^2;
    if rr<rr0, rr=rr+rr0;  end
    val = (rk_1^2+rk^2-ds^2)/(2*rr);
    if val>1, val=1; elseif val<-1, val=-1; end
    dlon = sign(xy(1))*acos(val);
    sy = xy(2);
    glvdata([rk_1;ds;syaw;xy]);
