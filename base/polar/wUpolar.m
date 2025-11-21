function [wU, ps] = wUpolar(r, vn, an, ts)
% wU = wUpolar(1, [-2;0;0], [0;0;0], 1);  wU/glv.deg
    s = vn*ts+an*0.5*ts^2;  ds = norm(s(1:2));
    syaw = atan2(s(1),s(2));
    ps = sqrt(r^2+ds^2 - 2*r*ds*cos(syaw));
    if ds<1e-6
        wU = s(1)/max(r,ps);  % must r>0 || ps>0
    else
        wU = asin(sin(syaw)/ps*ds);  % sin(syaw)/ps=sin(wU)/ds
        if ds>r,
            if s(1)>0, wU=pi-wU;
            elseif s(1)<0, wU=-pi-wU; end
        end
    end
    wU = wU/ts;