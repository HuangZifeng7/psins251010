function [Cmm, theta] = amclbt(am, isfig)
% Acc & Mag installation angles calibration.
%
% Prototype:  [Cmm, theta] = amclbt(am, isfig)
% Inputs: am - =[acc,mag,t] or [gyro,acc,mag,t] data
%         isfig - result figure flag
% Outputs: Cmm - C^mreal_mmeas 
%          theta - inclination angle
%
% See also  magellipfit, magplot, imuclbt.

% Copyright(c) 2009-2024, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 06/12/2024 
global glv
    % msplot(211,[acc,normv(acc)]); msplot(212,[mag,normv(mag)]);
    if nargin<2, isfig=1; end
    [m, n] = size(am);
    if n==7||n==10, t=am(:,10);  else, t=(1:length(am))';  end
    if n>=9,
        idx=normv(am(:,1:3))<0.5*glv.dps;  % delete non-static acc by gyro
        acc=am(idx,4:6); mag=am(idx,7:9);  t=t(idx);
    else
        acc=am(:,1:3); mag=am(:,4:6);
    end
    if isfig, msplot(321, t, [acc,normv(acc)]);  msplot(2, t, [mag,normv(mag)]); end
    [acc, idx] = abnormaldel(acc, 0.001); mag(idx,:) = []; t(idx)=[];
    level = [0.1, 0.01, 0.003];
    for k=1:3  % delete non-static acc by acc norm
        na = normv(acc); mna = mean(na);
        idx = (abs(na/mna-1)<level(k));
        acc=acc(idx,1:3)/mna; mag=mag(idx,1:3); t=t(idx);
    end
%     if isfig, msplot(3, t, acc);  msplot(4, t, mag); end
    [b, a] = ar1coefs(1, 100);
    nm = sum(abs(acc),2);   nm = filtfilt(b,a,nm)-nm;  % delete by high-pass filter
    idx = abs(nm)<0.005;  acc=acc(idx,:); mag=mag(idx,:); t=t(idx);
    if isfig, msplot(3, t, [acc,normv(acc)], 't','Acc');  msplot(4, t, [mag,normv(mag)], 't','Mag'); end
    H = vv2H(acc,mag);
    X = lscov(H, ones(length(acc),1));
    Cmm = reshape(X,3,3)';  dCmm=det(Cmm);  c=sign(dCmm)*1/abs(dCmm)^(1/3);
    Cmm = Cmm*c;  Cmm=foam(Cmm);
    theta = -asin(c);
    mag = mag*Cmm';
    if isfig, msplot(5, t, [acc,normv(acc)]);  msplot(6, t, [mag,normv(mag)]);
              title(sprintf('incl=%.2f\\circ',theta/glv.deg)); end
%     if isfig, msplot(5, t, [acc,normv(acc)]);  msplot(6, t, [mag,normv(mag(:,1:2)),normv(mag(:,2:3)),normv(mag(:,[3,1]))]); end

    
function H = vv2H(v1, v2)
    H = [ v1(:,1).*v2(:,1), v1(:,1).*v2(:,2), v1(:,1).*v2(:,3), ...
          v1(:,2).*v2(:,1), v1(:,2).*v2(:,2), v1(:,2).*v2(:,3), ...
          v1(:,3).*v2(:,1), v1(:,3).*v2(:,2), v1(:,3).*v2(:,3) ];
    return;
% syms ax ay az bx by bz c11 c12 c13 c21 c22 c23 c31 c32 c33
% A=[ax; ay; az]; B=[bx; by; bz]; C = [c11 c12 c13; c21 c22 c23; c31 c32 c33];
% A.'*(C*B)
%  ax*(bx*c11 + by*c12 + bz*c13) + ay*(bx*c21 + by*c22 + bz*c23) + az*(bx*c31 + by*c32 + bz*c33)