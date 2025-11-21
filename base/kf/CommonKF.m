function [Xk, Pk] = CommonKF(Phikk_1, Gammak, Qk, Xk_1, Pk_1, Hk, Rk, Zk, TM)
% Single/Doulbe-datatype Kalman filter,
% NOTE: Rk to be vector, Qk maybe vector with Gammak=1.
%       TM='T' for time update; ='M' for measure update; ='B' time&meas update
%       If input Qk,Rk to be single datatype, then process single-point KF
%
% See also PotterKF.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 11/04/2025
    if nargin<8, TM='B'; end
    if TM=='M' % measure update only
        Xkk_1 = Xk_1; Pkk_1 = Pk_1;
    else  % T & B
        Xkk_1 = Phikk_1*Xk_1;
        if isempty(Gammak), GQG=diag([Qk;zeros(length(Phikk_1)-length(Qk),1)]);
        else, GQG=Gammak*Qk*Gammak'; end
        Pkk_1 = Phikk_1*Pk_1*Phikk_1'+GQG;
        if nargin<6, Xk=Xkk_1; Pk=Pkk_1; return; end % if no measure input
    end
    [m, n] = size(Hk);
    for k=1:m  % sequential measure update
        Kk = Pkk_1*Hk(k,:)'/(Hk(k,:)*Pkk_1*Hk(k,:)'+Rk(k));
        Xk = Xkk_1+Kk*(Zk(k)-Hk(k,:)*Xkk_1);
        Pk = Pkk_1-Kk*Hk(k,:)*Pkk_1;   % Pk = (I-Kk*Hk(k,:))*Pkk_1;
        if k<m, Xkk_1=Xk; Pkk_1=Pk; end
    end
    Pk=(Pk+Pk')/2;
