function lv = inleverest(imu, avp0)
% SIMU inner lever estimation with single axis rotation just for 
% partial parameter.
%
% Prototype: lv = inleverest(imu, avp0)
% Inputs: imu - SIMU data
%         avp0 - initial parameters, avp0 = [att0,vn0,pos0], or see the code note
% Output: lv - SIMU inner lv parameter
%
% Example:
%   lv = inleverest(imu, [10;pos0]);
%   lv = inleverest(datacut(imu,330,350), [10;gps(1,1:3)']);
%
% See also  imulvest, imulever, imulvplot, sysclbt, inspure.

% Copyright(c) 2009-2024, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 12/12/2024
global glv
    [nn, ts, nts] = nnts(2, imu(:,end));
    if avp0(1)>pi, aT=fix(avp0(1)/ts); pos=avp0(2:4);
        avp0=[alignsb(imu(1:aT,:),pos);pos]; imu(1:aT,:)=[]; end % avp0=[alignT;pos]
    if abs(norm(avp0(1:4))-1)<1e-6, avp0(1:3)=q2att(avp0(1:4)); avp0(4)=[]; end % avp0=[qnb; ...]
    if length(avp0)<9, avp0=[avp0(1:3);zeros(3,1);avp0(4:end)]; end  % avp0=[att;pos]
    dotwf = imudot(imu, 5.0);  Cba = eye(3);  qnb = a2qua(avp0(1:3)); vn = zeros(3,1);
    [wnie, g, gn] = wnieg(avp0(7:9));
    len = length(imu);
    av = zeros(fix(len/nn),7);  kk = 1;
    timebar(nn, len, 'SIMU inner lever estimation.');
    H = zeros(3,9);   absphi = zeros(3,1);
    for k=1:nn:len-nn
        k1 = k+nn-1;
        wm = imu(k:k1,1:3); vm = imu(k:k1,4:6); t = imu(k1,end); dwb = mean(dotwf(k:k1,1:3),1)';
        [phim, dvbm] = cnscl([wm,vm]);  absphi = absphi+abs(phim);
        wb = phim/nts; fb = dvbm/nts;
        SS = lvS(Cba, wb, dwb); % fL = SS*[rx;ry;rz];  % lever arm
        fn = qmulv(qnb, fb); % - Cnb*fL;
        an = rotv(-wnie*nts/2, fn) + gn;
        vn = vn + an*nts;  % vel update
        Cnb = q2mat(qnb);
        H = H + Cnb*SS*nts;
        qnb = qupdt2(qnb, phim, wnie*nts);  % att update
        av(kk,:) = [q2att(qnb); vn; t];  kk=kk+1;
        timebar(nn);
    end
    [m, kaxis] = max(absphi);
    lv = zeros(3,1);
    if kaxis==3,     lv([1,2]) = lscov(H(:,[1,5]), vn);  % dvn = Cnb*SS * [rx;ry;rz]
    elseif kaxis==1, lv([1,3]) = lscov(H(:,[5,9]), vn);
    elseif kaxis==2, lv([2,3]) = lscov(H(:,[1,9]), vn); end
    imu = imulever(imu, lv, ts);
    avp = inspure(imu, avp0, 'f', 0);
    insplot(av, 'av');
    subplot(222), plot(avp(:,end), avp(:,4:5), '-.');  title(['LeverArm(cm): ',sprintf('%.2f, ', lv*100)]);
    subplot(224), plot(avp(:,end), avp(:,6), '-.');
    
function SS = lvS(Cba, wb, dotwb)
    U = (Cba')^-1; V1 = Cba(:,1)'; V2 = Cba(:,2)'; V3 = Cba(:,3)';
    Q11 = U(1,1)*V1; Q12 = U(1,2)*V2; Q13 = U(1,3)*V3;
    Q21 = U(2,1)*V1; Q22 = U(2,2)*V2; Q23 = U(2,3)*V3;
    Q31 = U(3,1)*V1; Q32 = U(3,2)*V2; Q33 = U(3,3)*V3;
    W = askew(dotwb)+askew(wb)^2;
    SS = [Q11*W, Q12*W, Q13*W; Q21*W, Q22*W, Q23*W; Q31*W, Q32*W, Q33*W];    