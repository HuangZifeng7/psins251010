function lcinserrcov(qp, ap0, A0, avperr, imuerr)
% Launch-centered inertial INS error covariances analysis.
% Ref. 'Yan G. 一种基于实测飞行轨迹的惯导误差分析方法研究, NPT,2024'.
%
% Example
%
% See also  inserrcov, insupdate, kfupdate, sinsgps, kfstat.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 12/05/2025
global glv;
    if nargin<3, A0=0; end
    if size(ap0,1)>7; ap0=ap0([1:3,7:9]); end
    ts = diff(qp(1:2,end));
    qvp = lciqp2qvp(qp, ts);
    [imu, qvp] = lciqvp2imu(qvp, ap0(4:6), A0); % imuplot(imu), insplot(qvp,'qvpi');
    %%
    len = length(imu);
    N = 33; L = 6;
    [nn,ts,nts,nts2] = nnts(1,diff(imu(1:2,end)));
    Cn0a = a2mat([0;A0;0]); % [F;U;R]
    [re0,Ce0n0] = blh2xyz(ap0(4:6));  Ce0n0 = CenNUE(Ce0n0);  Ce0a = Ce0n0*Cn0a;  ra0 = Ce0a'*re0;
    qab = a2qua1(ap0(1:3));  % ^a for launch inertial frame, _b for body frame
    va = Ce0a'*cross([0;0;glv.wie],re0);  pa = zeros(3,1);
    %%
    Ft = zeros(N);  Phik0 = eye(N);  Gamma = zeros(N,L);
    Pjk = []; for j=1:N, Pjk{j} = zeros(N); end
    Qjk = []; for j=1:L, Qjk{j} = zeros(N); end
    PQk = zeros(N,N+L);  Pk = zeros(N,N+L,len);
    ki = timebar(1, len, 'Trajectory-based INS error evaluation.');
    for k=1:len
        k1 = k+nn-1;  wvm = imu(k:k1,1:6);  t = imu(k1,end);
        [phim, dvbm] = cnscl(wvm);  wib = phim/nts; fb = dvbm/nts;
        Ceka = rv2m([0;0;-glv.wie*(t-nts2)])*Ce0a;  % C^ek_e0*C^e0_a
        re = Ceka*(ra0+pa+va*nts2);
        fe = gravj4(re);  % gravity calculate
        va1 = va + qmulv(qab,dvbm) + Ceka'*fe*nts;  % vel update
        pa = pa + (va+va1)*nts2;  va=va1;  % pos update
        qab = qupdt(qab, phim);  % att update
        %%
        Cab = q2mat(qab);  nre = norm(re);  beta2 = norm(fe)/nre;  ue = re/nre;
        Ft(1:3,10:12) = -Cab;
        Ft(4:6,[1:3,7:9,13:15]) = [askew(Cab*fb),-Ceka'*beta2/2*(glv.I33-3*ue*ue')*Ceka,Cab];
        Ft(7:9,4:6) = eye(3);
        Ft(1:3,16:24) = [-wib(1)*Cab, -wib(2)*Cab, -wib(3)*Cab];
        Ft(4:6,25:33) = [fb(1)*Cab, fb(2)*Cab(:,2:3), fb(3)*Cab(:,3), Cab*diag(fb.^2)];    
        Phikk_1 = eye(N)+Ft*ts;  Gamma(1:6,1:6) = eye(6);
        Phik0 = Phikk_1*Phik0;
        for j=1:N,  Pjk{j} = Phik0(:,j)*Phik0(:,j)';  PQk(:,j) = diag(Pjk{j});  end
        for j=1:L,  Qjk{j} = Phikk_1*Qjk{j}*Phikk_1' + Gamma(:,j)*Gamma(:,j)';  PQk(:,N+j) = diag(Qjk{j});  end
        Pk(:,:,k) = PQk;  % P,Q 引起状态误差
        ki = timebar;
    end
    %%
    if nargin<5,  imuerr=imuerrset(0.01,50,0.001,5,  0,0,0,0,  30,20, 5,5, 5);  end
    if nargin<4,  avperr=avperrset([0.3;3;0.3], 0.0, 0.0);  end
    err = [avperr; imuerr.eb; imuerr.db; imuerr.dKga; imuerr.Ka2; [imuerr.web; imuerr.wdb]*sqrt(ts)].^2;
    Pkt = zeros(len,N);  Pki = zeros(len,39,9);
    for k=1:len
        s = zeros(N,1);
        for j=1:N+L
            s = s+Pk(:,j,k)*err(j);
        end
        Pkt(k,:) = s';   % k 时刻各状态误差
        for j=1:9
            Pki(k,:,j) = Pk(j,:,k).*err';  % 姿态、速度、位置误差
        end
    end
    sPkt = sqrt(Pkt);  t = imu(:,end);
    myfig
    subplot(311), plot(t, sPkt(:,1:3)/glv.min); xygo('phi');  legend('\phi_x', '\phi_y', '\phi_z', 'Location','Best');  title('( a )');
    subplot(312), plot(t, sPkt(:,4:6)); xygo('dV');  legend('\deltaV_x', '\deltaV_y', '\deltaV_z', 'Location','Best');  title('( b )');
    subplot(313), plot(t, [sPkt(:,7:8)*glv.Re/glv.Re,sPkt(:,9)]); xygo('dP');  legend('\deltax', '\deltay', '\deltaz', 'Location','Best');  title('( c )');
    %
    err = [];
    for j=1:9,  err(j,:) = Pki(end,:,j);  end;  err(:,end+1) = sum(err,2);  err=sqrt(err);  % 最后时刻误差，总误差
    err(1:3,:)=err(1:3,:)/glv.min; err(7:8,:)=err(7:8,:)*glv.Re/glv.Re; 
    myfig;
    % subplot(311), plot(err(1:3,:)','-.o','linewidth',2); xlim([1,40]); grid on; ylabel('\it\phi\rm / ( \prime )'); legend('\phi_E', '\phi_N', '\phi_U', 'Location','Best');
    subplot(311), plot(err(1:3,:)','-.o','linewidth',2); xlim([1,40]); grid on; ylabel('\it\phi\rm / ( \prime )'); legend('\phi_x', '\phi_y', '\phi_z', 'Location','Best');
    % xtl = {'phiE/N/U', 'dVE/N/U', 'dLat/Lon/Hgt', 'ebx/y/z', 'dbx/y/z', 'dkgx/y/zx', 'dkgx/y/zy', 'dkgx/y/zz', 'dkax/y/zx', 'dkay/zy,zz', 'KA2x/y/z', 'wgx/y/z', 'wax/y/z', 'Total'};
    xtl = {'phix/y/z', 'dVx/y/z', 'dx/y/z', 'ebx/y/z', 'dbx/y/z', 'dkgx/y/zx', 'dkgx/y/zy', 'dkgx/y/zz', 'dkax/y/zx', 'dkay/zy,zz', 'KA2x/y/z', 'wgx/y/z', 'wax/y/z', 'Total'};
    set(gca, 'xtick', [2:3:N+L,40], 'XTicklabel', xtl); title('( a )');
    % subplot(312), plot(err(4:6,:)','-.o','linewidth',2); xlim([1,40]); grid on; ylabel('\delta\itV\rm / ( m/s )'); legend('\deltaV_E', '\deltaV_N', '\deltaV_U', 'Location','Best');
    subplot(312), plot(err(4:6,:)','-.o','linewidth',2); xlim([1,40]); grid on; ylabel('\delta\itV\rm / ( m/s )'); legend('\deltaV_x', '\deltaV_y', '\deltaV_z', 'Location','Best');
    set(gca, 'xtick', [2:3:N+L,40], 'XTicklabel', xtl); title('( b )');
    % subplot(313), plot(err(7:9,:)','-.o','linewidth',2); xlim([1,40]); grid on; ylabel('\delta\itP\rm / m'); legend('\deltaL', '\delta\lambda', '\deltaH', 'Location','Best');
    subplot(313), plot(err(7:9,:)','-.o','linewidth',2); xlim([1,40]); grid on; ylabel('\delta\itP\rm / m'); legend('\deltax', '\deltay', '\deltaz', 'Location','Best');
    set(gca, 'xtick', [2:3:N+L,40], 'XTicklabel', xtl); title('( c )');
    % [xkpk, kfs, trj] = tbinseval(ap1, avperr, imuerr, 1); 
    %
    nn = ceil(1/diff(t(1:2)));
    t1 = t(1:nn:end);
    % zlb={'phiE','phiN','phiU','dVE','dVN','dVU','dLat','dLon','dH'};
    zlb={'phix','phiy','phiz','dVx','dVy','dVz','dx','dy','dz'};
    % zlb1={'phiE / \prime','phiN / \prime','phiU / \prime','dVE / m/s','dVN / m/s','dVU / m/s','\delta\itL\rm / m','\delta\it\lambda\rm / m','\delta\itH\rm / m'};
    zlb1={'phix / \prime','phiy / \prime','phiz / \prime','dVx / m/s','dVy / m/s','dVz / m/s','\delta\itx\rm / m','\delta\ity\rm / m','\delta\itz\rm / m'};
    hNull = figure; set(hNull,'Visible','off');
    for k=1:9
        sPki = sqrt(Pki(:,:,k));
        if k<4, sPki = sPki/glv.min;  % 1,2,3 /glv.min
        elseif k==7 || k==8, sPki = sPki*glv.Re/glv.Re;  end % 7,8 *glv.Re
        z = interp2(1:39, t1, sPki(1:nn:end,:), 1:0.1:39, t1);
        h=figure(k+10); set(h, 'WindowStyle','docked','NumberTitle','off','Name',zlb{k});
        mesh(1:0.1:39, t1, z); view(59,53);
        % xtl = {'phiE/N/U', 'dVE/N/U', 'dLat/Lon/Hgt', 'ebx/y/z', 'dbx/y/z', 'dkgx/y/zx', 'dkgx/y/zy', 'dkgx/y/zz', 'dkax/y/zx', 'dkay/zy,zz', 'KA2x/y/z', 'wgx/y/z', 'wax/y/z'};
        xtl = {'phix/y/z', 'dVx/y/z', 'dx/y/z', 'ebx/y/z', 'dbx/y/z', 'dkgx/y/zx', 'dkgx/y/zy', 'dkgx/y/zz', 'dkax/y/zx', 'dkay/zy,zz', 'KA2x/y/z', 'wgx/y/z', 'wax/y/z'};
        set(gca, 'xtick', [2:3:N+L], 'XTicklabel', xtl);
        ylabel('t / s'); zlabel(zlb1{k});
    end
    close(hNull);
