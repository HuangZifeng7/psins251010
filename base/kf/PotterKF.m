function [Xk, Dk, Pk] = PotterKF(Phikk_1, Gammak, sQk, Xk_1, Dk_1, Hk, sRk, Zk, TM)
% Potter square-root Kalman filter,
% NOTE: Rk to be vector, sQk maybe vector with Gammak=[].
%       TM='T' for time update; ='M' for measure update; ='B' time&meas update
%       If input sQk,sRk to be single datatype, then process single-point Potter SR-KF
%
% See also CommonKF.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 11/04/2025
    if nargin<8, TM='B'; end
    if TM=='M' % measure update only
        Xkk_1 = Xk_1; Dkk_1 = Dk_1;
    else  % T & B
        Xkk_1 = Phikk_1*Xk_1;
        if isempty(Gammak), GsQ=[diag(sQk);zeros(length(Phikk_1)-length(sQk),length(sQk))];
        else, GsQ=Gammak*sQk; end
        %[~, Dkk_1] = myqr([Phikk_1*Dk_1, GsQ]');  Dkk_1=Dkk_1';  % Dkk_1 is lower tri-matrix
        %[~, Dkk_1] = qr_householder([Phikk_1*Dk_1, GsQ]');  Dkk_1=Dkk_1';
        [~, Dkk_1] = qr([Phikk_1*Dk_1, GsQ]');  Dkk_1=Dkk_1(1:length(Phikk_1),:)';
        %Dkk_1 = myhhqr([Phikk_1*Dk_1, GsQ]');  Dkk_1=Dkk_1';
        if nargin<6, Xk=Xkk_1; Dk=Dkk_1; return; end % if no measure input
    end
    [m, n] = size(Hk);  I = eye(n,'like',sRk);
    for k=1:m  % sequential measure update
        ak = (Hk(k,:)*Dkk_1)';
        rho2 = ak'*ak+sRk(k)^2; rho = sqrt(rho2);
        gk = rho2+rho*sRk(k);
        Dak = Dkk_1*ak;
        Dk = Dkk_1-1/gk*Dak*ak';  % Dk = Dkk_1*(I-1/gk*ak*ak');   % bad
        Kk = Dak/rho2;
        Xk = Xkk_1+Kk*(Zk(k)-Hk(k,:)*Xkk_1);
        if k<m, Xkk_1=Xk; Dkk_1=Dk; end
    end
    if nargout>2, Pk = double(Dk)*double(Dk'); end

function [Q, R] = myqr(A)  % QR分解，A=Q*R, 其中Q'*Q=I，R为上三角阵
    [m, n] = size(A);
    if n>m,  error('n must not less than m.'); end
    R = zeros(n,'like',A);
    for i=1:n
        R(i,i) = sqrt(A(:,i)'*A(:,i));
        A(:,i) = A(:,i)/R(i,i);
        j = i+1:n;
        R(i,j) = A(:,i)'*A(:,j);
        A(:,j) = A(:,j)-A(:,i)*R(i,j);
    end
    Q = A;
    
function R = myhhqr(A)  % HouseHolder QR
    [m, n] = size(A);
    for k=1:n
        v = A(k:m,k);
        sigma = v(2:end)'*v(2:end);
        if sigma==0, continue; end
        mu = sqrt(A(k,k)^2+sigma);
        if A(k,k)<=0
            v(1) = A(k,k)-mu;  % v-|v|*e1
        else
            v(1) = -sigma/(A(k,k)+mu);  %!!! -sigma/(A(k,k)+mu)==A(k,k)-mu
        end
        beta = 2/(v(1)^2+sigma);
        beta = beta*v(1)^2; v = v/v(1);
        A(k:m,k:n) = A(k:m,k:n)-beta*v*(v'*A(k:m,k:n));
    end
    R = triu(A(1:n,:));

function [v, beta] = householder_vector(x)
    % 计算Householder向量v和缩放系数beta，使得 (I - beta*v*v')x = ||x||_2 * e1
    % 输入：列向量x
    % 输出：Householder向量v，系数beta
    
    n = length(x);
    sigma = x(2:n)' * x(2:n);  % 计算x(2:end)的平方和
    v = [1; x(2:n)];           % 初始化v = [1; x(2:n)]
    
    if sigma == 0
        beta = 0;  % 若x已经是基向量，无需变换
    else
        mu = sqrt(x(1)^2 + sigma);
        if x(1) <= 0
            v(1) = x(1) - mu;
        else
            v(1) = -sigma / (x(1) + mu);
        end
        beta = 2 * v(1)^2 / (sigma + v(1)^2);
        v = v / v(1);  % 归一化，使v(1)=1
    end


function [A] = householder_apply(A, v, beta)
    % 应用Householder变换 (I - beta*v*v') 到矩阵A
    % 输入：矩阵A，Householder向量v，系数beta
    % 输出：变换后的矩阵A
    
    A = A - beta * v * (v' * A);  % 高效计算，避免显式构造H


function [Q, R] = qr_householder(A)
    % 使用Householder变换计算QR分解 A = QR
    % 输入：矩阵A (m×n, m >= n)
    % 输出：正交矩阵Q，上三角矩阵R
    
    [m, n] = size(A);
    Q = eye(m);    % 初始化Q为单位矩阵
    R = A;         % 初始化R为A
    
    for k = 1:n
        x = R(k:m, k);                  % 当前列的下部分
        [v, beta] = householder_vector(x);  % 计算Householder向量
        
        % 应用到R的右下子矩阵
        R(k:m, k:n) = householder_apply(R(k:m, k:n), v, beta);
        
        % 更新Q (累积Householder变换)
        Q(:, k:m) = Q(:, k:m) - beta * (Q(:, k:m) * v) * v';
    end
    
    R = triu(R(1:n, :));  % 确保R是上三角
