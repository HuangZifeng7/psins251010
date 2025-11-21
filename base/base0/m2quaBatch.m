function qnb = m2quaBatch(Cnb)
% Convert direction cosine matrix(DCM) to attitude quaternion.
%
% Prototype: qnb = m2quaBatch(Cnb)
% Input: Cnb - DCM from navigation-frame to body-frame
% Output: qnb - attitude quaternion
%
% Example
%   Cnb = a2matBatch(randn(10,3));  qnb = m2quaBatch(Cnb);
%
% See also  m2qua, a2mat, a2qua, m2att, q2att, q2mat, attsyn, m2rv.

% Copyright(c) 2009-2024, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 29/11/2024
    C11 = Cnb(:,1); C12 = Cnb(:,2); C13 = Cnb(:,3); 
    C21 = Cnb(:,4); C22 = Cnb(:,5); C23 = Cnb(:,6); 
    C31 = Cnb(:,7); C32 = Cnb(:,8); C33 = Cnb(:,9);
    qnb = zeros(length(C11),4);
    idx = C11>=C22+C33;  ii = idx;
    if ~isempty(idx)
        q1 = 0.5*sqrt(1+C11(idx)-C22(idx)-C33(idx));  qq4 = 4*q1;
        q0 = (C32(idx)-C23(idx))./qq4; q2 = (C12(idx)+C21(idx))./qq4; q3 = (C13(idx)+C31(idx))./qq4;
        qnb(idx,:) = [q0, q1, q2, q3];
    end
    idx = C22>=C11+C33;  ii = ii|idx;
    if ~isempty(idx)
        q2 = 0.5*sqrt(1-C11(idx)+C22(idx)-C33(idx));  qq4 = 4*q2;
        q0 = (C13(idx)-C31(idx))./qq4; q1 = (C12(idx)+C21(idx))./qq4; q3 = (C23(idx)+C32(idx))./qq4;
        qnb(idx,:) = [q0, q1, q2, q3];
    end
    idx = C33>=C11+C22;  ii = ii|idx;
    if ~isempty(idx)
        q3 = 0.5*sqrt(1-C11(idx)-C22(idx)+C33(idx));  qq4 = 4*q3;
        q0 = (C21(idx)-C12(idx))./qq4; q1 = (C13(idx)+C31(idx))./qq4; q2 = (C23(idx)+C32(idx))./qq4;
        qnb(idx,:) = [q0, q1, q2, q3];
    end
    idx = ~ii;
    if ~isempty(idx)
        q0 = 0.5*sqrt(1+C11(idx)+C22(idx)+C33(idx));  qq4 = 4*q0;
        q1 = (C32(idx)-C23(idx))./qq4; q2 = (C13(idx)-C31(idx))./qq4; q3 = (C21(idx)-C12(idx))./qq4;
        qnb(idx,:) = [q0, q1, q2, q3];
    end
    idx = qnb(:,1)<0;
    if ~isempty(idx)
        qnb(idx,:) = -qnb(idx,:);
    end
    
