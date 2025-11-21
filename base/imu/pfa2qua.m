function qrb = pfa2qua(pfa, seqstr)
% Platform frame angles to attitude quaternion.
%
% Prototype: qrb = pfa2qua(pfa, seqstr)
% Inputs: pfa - platform frame angle array
%         seqstr - rotation sequence string
% Output: qrb - attitude quaternion output
%
% Example:
%    qrb = pfa2qua(randn(10,3), 'zxy'); 
%
% See also  attrottt, trjattrot, rot2ang, a2qua, a2qua1.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 18/05/2025
    qrb = [ones(length(pfa),1), zeros(length(pfa),3)];
    for k=1:length(seqstr)
        [C, qk] =  rxyz(pfa(:,k),seqstr(k));
        qrb = qmulBatch(qrb, qk);
    end
