function q = rv2qqBatch(rv, q0)
% Calculate quaternion rotation vector by patch processing.
%
% Prototype: q = rv2qqBatch(rv, q0)
% Inputs: rv - rotation vector serial between [qnbk_1; qnbk]
%         q0 - intital attitude quaternion
% Output: q - quaternion serial
%
% See also  qq2rvBatch, qq2rv, qq2phi, qq2afa.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 21/05/2025
    if nargin<2, q0=[1;0;0;0]; end
    q = zeros(length(rv)+1,4);  q(1,:) = q0';
    for k=2:length(q)
        q(k,:) = qupdt(q(k-1,:)', rv(k-1,:)');
    end
    