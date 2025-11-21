function t = constts(t, ts, t0)
% Set constant time interval to time sequence t.
%
% Prototype: t = constts(t, ts, t0)
% Inputs: t - input t with un-constant time interval
%         ts - time interval
%         t0 - start time
% Output: t - output t with constant time interval
%
% See also  cnt2t.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 06/01/2025
    if nargin<3, t0=t(1,end); end
    if nargin<2, ts=mean(diff(t(:,end))); end
    t(:,end) = t0+(0:length(t)-1)'*ts;
