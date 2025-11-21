function [ang, w, a] = rot2ang(ang, w, ts, smn)
% Rotate to angle 'ang' with rate 'w' and smooth point 'smn'.
%
% Prototype: [ang, w, a] = rot2ang(ang, w, ts, smn)
% Inputs: ang - angle to rotation
%         w - angular rate
%         ts - sampling interval
%         smn - smooth points
% Outputs: ang - [angle,t]
%          w/a - angular rate/jert
%
% Example
%   [ang, w, a] = rot2ang(90*glv.deg, 20*glv.dps, 0.01, 100);
%   myfig, plot(ang(:,2), [ang(:,1)/glv.deg, w/glv.dps, a/glv.dpss]); grid on
%
% See also  attrottt, rxyz, trjattrot.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 13/05/2025
    if nargin<4, smn=10; end
    T = abs(ang/w);  Ti = fix(T/ts);  w = ang/(Ti*ts);
    w = [zeros(smn,1); repmat(w,Ti,1); zeros(smn,1)];
    w = smooth(w,smn);
    t = (1:length(w))'*ts;
    ang = [cumsum(w)*ts, t];
    a = [0; diff(w)]/ts;
    