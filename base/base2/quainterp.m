function q1 = quainterp(q, t)
% Attitude linear interpolation. 
% It can be denoted as att = att1+(att2-att1)*ratio
%
% Prototype: att = attinterp(att1, att2, ratio)
% Inputs: q - [q0,q1,q2,q3,t] array 
%         t - interpolation time array t or interval ts
% Output: q1 - interpolated quaternion out
%
% Example
%    q = appendt(rv2qqBatch(randn(100,3)/5),0.1); q1 = quainterp(q, 0.01);
%    myfig, plot(q1(:,5),[q1(:,1:4),normv(q1(:,1:4))],'-', q(:,5),q(:,1:4),'o');  grid on
%
% See also  avpinterp, attinterp1, qq2phi, qaddphi, q2att.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 21/05/2025
    if length(t)==1;  t=(q(1,end):t:q(end,end))';  end  % ts->t
    [rv, q0] = qq2rvBatch(q(:,1:4)); 
    rv = [[0,0,0];rv];
    rv(:,1:3) = cumsum(rv);
    rv = interp1(q(:,end), rv, t, 'spline');
    % rv = interp1(q(:,end), rv, t, 'linear');
    rv = diff(rv);
    q1 = [rv2qqBatch(rv, q0), t];
    
    