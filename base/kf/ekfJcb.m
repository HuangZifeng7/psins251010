function [Jcb, y] = ekfJcb(hfx, x, tpara)
% Jacobian matrix & state transformation for nonlinear equation.
%
% Prototype: [Jcb, y] = ekfJcb(hfx, x, tpara)
% Inputs: hfx - a handle for nonlinear state equation
%         x - state vector
%         tpara - some time-variant parameter pass to hfx
% Outputs: Jcb - Jacobian matrix of hfx
%          y - y = fhx(x)
%
% See also  ekf, ckf, ukf, Jacob5, ukfUT, alignvn_ekf.

% Copyright(c) 2009-2022, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 04/03/2022
    [Jcb, y] = feval(hfx, x, tpara);
