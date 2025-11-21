function [sigma, tau] = tvar(y0, tau0, isfig)
% Calculate traditional variance series.
%
% Prototype: [sigma, tau] = tvar(y0, tau0)
% Inputs: y0 - data 
%         tau0 - sampling interval
%         isfig - figure flag
% Outputs: sigma - variance
%          tau - variance correlated time
%
% Example: 
%     y = randn(10000,1) + 0.00005*[1:10000]';
%     [sigma, tau] = tvar(y, 1);
%
% See also  avar, oavar, tavar.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 22/01/2025
    if nargin<3, isfig=1; end
    N = length(y0);
    y = y0; NL = N;
    for k = 1:log2(N)
        sigma(k,1) = std(y); % std
        tau(k,1) = 2^(k-1)*tau0;      % correlated time
        NL = floor(NL/2);
        if NL<3
            break;
        end
        y = 1/2*(y(1:2:2*NL) + y(2:2:2*NL));  % mean & half data length
    end
    if isfig==1
        myfig;
        [s, t] = avar(y0, tau0, [], 0);
        loglog(tau, [sigma,s]); xygo('\tau / s', '\sigma')
    end
