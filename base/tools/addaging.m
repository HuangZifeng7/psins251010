function data = addaging(data, agMag, agTime, agDeg0, nclm)
% Add aging to n-column data.
%
% Prototype: data = addaging(data, agMag, agTime, agDeg0, nclm)
% Inputs: data - data input
%         agMag, agTime, agDeg0 - aging magnitude, time length, 
%            sin begin degree (0~90 in deg, not rad); 
%         nclm - n-columns index
%
% Example
%   addaging(appendt(randn(1000,3),0.1), [50;-100], [50;80], 0, [1;2]);
%   addaging(appendt(randn(1000,3),0.1), -10, 100, 80, [1,3]);
%
% See also  addslope, adddt.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 03/09/2025
    if nargin<5, nclm=1:size(data,2)-1; end;  M=length(nclm);
    if nargin<4, agDeg0=60; end;  if length(agDeg0)==1, agDeg0=repmat(agDeg0,M,1); end
    if nargin<3, agTime=60; end;  if length(agTime)==1, agTime=repmat(agTime,M,1); end
    if length(agMag)==1, agMag=repmat(agMag,M,1); end
    ts = mean(diff(data(:,end)));
    for m=1:M
        afa = agDeg0(m)+(0:ts/agTime(m):1)'*(90-agDeg0(m));  % make sure ts<<agTime
        y = sind(afa);  y=y*agMag(m)/(y(end)-y(1));  y=y-y(end);
        n = min(size(data,1),length(y));
        data(1:n,nclm(m)) = data(1:n,nclm(m)) + y(1:n,1);
    end
    if nargout<1
        myfig, plot(data(:,end), data(:,1:end-1)); xygo('val');
    end
