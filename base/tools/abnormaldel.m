function [data, idx] = abnormaldel(data, ratio)
% Delete abnormal value(outlier) from the data.
%
% Prototype: data = abnormaldel(data, ratio)
% Inputs: data - data input
%         ratio - abnormal ratio = 0~0.5, if >1 for count
% Output: data - data output without abnormal
%         idx - idnex for abnormal data
%
% Example:
%   d0=appendt(randn(100,1),1);  [d1,idx]=abnormaldel(d0(:,1)); d1=d0; d1(idx,:)=[];
%   myfig, plot(d0(:,2),d0(:,1), '-', d1(:,2), d1(:,1), 'o')
%
% See also abnormaladd, smoothol.

% Copyright(c) 2009-2024, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 07/12/2024
    if nargin<2, ratio=0.1; end
    if size(data,2)>1,  nm = normv(data);
    else,               nm = data;    end
    [~, idx] = sort(nm);
    if ratio<0.5,  n = fix(length(idx)*ratio);
    else,          n = fix(ratio);  end
    if n<1, n=1; end
    idx = [idx(1:n);idx(end-n+1:end)];
    data(idx,:)=[];
    