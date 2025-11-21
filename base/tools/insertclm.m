function data = insertclm(data0, clms, data1)
% Insert data1 to data0 by specified column.
%
% Prototype: data = insertclm(data0, clms, data1)
% Inputs: data0 - input data
%         clms - specified column
%         data1 - data to be inserted
% Output: data, output data
%
% Example
%   dd = insert(randn(10,3), 2);
%
% See also  datacut.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 12/01/2025
    if nargin<3, data1=0;  end
    clm = 1:(size(data0,2)+length(clms));
    data = zeros(size(data0,1),clm(end));
    data(:,clms) = data1;
    data(:,setdiff(clm,clms)) = data0;

