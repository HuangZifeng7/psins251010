function [dd, lns] = bffe(fns, frl, intv)
% Batch bin-file frame extract, all in fixed column length.
%
% Prototype: [dd, lns] = bffe(fns, frl, intv)
% Inputs: fns = file-names from dir function
%         frl - frame length
%         intv - every 'intv' row read 1 row, intv=-1 for read all rows
% Outputs: dd - check result
%          lns - lines read for each file
%
% See also  bffr, bffc, frm2num.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 28/09/2025
    if nargin<3, intv=100*60; end  % 1min@100Hz
    if ischar(fns)
        fdir = strfind(fns,'\');
        if ~isempty(fdir), fdir=fns(1:fdir(end)); else, fdir='.\'; end
        fns = dir(fns);
    end
    dd = [];  lns = zeros(size(fns));
    timebar(1, length(fns));
    for k=1:length(fns)
        ddi = bffr([fdir,fns(k).name], frl, intv);
        dd = [dd; ddi];  lns(k) = size(ddi,1);
        timebar(1);
    end
