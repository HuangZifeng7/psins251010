function [res, fns] = bffc(fns, frl, cntIdx, ovf)
% Batch bin-file frame count check, all in fixed column length.
%
% Prototype: [res, fns] = bffc(fns, frl, cntIdx, ovf)
% Inputs: fns = file-names from dir function
%         frl - frame length
%         cntIdx - count index
%         ovf - overflow number of cntIdx
% Outputs: res - check result
%          fns - file-names and result
%
% See also  bffr, bffe.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 28/09/2025
    if nargin<4, ovf=-255; end
    if nargin<3, cntIdx=frl; end
    if isempty(cntIdx), cntIdx=frl; end 
    ovf1 = 1-ovf;
    if ischar(fns)
        fdir = strfind(fns,'\');
        if ~isempty(fdir), fdir=fns(1:fdir(end)); else, fdir='.\'; end
        fns=dir(fns);
    end
    res = zeros(length(fns),4);
    timebar(1, length(fns));
    for k=1:length(fns)
        dd = bffr([fdir,fns(k).name], frl, 2);
        frm = fns(k).bytes/frl;
        res(k,:) = [dd(1,cntIdx),dd(2,cntIdx), frm, mod(frm+dd(1,cntIdx)-1,ovf1)];
        fns(k).frmlen = frm;  fns(k).lost = res(k,2)-res(k,4);
        timebar(1);
    end
    myfig, plot([fns.lost]); xygo('k', 'diff');
