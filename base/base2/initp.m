function [ts, pos] = initp(tsfs, lat, lon, hgt)
% Sampling interval & geographic position setting.
%
% Prototype: [ts, pos] = initp(tsfs, lat, lon, hgt)
% Inputs: tsfs - sampling interval, >1 for sampling frequency
%         lat,lon,hgt - lat and lon are always in arcdeg,
%             & height is in m.
% Outputs: ts - sampling interval
%          pos - =[lat*arcdeg; lon*arcdeg; hgt]
% 
% See also  posset, avpset, llh, nnts.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 12/09/2025
global glv
    if nargin<4, hgt=0; end 
    if nargin<3, lon=0; end 
    if tsfs>1,  tsfs=1/tsfs;  end
    ts = tsfs;
    pos = [lat*glv.deg; lon*glv.deg; hgt];
