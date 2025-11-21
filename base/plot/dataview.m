function dataview(varargin)
% View data in the Simulation Data Inspector.
%
% Prototype: dataview(data1,data2,...)
% Inputs: data1,data2,... - input data, make sure the last column is time tag
%          
% See also  imuview, setvals.
%
% Example
%   dataview(imu, gps);

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 03/09/2025
    Simulink.sdi.view; Simulink.sdi.clear;
    viewRun = Simulink.sdi.Run.create;
    for n=1:nargin
        N = size(varargin{1},2);
        for k=1:N-1
            tsk = timeseries(varargin{n}(:,k), varargin{n}(:,N), 'Name', sprintf('%03d',n*100+k));
            add(viewRun,'vars',tsk);
        end
    end
    % Simulink.sdi.close