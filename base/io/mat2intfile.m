function data = mat2intfile(intfile, matfile)
% Convert '.mat' to int-type file '.izp' (IMU zip).
%
% Prototype: data = mat2intfile(intfile, matfile)
% Inputs: intfile - int file name, with extension '.izp'
%         matfile - mat file name (will be) in mat file, or array
% Output: data - data array read from the '.mat' file

% See also  matbinfile.

% Copyright(c) 2009-2024, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 11/12/2024
    if nargin<2, matfile=[intfile(1:end-4),'.mat']; end
    if ischar(matfile)
        S = load(matfile, '-mat');
        names = fieldnames(S);
        if length(names)>1, error('too many variables in .mat file!'); end
        data = getfield(S, names{1});
    else
        data = matfile;
    end
    fid = fopen(intfile,'w');
    fwrite(fid, data', class(data));
    fclose(fid);

