function imu = imuaxis(imu, ad)
% IMU axis order translation.
%
% Prototype: imu = imuaxis(imu, ad)
% Inputs: imu - the user's raw SIMU data
%         ad - axis order
% Output: imu - SIMU data after modification.
%
% Examples:
%    imu = imuaxis(imu, [-1;2;3]);
%    imu = imuaxis(imu, [-2;1;3; 4;-5;6]);
%    imu = imuaxis(imu, [5,4,6; -2;1;3]);
%
% See also  imurfu, imuidx, imurot.

% Copyright(c) 2009-2024, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 02/12/2024
    if length(ad)==3,
        if max(abs(ad))<4, ad(4:6)=sign(ad(1:3)).*(abs(ad(1:3))+3);
        else,              ad(4:6)=sign(ad(1:3)).*(abs(ad(1:3))-3);  end
    end
    imu(:,1:6) = imu(:,abs(ad)); % change order
    for k=1:6  % negative
        if ad(k)<0, imu(:,k)=-imu(:,k); end
    end
