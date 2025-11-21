function lost = fctest(cnt)
% Frame count continuity test.
%
% Prototype: fctest(cnt)
% Inputs: t - time stamp
%         xist - x-axis is t, default 0
%
% See also  ttest, cnt2t.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 28/09/2025
    t = cnt2t(cnt(:,end));  dt = diff(t);
    lost = sum(dt)-length(dt);
    myfig,
    subplot(311), plot(cnt(:,end));  xygo('k', 'k / frame count');
    subplot(312), plot(t);  xygo('k', 'k');
    subplot(313), plot(dt);  xygo('k', 'count diff');
    if lost>0
        idx = find(dt>1);
        subplot(311), hold on, plot(idx,cnt(idx,end),'o');
        subplot(312), hold on, plot(idx,t(idx),'o');
        subplot(313), hold on, plot(idx,dt(idx),'o');
        if length(idx)<10
            title([sprintf('lost=%d',lost),' @ ', sprintf('%d,',idx)]);
        else
            title([sprintf('lost=%d',lost),' @ ', sprintf('%d,',idx(1:9)),'...']);
        end
    end

