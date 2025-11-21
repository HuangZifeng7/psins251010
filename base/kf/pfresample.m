function idx = pfresample(w)
% PF weight resample
    len = length(w);
    cdf = cumsum(w/sum(w));
    uni = rand(len,1);  idx = ones(len,1);
    for k=1:len
        unik = uni(k);
        for m=2:len
            if cdf(m-1)<unik && unik<=cdf(m)
                idx(k)=m; break;
            end
        end
    end
    return;

%% test

len=100; w = abs(randn(len,1)); idx = pfresample(w);  myfig, plot([w/sum(w)*len^2, idx]);

w = [1 2 1 20 3 2 3 2 19 3]';  w=w/sum(w);
a = cumsum(w);  u = rand(10,1);
myfig, plot([w,a,u],'-o','linewidth',2); grid on; legend('w(i)','a(j)','\xi(n)');
xlabel('particles-i'); ylabel('weight / cumsum(w) / rand'); xlim([1,10]);

