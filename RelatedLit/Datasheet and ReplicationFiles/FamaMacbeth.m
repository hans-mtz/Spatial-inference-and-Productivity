function [b,se,btemp] = FamaMacbeth(X,Y,index)

p = max(index);
k = size(X,2);

btemp = zeros(max(index),k);
pmiss = 0;
for ii = 1:max(index)
    fii = index == ii;
    if sum(fii) > k
        btemp(ii,:) = (X(fii,:)\Y(fii))';
    else
        btemp(ii,:) = NaN*ones(1,k);
        pmiss = pmiss+1;
    end
end
b = nanmean(btemp);
se = nanstd(btemp)/sqrt(p-pmiss);
