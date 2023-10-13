function f = FRQFillLagGroupMean(VarOfInterest,csindex,tsindex)

n = size(VarOfInterest,1);
f = zeros(n,1); 
for ii = 1:max(csindex)
    for jj = 1:max(tsindex)
        ind1 = find(csindex == ii & tsindex == jj);
        ind0 = find(csindex == ii & tsindex == jj-1);
        if ~isempty(ind1) && ~isempty(ind0)
            f(ind1) = mean(VarOfInterest(ind0));            
        end
    end
end
