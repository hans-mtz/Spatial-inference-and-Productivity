function f = FRQFillLagGroupMeanTT(VarOfInterest,csindex,tsindex,tt)

n = size(VarOfInterest,1);
f = zeros(n,1); 
for ii = 1:max(csindex)
    ind1 = find(csindex == ii & tsindex == tt);
    ind0 = find(csindex == ii & tsindex == tt-1);
    if ~isempty(ind1) && ~isempty(ind0)
        f(ind1) = mean(VarOfInterest(ind0));
    end
end
