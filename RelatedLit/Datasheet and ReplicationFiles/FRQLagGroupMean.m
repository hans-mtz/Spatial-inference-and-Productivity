function f = FRQLagGroupMean(VarOfInterest,csindex,tsindex,p)

if nargin == 3
    p = 1;
end
if isempty(p)
    p = 1;
end

n = size(VarOfInterest,1);
f = zeros(n,1); 
for ii = 1:max(csindex)
    for jj = 1:max(tsindex)
        ind1 = find(csindex == ii & tsindex == jj);
        ind0 = find(csindex == ii & tsindex == jj-p);
        if ~isempty(ind1) && ~isempty(ind0)
            f(ind1) = mean(VarOfInterest(ind0)); 
        else
            f(ind1) = NaN;
        end
    end
end
