function pvals = FamaMacbethRand(b,h0)

b = b(~isnan(b(:,1)),:);
nb = size(b,1);

if nb < 12
    sets = 'setprod([-1,1],';
    for jj = 2:nb-1
        sets = strcat(sets,'[-1,1],');
    end
    sets = strcat(sets,'[-1,1]);');
    signs = eval(sets);
else
    signs = 2*(rand(1000,size(b,1)) < .5)-1; 
end

tests = zeros(size(signs,1),size(b,2));
for ii = 1:size(signs,1)
    tests(ii,:) = (signs(ii,:)*(b-ones(nb,1)*h0)/nb)./(std((signs(ii,:)'*ones(1,size(b,2))).*(b-ones(nb,1)*h0))/sqrt(nb));
end

tests = abs(tests);        
FM = abs((ones(1,nb)*(b-ones(nb,1)*h0)/nb)./(std((ones(1,nb)'*ones(1,size(b,2))).*(b-ones(nb,1)*h0))/sqrt(nb)));

pvals = mean(ones(size(tests,1),1)*FM <= tests);


