function r = egb2rnd(param,n)

if nargin == 1
    n = 1;
end
if isempty(n)
    n = 1;
end

prob = rand(n,1);

r = zeros(n,1);

bi1 = betaincinv(prob,param(3),param(4));
bi2 = betaincinv(1-prob,param(4),param(3));

r(prob < .9) = param(1) - param(2)*log((1-bi1(prob < .9))./bi1(prob < .9));
r(prob >= .9) = param(1) + param(2)*log((1-bi2(prob >= .9))./bi2(prob >= .9));

