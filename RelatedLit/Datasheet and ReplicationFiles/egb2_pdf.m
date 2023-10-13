function f=egb2_pdf(x,param)
%pdf of the egb2
%Y = egb2_pdf(x,param = [delta,sigma,p,q]) returns the value of the EGB2 pdf with parameters 
%delta,sigma,p,q evaluated at x.  Note that delta,sigma,p,q must be scalars.

d = param(1);
s = param(2);
p = param(3);
q = param(4);

if p > 0 && q > 0
    f = (exp((p.*(x-d))./s))./(abs(s).*beta(p,q).*((1+(exp((x-d)./s))).^(p+q)));
else
    f = NaN;
end

