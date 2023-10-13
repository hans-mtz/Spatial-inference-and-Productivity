function Iplus = index_mbb( index, bsize, T, J )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

iconsecutive = repmat((1:bsize-1)',1,ceil(T/bsize));

imat= repmat(index,bsize-1,1) + iconsecutive;

Iaux= [index ; imat];

n1=size(Iaux,1); 
n2=size(Iaux,2);

Iaux_=reshape(Iaux,[n1*n2 1]);

Iaux_=Iaux_(1:T,1);

Iplus=repmat(Iaux_,J,1);


end

