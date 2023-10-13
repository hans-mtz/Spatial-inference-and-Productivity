function err = FRQSim_GenErr(distparams,Sigma,lagparams,indX,indT)
% lagparams are parameters for lagged group means with groups defined by
% each column of indX cross indT.  indT is the time index.  Need to have
% coefficients in lagparams line up with columns of indX.

n = numel(indT);
kX = size(indX,2);
err = zeros(n,1);

%% t = 1
S1 = Sigma{1};
n1 = size(S1,1);
v1 = (chol(S1)')*egb2rnd(distparams,n1);
err(indT == 1) = v1;

%% t > 1

for tt = 2:max(indT)
    Stt = Sigma{tt};
    ntt = size(Stt,1);
    vtt = (chol(Stt)')*egb2rnd(distparams,ntt);
    
    % Need to construct appropriate functions of lagged errors
    
    lagErr = zeros(n,kX);
    for ii = 1:kX
        lagErr(:,ii) = FRQFillLagGroupMeanTT(err,indX(:,ii),indT,tt);
    end
    err(indT == tt) = lagErr(indT == tt,:)*lagparams + vtt;

end

