kse = -4:.1:4;
indSize = kse == 0;
nse = numel(kse);
s = std(b);

warning('off','MATLAB:rankDeficientMatrix');
crsS = zeros(nSim,kx);
crsSIC = zeros(nSim,kx);
p4 = zeros(nSim,kx);
p2 = zeros(nSim,kx);
pCRSS = zeros(nse,kx);
pCRSSIC = zeros(nse,kx);
pCRS4 = zeros(nse,kx);
pCRS2 = zeros(nse,kx);
for jj = 1:nse
    disp(jj);
    balt = bTest3' + kse(jj)*s;
    
    for ii = 1:nSim
        disp(ii);
        filename = strcat(yDir,'\ySim',num2str(ii),'.mat');
        load(filename);
        YFMFxT4 = y - DFMFxT4*(DFMFxT4\y);
        YFMFxT2 = y - DFMFxT2*(DFMFxT2\y);
        YFMSIC1xT = y - DFMSIC1xT*(DFMSIC1xT\y);
        YFMSxT = y - DFMSxT*(DFMSxT\y);
        
        % T4 groups
        [~,~,tmp3] = FamaMacbeth(XFMFxT4,YFMFxT4,tind4);
        p4(ii,:) = FamaMacbethRand(tmp3,balt);
        
        % T2 groups
        [~,~,tmp3] = FamaMacbeth(XFMFxT2,YFMFxT2,tind2);
        p2(ii,:) = FamaMacbethRand(tmp3,balt);
        
        % size groups
        [~,~,tmp3] = FamaMacbeth(XFMSIC1xT,YFMSIC1xT,sic1);
        crsSIC(ii,:) = FamaMacbethRand(tmp3,balt);
        
        % state groups
        [~,~,tmp3] = FamaMacbeth(XFMSxT,YFMSxT,state);
        crsS(ii,:) = FamaMacbethRand(tmp3,balt);
    end
    pCRS4(jj,:) = mean(p4 < .05);
    pCRS2(jj,:) = mean(p2 < .05);
    pCRSSIC(jj,:) = mean(crsSIC < .05);
    pCRSS(jj,:) = mean(crsS < .05);
end
