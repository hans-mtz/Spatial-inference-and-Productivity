% Final set of simulation results for JAR review - Wild Cluster Bootstrap
% with Mammen weights.  Only look at one-way clustering.

clear ;
warning('off','MATLAB:rankDeficientMatrix');


try %#ok<TRYNC>
    parpool();
end
pctRunOnAll warning('off','MATLAB:rankDeficientMatrix');

% MainDir = 'C:\Users\chansen1\Dropbox\JARReview\Example_VerdiInvestment' ;
% cd(MainDir);
% yDir = 'C:\Users\chansen1\Dropbox\JARReview\Example_VerdiInvestment\SimulatedOutcomes' ;  
MainDir = 'D:\Dropbox\JARReview\Example_VerdiInvestment' ;
cd(MainDir);
yDir = 'D:\Dropbox\JARReview\Example_VerdiInvestment\SimulatedOutcomes' ;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load x data, index variables, and true regression coefficient values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load XCatIndices ;   % Indices of categories used to form lagged variables
                     % indX = [f_ind state sic2 sizecat cashcat rescat frqcat intcat]
load TimeIndices ;   % Time indices - t_ind
load FirmIndices ;   % Firm indices - f_ind
load xMat ;          % Design matrix from actual data
                     % x = [res ind frq frqint cash tq mve lev age]
load olsparams ;     % OLS estimates of parameters on x variables
                     % bTest3 = coefficients on x used to generate y
load BCVResidualRegressions ; % Results from residual analysis regressions 
                     % from baseline two-way FE model.  Included variables are
                     % pac1 spac1 - time series model with three lags of own residual and cross-sectional averages
                     % acf1 acf2 acf3 sacf1 sacf2 sacf3 - "autocorrelations" as above
                     % spcoef se_spcoef - coefficients from spatial cross-products regression
                     % sptest - per variable test statistics from spatial cross-products regression
                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data preprocessing - Partial FE out of design matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Baseline - firm and year
Di = sparse(dummyvar(f_ind));
Dt = sparse(dummyvar(t_ind));
D = [Di Dt];
X = x - D*(D\x);
M = inv(X'*X);

n = size(X,1);
kx = size(X,2);
k = kx+trace(D\D);

% firm and state cross year
state = indX(:,2);
SxT = groupcross(state,t_ind);
DSxT = sparse(dummyvar(SxT));
DFMSxT = [Di DSxT];
XFMSxT = x - DFMSxT*(DFMSxT\x);
MSxT = inv(XFMSxT'*XFMSxT);
kSxT = kx + trace(DFMSxT\DFMSxT);

% firm and 1 digit SIC cross year
sic2 = indX(:,3);
sic1 = recode(floor(sic2/10)+1);
SIC1xT = groupcross(sic1,t_ind);
DSIC1xT = sparse(dummyvar(SIC1xT));
DFMSIC1xT = [Di DSIC1xT];
XFMSIC1xT = x - DFMSIC1xT*(DFMSIC1xT\x);
MSIC1xT = inv(XFMSIC1xT'*XFMSIC1xT);
kSIC1xT = kx + trace(DFMSIC1xT\DFMSIC1xT);

% firm and 2 digit SIC cross year
sic2 = recode(sic2);
SIC2xT = groupcross(sic2,t_ind);
DSIC2xT = sparse(dummyvar(SIC2xT));
DFMSIC2xT = [Di DSIC2xT];
XFMSIC2xT = x - DFMSIC2xT*(DFMSIC2xT\x);
MSIC2xT = inv(XFMSIC2xT'*XFMSIC2xT);
kSIC2xT = kx + trace(DFMSIC2xT\DFMSIC2xT);

% firm and size groups cross year
sizecat = indX(:,4);
SIZExT = groupcross(sizecat,t_ind);
DSIZExT = sparse(dummyvar(SIZExT));
DFMSIZExT = [Di DSIZExT];
XFMSIZExT = x - DFMSIZExT*(DFMSIZExT\x);
MSIZExT = inv(XFMSIZExT'*XFMSIZExT);
kSIZExT = kx + trace(DFMSIZExT\DFMSIZExT);

% firm cross eight year time blocks and year
tind8 = ceil(t_ind/8);
tind8(tind8 == 3) = 2;
FxT8 = groupcross(f_ind,tind8);
for ii = 1:max(FxT8)
    DFxT8(:,ii) = sparse(FxT8 == ii); %#ok<SPRIX,SAGROW>
end
DFMFxT8 = [DFxT8 Dt];
XFMFxT8 = x - DFMFxT8*(DFMFxT8\x);
MFxT8 = inv(XFMFxT8'*XFMFxT8);
kFxT8 = kx + trace(DFMFxT8\DFMFxT8);

% firm cross six year time blocks and year
tind6 = ceil(t_ind/6);
FxT6 = groupcross(f_ind,tind6);
for ii = 1:max(FxT6)
    DFxT6(:,ii) = sparse(FxT6 == ii); %#ok<SPRIX,SAGROW>
end
DFMFxT6 = [DFxT6 Dt];
XFMFxT6 = x - DFMFxT6*(DFMFxT6\x);
MFxT6 = inv(XFMFxT6'*XFMFxT6);
kFxT6 = kx + trace(DFMFxT6\DFMFxT6);

% firm cross four year time blocks and year
tind4 = ceil(t_ind/4);
tind4(tind4 == 5) = 4;
FxT4 = groupcross(f_ind,tind4);
for ii = 1:max(FxT4)
    DFxT4(:,ii) = sparse(FxT4 == ii); %#ok<SPRIX,SAGROW>
end
DFMFxT4 = [DFxT4 Dt];
XFMFxT4 = x - DFMFxT4*(DFMFxT4\x);
MFxT4 = inv(XFMFxT4'*XFMFxT4);
kFxT4 = kx + trace(DFMFxT4\DFMFxT4);

% firm cross two year time blocks and year
tind2 = ceil(t_ind/2);
tind2(tind2 == 9) = 8;
FxT2 = groupcross(f_ind,tind2);
for ii = 1:max(FxT2)
    DFxT2(:,ii) = sparse(FxT2 == ii); %#ok<SPRIX,SAGROW>
end
DFMFxT2 = [DFxT2 Dt];
XFMFxT2 = x - DFMFxT2*(DFMFxT2\x);
MFxT2 = inv(XFMFxT2'*XFMFxT2);
kFxT2 = kx + trace(DFMFxT2\DFMFxT2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Predefine variables for simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of simulation replications
nSim = 1000 ;

% Point estimates - FE
b = zeros(nSim,kx); % Baseline estimator (Firm and Time fixed effects)
bSxT = zeros(nSim,kx); % Firm and State X Time effects
bSIC1xT = zeros(nSim,kx); % Firm and SIC1 X Time effects
bSIC2xT = zeros(nSim,kx); % Firm and SIC2 X Time effects
bSIZExT = zeros(nSim,kx); % Firm and Size X Time effects
bFxT8 = zeros(nSim,kx); % Firm X 8 year time block and Time effects
bFxT6 = zeros(nSim,kx); % Firm X 6 year time block and Time effects
bFxT4 = zeros(nSim,kx); % Firm X 4 year time block and Time effects
bFxT2 = zeros(nSim,kx); % Firm X 2 year time block and Time effects

% Standard errors - clustered
sefirm = zeros(nSim,kx);  % baseline cluster by firm
sestate = zeros(nSim,kx);  % baseline cluster by state
sesic1 = zeros(nSim,kx);  % baseline cluster by sic1
sesic2 = zeros(nSim,kx);  % baseline cluster by sic2
sesize = zeros(nSim,kx);  % baseline cluster by size
seT8 = zeros(nSim,kx);  % baseline cluster by 8 year time block
seT6 = zeros(nSim,kx);  % baseline cluster by 6 year time block
seT4 = zeros(nSim,kx);  % baseline cluster by 4 year time block
seT2 = zeros(nSim,kx);  % baseline cluster by 2 year time block

seSxTstate = zeros(nSim,kx);  % state X time effects, cluster by state

seSIC1xTsic1 = zeros(nSim,kx);  % sic1 X time effects, cluster by sic1

seSIC2xTsic1 = zeros(nSim,kx);  % sic2 X time effects, cluster by sic1
seSIC2xTsic2 = zeros(nSim,kx);  % sic2 X time effects, cluster by sic2

seSIZExTsize = zeros(nSim,kx);  % size X time effects, cluster by size

seFxT8 = zeros(nSim,kx);  % t8 X firm effect, cluster by T8

seFxT6 = zeros(nSim,kx);  % t6 X firm effects, cluster by T6

seFxT4 = zeros(nSim,kx);  % t4 X firm effects, cluster by T4

seFxT2 = zeros(nSim,kx);  % t2 X firm effects, cluster by T2

% Bootstrap critical values
cvfirm = zeros(nSim,kx);  % baseline cluster by firm
cvstate = zeros(nSim,kx);  % baseline cluster by state
cvsic1 = zeros(nSim,kx);  % baseline cluster by sic1
cvsic2 = zeros(nSim,kx);  % baseline cluster by sic2
cvsize = zeros(nSim,kx);  % baseline cluster by size
cvT8 = zeros(nSim,kx);  % baseline cluster by 8 year time block
cvT6 = zeros(nSim,kx);  % baseline cluster by 6 year time block
cvT4 = zeros(nSim,kx);  % baseline cluster by 4 year time block
cvT2 = zeros(nSim,kx);  % baseline cluster by 2 year time block

cvSxTstate = zeros(nSim,kx);  % state X time effects, cluster by state

cvSIC1xTsic1 = zeros(nSim,kx);  % sic1 X time effects, cluster by sic1

cvSIC2xTsic1 = zeros(nSim,kx);  % sic2 X time effects, cluster by sic1
cvSIC2xTsic2 = zeros(nSim,kx);  % sic2 X time effects, cluster by sic2

cvSIZExTsize = zeros(nSim,kx);  % size X time effects, cluster by size

cvFxT8 = zeros(nSim,kx);  % t8 X firm effect, cluster by T8

cvFxT6 = zeros(nSim,kx);  % t6 X firm effects, cluster by T6

cvFxT4 = zeros(nSim,kx);  % t4 X firm effects, cluster by T4

cvFxT2 = zeros(nSim,kx);  % t2 X firm effects, cluster by T2

nb = 500;  % Number of bootstrap replications

%% Main simulation step
for ii = 1:nSim
    disp(ii);
    
    % Load presimulated outcome data
    filename = strcat(yDir,'\ySim',num2str(ii),'.mat');
    load(filename);

    % OLS/FE estimates

    % Baseline estimate
    Y = y - D*(D\y);
    bii = X\Y;
    b(ii,:) = bii';
    e = Y-X*bii;
 
    sefirm(ii,:) = cluster_se(X,e,M,f_ind,k)';
    sestate(ii,:) = cluster_se(X,e,M,state,k)';
    sesic1(ii,:) = cluster_se(X,e,M,sic1,k)';
    sesic2(ii,:) = cluster_se(X,e,M,sic2,k)';
    sesize(ii,:) = cluster_se(X,e,M,sizecat,k)';
    seT8(ii,:) = cluster_se(X,e,M,tind8,k)';
    seT6(ii,:) = cluster_se(X,e,M,tind6,k)';
    seT4(ii,:) = cluster_se(X,e,M,tind4,k)';
    seT2(ii,:) = cluster_se(X,e,M,tind2,k)';
    
    cvfirm(ii,:) = cluster_wildboot(X,e,M,D,f_ind,k,bii,nb);
    cvstate(ii,:) = cluster_wildboot(X,e,M,D,state,k,bii,nb);
    cvsic1(ii,:) = cluster_wildboot(X,e,M,D,sic1,k,bii,nb);
    cvsic2(ii,:) = cluster_wildboot(X,e,M,D,sic2,k,bii,nb);
    cvsize(ii,:) = cluster_wildboot(X,e,M,D,sizecat,k,bii,nb);
    cvT8(ii,:) = cluster_wildboot(X,e,M,D,tind8,k,bii,nb);
    cvT6(ii,:) = cluster_wildboot(X,e,M,D,tind6,k,bii,nb);
    cvT4(ii,:) = cluster_wildboot(X,e,M,D,tind4,k,bii,nb);
    cvT2(ii,:) = cluster_wildboot(X,e,M,D,tind2,k,bii,nb);

    % Firm and State X Time effects
    YFMSxT = y - DFMSxT*(DFMSxT\y);
    bSxTii = XFMSxT\YFMSxT;
    bSxT(ii,:) = bSxTii';
    eSxT = YFMSxT - XFMSxT*bSxTii;

    seSxTstate(ii,:) = cluster_se(XFMSxT,eSxT,MSxT,state,kSxT)';

    cvSxTstate(ii,:) = cluster_wildboot(XFMSxT,eSxT,MSxT,DFMSxT,state,kSxT,bSxTii,nb);
    
    % Firm and SIC1 X Time effects
    YFMSIC1xT = y - DFMSIC1xT*(DFMSIC1xT\y);
    bSIC1xTii = XFMSIC1xT\YFMSIC1xT;
    bSIC1xT(ii,:) = bSIC1xTii';
    eSIC1xT = YFMSIC1xT - XFMSIC1xT*bSIC1xTii;
    
    seSIC1xTsic1(ii,:) = cluster_se(XFMSIC1xT,eSIC1xT,MSIC1xT,sic1,kSIC1xT)';
    
    cvSIC1xTsic1(ii,:) = cluster_wildboot(XFMSIC1xT,eSIC1xT,MSIC1xT,DFMSIC1xT,sic1,kSIC1xT,bSIC1xTii,nb);
    
    % Firm and SIC2 X Time effects
    YFMSIC2xT = y - DFMSIC2xT*(DFMSIC2xT\y);
    bSIC2xTii = XFMSIC2xT\YFMSIC2xT;
    bSIC2xT(ii,:) = bSIC2xTii';
    eSIC2xT = YFMSIC2xT - XFMSIC2xT*bSIC2xTii;

    seSIC2xTsic1(ii,:) = cluster_se(XFMSIC2xT,eSIC2xT,MSIC2xT,sic1,kSIC2xT)';
    seSIC2xTsic2(ii,:) = cluster_se(XFMSIC2xT,eSIC2xT,MSIC2xT,sic2,kSIC2xT)';
        
    cvSIC2xTsic1(ii,:) = cluster_wildboot(XFMSIC2xT,eSIC2xT,MSIC2xT,DFMSIC2xT,sic1,kSIC2xT,bSIC2xTii,nb);
    cvSIC2xTsic2(ii,:) = cluster_wildboot(XFMSIC2xT,eSIC2xT,MSIC2xT,DFMSIC2xT,sic2,kSIC2xT,bSIC2xTii,nb);
        
    % Firm and Size X Time effects
    YFMSIZExT = y - DFMSIZExT*(DFMSIZExT\y);
    bSIZExTii = XFMSIZExT\YFMSIZExT;
    bSIZExT(ii,:) = bSIZExTii';
    eSIZExT = YFMSIZExT - XFMSIZExT*bSIZExTii;

    seSIZExTsize(ii,:) = cluster_se(XFMSIZExT,eSIZExT,MSIZExT,sizecat,kSIZExT)';

    cvSIZExTsize(ii,:) = cluster_wildboot(XFMSIZExT,eSIZExT,MSIZExT,DFMSIZExT,sizecat,kSIZExT,bSIZExTii,nb);

    % Firm X 8 year block and Time effects
    YFMFxT8 = y - DFMFxT8*(DFMFxT8\y);
    bFxT8ii = XFMFxT8\YFMFxT8;
    bFxT8(ii,:) = bFxT8ii;
    eFxT8 = YFMFxT8 - XFMFxT8*bFxT8ii;
    
    seFxT8(ii,:) = cluster_se(XFMFxT8,eFxT8,MFxT8,tind8,kFxT8)';
    
    cvFxT8(ii,:) = cluster_wildboot(XFMFxT8,eFxT8,MFxT8,DFMFxT8,tind8,kFxT8,bFxT8ii,nb);
    
    % Firm X 6 year block and Time effects
    YFMFxT6 = y - DFMFxT6*(DFMFxT6\y);
    bFxT6ii = XFMFxT6\YFMFxT6;
    bFxT6(ii,:) = bFxT6ii;
    eFxT6 = YFMFxT6 - XFMFxT6*bFxT6ii;
    
    seFxT6(ii,:) = cluster_se(XFMFxT6,eFxT6,MFxT6,tind6,kFxT6)';

    cvFxT6(ii,:) = cluster_wildboot(XFMFxT6,eFxT6,MFxT6,DFMFxT6,tind6,kFxT6,bFxT6ii,nb);

    % Firm X 4 year block and Time effects
    YFMFxT4 = y - DFMFxT4*(DFMFxT4\y);
    bFxT4ii = XFMFxT4\YFMFxT4;
    bFxT4(ii,:) = bFxT4ii;
    eFxT4 = YFMFxT4 - XFMFxT4*bFxT4ii;
    
    seFxT4(ii,:) = cluster_se(XFMFxT4,eFxT4,MFxT4,tind4,kFxT4)';

    cvFxT4(ii,:) = cluster_wildboot(XFMFxT4,eFxT4,MFxT4,DFMFxT4,tind4,kFxT4,bFxT4ii,nb);

    % Firm X 2 year block and Time effects
    YFMFxT2 = y - DFMFxT2*(DFMFxT2\y);
    bFxT2ii = XFMFxT2\YFMFxT2;
    bFxT2(ii,:) = bFxT2ii;
    eFxT2 = YFMFxT2 - XFMFxT2*bFxT2ii;

    seFxT2(ii,:) = cluster_se(XFMFxT2,eFxT2,MFxT2,tind2,kFxT2)';
        
    cvFxT2(ii,:) = cluster_wildboot(XFMFxT2,eFxT2,MFxT2,DFMFxT2,tind2,kFxT2,bFxT2ii,nb);

end

save FRQSim_FINAL_BOOT ;
