% Final set of simulation results for JAR review - trying firm x year two
% way clustering per referee request

clear ;
warning('off','MATLAB:rankDeficientMatrix');


% try %#ok<TRYNC>
%     parpool();
% end

MainDir = 'C:\Users\chansen1\Dropbox\JARReview\Example_VerdiInvestment' ;
cd(MainDir);
yDir = 'C:\Users\chansen1\Dropbox\JARReview\Example_VerdiInvestment\SimulatedOutcomes' ;  

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
FxT = groupcross(f_ind,t_ind);
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

% standard errors - two-way clustered
sefirmtime = zeros(nSim,kx);  % baseline, TW by firm and time
nfirmtime = 0; % number not positive definite, state and time

seSxTfirmtime = zeros(nSim,kx);  % state X time effects, TW by firm and time
nSxTfirmtime = 0;  % number not positive definite, state X time effects, TW by firm and time

seSIC1xTfirmtime = zeros(nSim,kx);  % sic1 X time effects, TW by firm and time
nSIC1xTfirmtime = 0;  % number not positive definite, sic1 X time effects, TW by firm and time

seSIC2xTfirmtime = zeros(nSim,kx);  % sic2 X time effects, TW by firm and time
nSIC2xTfirmtime = 0;  % number not positive definite, sic2 X time effects, TW by firm and time


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
     
    [~,VF] = cluster_se(X,e,M,f_ind,k);
    [~,VT] = cluster_se(X,e,M,t_ind,k);
    [~,VFT] = cluster_se(X,e,M,FxT,k);
    V = VF + VT - VFT;
    if min(eig(V)) > 0
        sefirmtime(ii,:) = sqrt(diag(V))';
    else
        [temp1,temp2] = eig(V);
        Dstar = max(temp2,0);
        sefirmtime(ii,:) = sqrt(diag(temp1*Dstar*temp1'))';
        nfirmtime = nfirmtime + 1;
    end
        
    % Firm and State X Time effects
    YFMSxT = y - DFMSxT*(DFMSxT\y);
    bSxTii = XFMSxT\YFMSxT;
    bSxT(ii,:) = bSxTii';
    eSxT = YFMSxT - XFMSxT*bSxTii;


    [~,VF] = cluster_se(XFMSxT,eSxT,MSxT,f_ind,kSxT);
    [~,VT] = cluster_se(XFMSxT,eSxT,MSxT,t_ind,kSxT);
    [~,VFT] = cluster_se(XFMSxT,eSxT,MSxT,FxT,kSxT);
    V = VF + VT - VFT;
    if min(eig(V)) > 0
        seSxTfirmtime(ii,:) = sqrt(diag(V))';
    else
        [temp1,temp2] = eig(V);
        Dstar = max(temp2,0);
        seSxTfirmtime(ii,:) = sqrt(diag(temp1*Dstar*temp1'))';
        nSxTfirmtime = nSxTfirmtime + 1;
    end
    
    % Firm and SIC1 X Time effects
    YFMSIC1xT = y - DFMSIC1xT*(DFMSIC1xT\y);
    bSIC1xTii = XFMSIC1xT\YFMSIC1xT;
    bSIC1xT(ii,:) = bSIC1xTii';
    eSIC1xT = YFMSIC1xT - XFMSIC1xT*bSIC1xTii;

    [~,VF] = cluster_se(XFMSIC1xT,eSIC1xT,MSIC1xT,f_ind,kSIC1xT);
    [~,VT] = cluster_se(XFMSIC1xT,eSIC1xT,MSIC1xT,t_ind,kSIC1xT);
    [~,VFT] = cluster_se(XFMSIC1xT,eSIC1xT,MSIC1xT,FxT,kSIC1xT);
    V = VF + VT - VFT;
    if min(eig(V)) > 0
        seSIC1xTfirmtime(ii,:) = sqrt(diag(V))';
    else
        [temp1,temp2] = eig(V);
        Dstar = max(temp2,0);
        seSIC1xTfirmtime(ii,:) = sqrt(diag(temp1*Dstar*temp1'))';
        nSIC1xTfirmtime = nSIC1xTfirmtime + 1;
    end
    
    % Firm and SIC2 X Time effects
    YFMSIC2xT = y - DFMSIC2xT*(DFMSIC2xT\y);
    bSIC2xTii = XFMSIC2xT\YFMSIC2xT;
    bSIC2xT(ii,:) = bSIC2xTii';
    eSIC2xT = YFMSIC2xT - XFMSIC2xT*bSIC2xTii;

    [~,VF] = cluster_se(XFMSIC2xT,eSIC2xT,MSIC2xT,f_ind,kSIC2xT);
    [~,VT] = cluster_se(XFMSIC2xT,eSIC2xT,MSIC2xT,t_ind,kSIC2xT);
    [~,VFT] = cluster_se(XFMSIC2xT,eSIC2xT,MSIC2xT,FxT,kSIC2xT);
    V = VF + VT - VFT;
    if min(eig(V)) > 0
        seSIC2xTfirmtime(ii,:) = sqrt(diag(V))';
    else
        [temp1,temp2] = eig(V);
        Dstar = max(temp2,0);
        seSIC2xTfirmtime(ii,:) = sqrt(diag(temp1*Dstar*temp1'))';
        nSIC2xTfirmtime = nSIC2xTfirmtime + 1;
    end
    
end

save FRQSim_FINAL_FIRMTIME ;
