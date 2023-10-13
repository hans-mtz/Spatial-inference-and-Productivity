% Final set of simulation results for JAR review - size, no Canay, Romano,
% Shaikh - need to do that in a different file - compute size and power
% simultaneously 
% Tabulate size and power elsewhere for all these estimators (and other
% stats, e.g. RMSE)

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
                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Look at results from residual regressions to benchmark cluster structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nrej = sum(sptest > [chi2inv(.95,1)*ones(3,1) ; chi2inv(.95,4)*ones(9,1)]*ones(1,1000)); % How many variables are "important" in space
nrejacf1 = sum(abs(acf1 - [-.1129;zeros(12,1)]*ones(1,1000))./sacf1 > tinv(.975,13));    % How many variables are "important" in time (first order)
nrejacf2 = sum(abs(acf2 - [-.1129;zeros(12,1)]*ones(1,1000))./sacf2 > tinv(.975,13));    % How many variables are "important" in time (second order)
nrejacf3 = sum(abs(acf3 - [-.1129;zeros(12,1)]*ones(1,1000))./sacf3 > tinv(.975,13));    % How many variables are "important" in time (third order)

sum(nrejacf1 == 1 & nrej == 1)  % Number of times only one spatial dimension loads
                               % Never!  Will make splitting on
                               % cross-section hard.  Only consider time
                               % series splits.
                               
Trep = 2*(nrejacf1 >= 1 & nrejacf2 == 0 & nrejacf3 == 0)...
    + 4*(nrejacf2 >= 1 & nrejacf3 == 0) + 6*(nrejacf3 >= 1);     % How large to make time blocks based on autocorrelation evidence

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check design homogeneity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Form partialed out matrices
pX1 = zeros(size(X)); % firm and time
pX2 = zeros(size(X)); % firm and state cross year
pX3 = zeros(size(X)); % firm and 1 digit SIC cross year
pX4 = zeros(size(X)); % firm and 2 digit SIC cross year
pX5 = zeros(size(X)); % firm and size groups cross year
pX6 = zeros(size(X)); % firm cross eight year time blocks and year
pX7 = zeros(size(X)); % firm cross six year time blocks and year
pX8 = zeros(size(X)); % firm cross four year time blocks and year
pX9 = zeros(size(X)); % firm cross two year time blocks and year
for kk = 1:kx
    pX1(:,kk) = X(:,kk) - X(:,setdiff(1:kx,kk))*(X(:,setdiff(1:kx,kk))\X(:,kk));
    pX2(:,kk) = XFMSxT(:,kk) - XFMSxT(:,setdiff(1:kx,kk))*(XFMSxT(:,setdiff(1:kx,kk))\X(:,kk));
    pX3(:,kk) = XFMSIC1xT(:,kk) - XFMSIC1xT(:,setdiff(1:kx,kk))*(XFMSIC1xT(:,setdiff(1:kx,kk))\XFMSIC1xT(:,kk));
    pX4(:,kk) = XFMSIC2xT(:,kk) - XFMSIC2xT(:,setdiff(1:kx,kk))*(XFMSIC2xT(:,setdiff(1:kx,kk))\XFMSIC2xT(:,kk));
    pX5(:,kk) = XFMSIZExT(:,kk) - XFMSIZExT(:,setdiff(1:kx,kk))*(XFMSIZExT(:,setdiff(1:kx,kk))\XFMSIZExT(:,kk));
    pX6(:,kk) = XFMFxT8(:,kk) - XFMFxT8(:,setdiff(1:kx,kk))*(XFMFxT8(:,setdiff(1:kx,kk))\XFMFxT8(:,kk));    
    pX7(:,kk) = XFMFxT6(:,kk) - XFMFxT6(:,setdiff(1:kx,kk))*(XFMFxT6(:,setdiff(1:kx,kk))\XFMFxT6(:,kk));    
    pX8(:,kk) = XFMFxT4(:,kk) - XFMFxT4(:,setdiff(1:kx,kk))*(XFMFxT4(:,setdiff(1:kx,kk))\XFMFxT4(:,kk));    
    pX9(:,kk) = XFMFxT2(:,kk) - XFMFxT2(:,setdiff(1:kx,kk))*(XFMFxT2(:,setdiff(1:kx,kk))\XFMFxT2(:,kk));    
end    

% Test equality of squares across groups assuming iid
disp('State');
Dstate = sparse(dummyvar(state));
Dstate = [ones(n,1) Dstate(:,2:end)];
iDstate = inv(Dstate'*Dstate);

bb1state = Dstate\(pX1.^2);
vv1state = pX1.^2 - Dstate*bb1state;
s21state = mean(vv1state.^2);
x21state = diag(bb1state(2:end,:)'*inv(iDstate(2:end,2:end))*bb1state(2:end,:))./s21state'; %#ok<*MINV>
pstate(:,1) = 1-chi2cdf(x21state,size(Dstate,2));

bb2state = Dstate\(pX2.^2);
vv2state = pX2.^2 - Dstate*bb2state;
s22state = mean(vv2state.^2);
x22state = diag(bb2state(2:end,:)'*inv(iDstate(2:end,2:end))*bb2state(2:end,:))./s22state';
pstate(:,2) = 1-chi2cdf(x22state,size(Dstate,2));

disp('sic1');
Dsic1 = sparse(dummyvar(sic1));
Dsic1 = [ones(n,1) Dsic1(:,2:end)];
iDsic1 = inv(Dsic1'*Dsic1);

bb1sic1 = Dsic1\(pX1.^2);
vv1sic1 = pX1.^2 - Dsic1*bb1sic1;
s21sic1 = mean(vv1sic1.^2);
x21sic1 = diag(bb1sic1(2:end,:)'*inv(iDsic1(2:end,2:end))*bb1sic1(2:end,:))./s21sic1'; %#ok<*MINV>
psic1(:,1) = 1-chi2cdf(x21sic1,size(Dsic1,2)); %#ok<*MNEFF,*NOPTS>

bb3sic1 = Dsic1\(pX3.^2);
vv3sic1 = pX3.^2 - Dsic1*bb3sic1;
s23sic1 = mean(vv3sic1.^2);
x23sic1 = diag(bb3sic1(2:end,:)'*inv(iDsic1(2:end,2:end))*bb3sic1(2:end,:))./s23sic1';
psic1(:,2) = 1-chi2cdf(x23sic1,size(Dsic1,2));

bb4sic1 = Dsic1\(pX4.^2);
vv4sic1 = pX4.^2 - Dsic1*bb4sic1;
s24sic1 = mean(vv4sic1.^2);
x24sic1 = diag(bb4sic1(2:end,:)'*inv(iDsic1(2:end,2:end))*bb4sic1(2:end,:))./s24sic1';
psic1(:,3) = 1-chi2cdf(x24sic1,size(Dsic1,2));

disp('sic2');
Dsic2 = sparse(dummyvar(sic2));
Dsic2 = [ones(n,1) Dsic2(:,2:end)];
iDsic2 = inv(Dsic2'*Dsic2);

bb1sic2 = Dsic2\(pX1.^2);
vv1sic2 = pX1.^2 - Dsic2*bb1sic2;
s21sic2 = mean(vv1sic2.^2);
x21sic2 = diag(bb1sic2(2:end,:)'*inv(iDsic2(2:end,2:end))*bb1sic2(2:end,:))./s21sic2'; %#ok<*MINV>
psic2(:,1) = 1-chi2cdf(x21sic2,size(Dsic2,2)); %#ok<*MNEFF,*NOPTS>

bb3sic2 = Dsic2\(pX3.^2);
vv3sic2 = pX3.^2 - Dsic2*bb3sic2;
s23sic2 = mean(vv3sic2.^2);
x23sic2 = diag(bb3sic2(2:end,:)'*inv(iDsic2(2:end,2:end))*bb3sic2(2:end,:))./s23sic2';
psic2(:,2) = 1-chi2cdf(x23sic2,size(Dsic2,2));

bb4sic2 = Dsic2\(pX4.^2);
vv4sic2 = pX4.^2 - Dsic2*bb4sic2;
s24sic2 = mean(vv4sic2.^2);
x24sic2 = diag(bb4sic2(2:end,:)'*inv(iDsic2(2:end,2:end))*bb4sic2(2:end,:))./s24sic2';
psic2(:,3) = 1-chi2cdf(x24sic2,size(Dsic2,2));

disp('sizecat');
Dsizecat = sparse(dummyvar(sizecat));
Dsizecat = [ones(n,1) Dsizecat(:,2:end)];
iDsizecat = inv(Dsizecat'*Dsizecat);

bb1sizecat = Dsizecat\(pX1.^2);
vv1sizecat = pX1.^2 - Dsizecat*bb1sizecat;
s21sizecat = mean(vv1sizecat.^2);
x21sizecat = diag(bb1sizecat(2:end,:)'*inv(iDsizecat(2:end,2:end))*bb1sizecat(2:end,:))./s21sizecat'; %#ok<*MINV>
psizecat(:,1) = 1-chi2cdf(x21sizecat,size(Dsizecat,2)); %#ok<*MNEFF,*NOPTS>

bb5sizecat = Dsizecat\(pX5.^2);
vv5sizecat = pX5.^2 - Dsizecat*bb5sizecat;
s25sizecat = mean(vv5sizecat.^2);
x25sizecat = diag(bb5sizecat(2:end,:)'*inv(iDsizecat(2:end,2:end))*bb5sizecat(2:end,:))./s25sizecat';
psizecat(:,2) = 1-chi2cdf(x25sizecat,size(Dsizecat,2));

disp('tind8');
Dtind8 = sparse(dummyvar(tind8));
Dtind8 = [ones(n,1) Dtind8(:,2:end)];
iDtind8 = inv(Dtind8'*Dtind8);

bb1tind8 = Dtind8\(pX1.^2);
vv1tind8 = pX1.^2 - Dtind8*bb1tind8;
s21tind8 = mean(vv1tind8.^2);
x21tind8 = diag(bb1tind8(2:end,:)'*inv(iDtind8(2:end,2:end))*bb1tind8(2:end,:))./s21tind8'; %#ok<*MINV>
ptind8(:,1) = 1-chi2cdf(x21tind8,size(Dtind8,2)); %#ok<*MNEFF,*NOPTS>

bb6tind8 = Dtind8\(pX6.^2);
vv6tind8 = pX6.^2 - Dtind8*bb6tind8;
s26tind8 = mean(vv6tind8.^2);
x26tind8 = diag(bb6tind8(2:end,:)'*inv(iDtind8(2:end,2:end))*bb6tind8(2:end,:))./s26tind8';
ptind8(:,2) = 1-chi2cdf(x26tind8,size(Dtind8,2));

disp('tind6');
Dtind6 = sparse(dummyvar(tind6));
Dtind6 = [ones(n,1) Dtind6(:,2:end)];
iDtind6 = inv(Dtind6'*Dtind6);

bb1tind6 = Dtind6\(pX1.^2);
vv1tind6 = pX1.^2 - Dtind6*bb1tind6;
s21tind6 = mean(vv1tind6.^2);
x21tind6 = diag(bb1tind6(2:end,:)'*inv(iDtind6(2:end,2:end))*bb1tind6(2:end,:))./s21tind6'; %#ok<*MINV>
ptind6(:,1) = 1-chi2cdf(x21tind6,size(Dtind6,2)); %#ok<*MNEFF,*NOPTS>

bb7tind6 = Dtind6\(pX7.^2);
vv7tind6 = pX7.^2 - Dtind6*bb7tind6;
s27tind6 = mean(vv7tind6.^2);
x27tind6 = diag(bb7tind6(2:end,:)'*inv(iDtind6(2:end,2:end))*bb7tind6(2:end,:))./s27tind6';
ptind6(:,2) = 1-chi2cdf(x27tind6,size(Dtind6,2));

disp('tind4');
Dtind4 = sparse(dummyvar(tind4));
Dtind4 = [ones(n,1) Dtind4(:,2:end)];
iDtind4 = inv(Dtind4'*Dtind4);

bb1tind4 = Dtind4\(pX1.^2);
vv1tind4 = pX1.^2 - Dtind4*bb1tind4;
s21tind4 = mean(vv1tind4.^2);
x21tind4 = diag(bb1tind4(2:end,:)'*inv(iDtind4(2:end,2:end))*bb1tind4(2:end,:))./s21tind4'; %#ok<*MINV>
ptind4(:,1) = 1-chi2cdf(x21tind4,size(Dtind4,2)); %#ok<*MNEFF,*NOPTS>

bb8tind4 = Dtind4\(pX8.^2);
vv8tind4 = pX8.^2 - Dtind4*bb8tind4;
s28tind4 = mean(vv8tind4.^2);
x28tind4 = diag(bb8tind4(2:end,:)'*inv(iDtind4(2:end,2:end))*bb8tind4(2:end,:))./s28tind4';
ptind4(:,2) = 1-chi2cdf(x28tind4,size(Dtind4,2));

disp('tind2');
Dtind2 = sparse(dummyvar(tind2));
Dtind2 = [ones(n,1) Dtind2(:,2:end)];
iDtind2 = inv(Dtind2'*Dtind2);

bb1tind2 = Dtind2\(pX1.^2);
vv1tind2 = pX1.^2 - Dtind2*bb1tind2;
s21tind2 = mean(vv1tind2.^2);
x21tind2 = diag(bb1tind2(2:end,:)'*inv(iDtind2(2:end,2:end))*bb1tind2(2:end,:))./s21tind2'; %#ok<*MINV>
ptind2(:,1) = 1-chi2cdf(x21tind2,size(Dtind2,2)); %#ok<*MNEFF,*NOPTS>

bb9tind2 = Dtind2\(pX9.^2);
vv9tind2 = pX9.^2 - Dtind2*bb9tind2;
s29tind2 = mean(vv9tind2.^2);
x29tind2 = diag(bb9tind2(2:end,:)'*inv(iDtind2(2:end,2:end))*bb9tind2(2:end,:))./s29tind2';
ptind2(:,2) = 1-chi2cdf(x29tind2,size(Dtind2,2));

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

% Point estimate - FM
bFMstate = zeros(nSim,kx); % FM - state
bFMSIC1 = zeros(nSim,kx);  % FM - SIC1
bFMsize = zeros(nSim,kx);  % FM - size
bFMT8 = zeros(nSim,kx); % FM - 8 year time block
bFMT6 = zeros(nSim,kx); % FM - 6 year time block
bFMT4 = zeros(nSim,kx); % FM - 4 year time block
bFMT2 = zeros(nSim,kx); % FM - 2 year time block

badapt = zeros(nSim,kx);  % Adaptive, based on residual regressions

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

% standard errors - two-way clustered
sestatetime = zeros(nSim,kx);  % baseline, TW by state and time
nstatetime = 0; % number not positive definite, state and time
sesic1time = zeros(nSim,kx);   % baseline, TW by sic1 and time
nsic1time = 0; % number not positive definite, sic1 and time
sesic2time = zeros(nSim,kx);   % baseline, TW by sic2 and time
nsic2time = 0; % number not positive definite, sic2 and time

seSxTstatetime = zeros(nSim,kx);  % state X time effects, TW by state and time
nSxTstatetime = 0;  % number not positive definite, state X time effects, TW by state and time

seSIC1xTsic1time = zeros(nSim,kx);  % sic1 X time effects, TW by sic1 and time
nSIC1xTsic1time = 0;  % number not positive definite, sic1 X time effects, TW by sic1 and time

seSIC2xTsic1time = zeros(nSim,kx);  % sic2 X time effects, TW by sic1 and time
nSIC2xTsic1time = 0;  % number not positive definite, sic2 X time effects, TW by sic1 and time
seSIC2xTsic2time = zeros(nSim,kx);  % sic2 X time effects, TW by sic2 and time
nSIC2xTsic2time = 0;  % number not positive definite, sic2 X time effects, TW by sic2 and time

% standard errors - FM
seFMstate = zeros(nSim,kx); % FM - state
seFMSIC1 = zeros(nSim,kx);  % FM - SIC1
seFMsize = zeros(nSim,kx);  % FM - size
seFMT8 = zeros(nSim,kx); % FM - 8 year time block
seFMT6 = zeros(nSim,kx); % FM - 6 year time block
seFMT4 = zeros(nSim,kx); % FM - 4 year time block
seFMT2 = zeros(nSim,kx); % FM - 2 year time block

seadapt = zeros(nSim,kx);  % Adaptive, based on residual regressions
dfadapt = zeros(nSim,1);

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
    
    [~,VS] = cluster_se(X,e,M,state,k);
    [~,VT] = cluster_se(X,e,M,t_ind,k);
    [~,VST] = cluster_se(X,e,M,SxT,k);
    V = VS + VT - VST;
    if min(eig(V)) > 0
        sestatetime(ii,:) = sqrt(diag(V))';
    else
        [temp1,temp2] = eig(V);
        Dstar = max(temp2,0);
        sestatetime(ii,:) = sqrt(diag(temp1*Dstar*temp1'))';
        nstatetime = nstatetime + 1;
    end
    
    [~,VS] = cluster_se(X,e,M,sic1,k);
    [~,VT] = cluster_se(X,e,M,t_ind,k);
    [~,VST] = cluster_se(X,e,M,SIC1xT,k);
    V = VS + VT - VST;
    if min(eig(V)) > 0
        sesic1time(ii,:) = sqrt(diag(V))';
    else
        [temp1,temp2] = eig(V);
        Dstar = max(temp2,0);
        sesic1time(ii,:) = sqrt(diag(temp1*Dstar*temp1'))';
        nsic1time = nsic1time + 1;
    end

    [~,VS] = cluster_se(X,e,M,sic2,k);
    [~,VT] = cluster_se(X,e,M,t_ind,k);
    [~,VST] = cluster_se(X,e,M,SIC2xT,k);
    V = VS + VT - VST;
    if min(eig(V)) > 0
        sesic2time(ii,:) = sqrt(diag(V))';
    else
        [temp1,temp2] = eig(V);
        Dstar = max(temp2,0);
        sesic2time(ii,:) = sqrt(diag(temp1*Dstar*temp1'))';
        nsic2time = nsic2time + 1;
    end
        
    % Firm and State X Time effects
    YFMSxT = y - DFMSxT*(DFMSxT\y);
    bSxTii = XFMSxT\YFMSxT;
    bSxT(ii,:) = bSxTii';
    eSxT = YFMSxT - XFMSxT*bSxTii;

    seSxTstate(ii,:) = cluster_se(XFMSxT,eSxT,MSxT,state,kSxT)';

    [~,VS] = cluster_se(XFMSxT,eSxT,MSxT,state,kSxT);
    [~,VT] = cluster_se(XFMSxT,eSxT,MSxT,t_ind,kSxT);
    [~,VST] = cluster_se(XFMSxT,eSxT,MSxT,SxT,kSxT);
    V = VS + VT - VST;
    if min(eig(V)) > 0
        seSxTstatetime(ii,:) = sqrt(diag(V))';
    else
        [temp1,temp2] = eig(V);
        Dstar = max(temp2,0);
        seSxTstatetime(ii,:) = sqrt(diag(temp1*Dstar*temp1'))';
        nSxTstatetime = nSxTstatetime + 1;
    end
    
    % Firm and SIC1 X Time effects
    YFMSIC1xT = y - DFMSIC1xT*(DFMSIC1xT\y);
    bSIC1xTii = XFMSIC1xT\YFMSIC1xT;
    bSIC1xT(ii,:) = bSIC1xTii';
    eSIC1xT = YFMSIC1xT - XFMSIC1xT*bSIC1xTii;
    
    seSIC1xTsic1(ii,:) = cluster_se(XFMSIC1xT,eSIC1xT,MSIC1xT,sic1,kSIC1xT)';

    [~,VS] = cluster_se(XFMSIC1xT,eSIC1xT,MSIC1xT,sic1,kSIC1xT);
    [~,VT] = cluster_se(XFMSIC1xT,eSIC1xT,MSIC1xT,t_ind,kSIC1xT);
    [~,VST] = cluster_se(XFMSIC1xT,eSIC1xT,MSIC1xT,SIC1xT,kSIC1xT);
    V = VS + VT - VST;
    if min(eig(V)) > 0
        seSIC1xTsic1time(ii,:) = sqrt(diag(V))';
    else
        [temp1,temp2] = eig(V);
        Dstar = max(temp2,0);
        seSIC1xTsic1time(ii,:) = sqrt(diag(temp1*Dstar*temp1'))';
        nSIC1xTsic1time = nSIC1xTsic1time + 1;
    end
    
    % Firm and SIC2 X Time effects
    YFMSIC2xT = y - DFMSIC2xT*(DFMSIC2xT\y);
    bSIC2xTii = XFMSIC2xT\YFMSIC2xT;
    bSIC2xT(ii,:) = bSIC2xTii';
    eSIC2xT = YFMSIC2xT - XFMSIC2xT*bSIC2xTii;

    seSIC2xTsic1(ii,:) = cluster_se(XFMSIC2xT,eSIC2xT,MSIC2xT,sic1,kSIC2xT)';
    seSIC2xTsic2(ii,:) = cluster_se(XFMSIC2xT,eSIC2xT,MSIC2xT,sic2,kSIC2xT)';
    
    [~,VS] = cluster_se(XFMSIC2xT,eSIC2xT,MSIC2xT,sic1,kSIC2xT);
    [~,VT] = cluster_se(XFMSIC2xT,eSIC2xT,MSIC2xT,t_ind,kSIC2xT);
    [~,VST] = cluster_se(XFMSIC2xT,eSIC2xT,MSIC2xT,SIC1xT,kSIC2xT);
    V = VS + VT - VST;
    if min(eig(V)) > 0
        seSIC2xTsic1time(ii,:) = sqrt(diag(V))';
    else
        [temp1,temp2] = eig(V);
        Dstar = max(temp2,0);
        seSIC2xTsic1time(ii,:) = sqrt(diag(temp1*Dstar*temp1'))';
        nSIC2xTsic1time = nSIC2xTsic1time + 1;
    end
    
    [~,VS] = cluster_se(XFMSIC2xT,eSIC2xT,MSIC2xT,sic2,kSIC2xT);
    [~,VT] = cluster_se(XFMSIC2xT,eSIC2xT,MSIC2xT,t_ind,kSIC2xT);
    [~,VST] = cluster_se(XFMSIC2xT,eSIC2xT,MSIC2xT,SIC2xT,kSIC2xT);
    V = VS + VT - VST;
    if min(eig(V)) > 0
        seSIC2xTsic2time(ii,:) = sqrt(diag(V))';
    else
        [temp1,temp2] = eig(V);
        Dstar = max(temp2,0);
        seSIC2xTsic2time(ii,:) = sqrt(diag(temp1*Dstar*temp1'))';
        nSIC2xTsic2time = nSIC2xTsic2time + 1;
    end
    
    % Firm and Size X Time effects
    YFMSIZExT = y - DFMSIZExT*(DFMSIZExT\y);
    bSIZExTii = XFMSIZExT\YFMSIZExT;
    bSIZExT(ii,:) = bSIZExTii';
    eSIZExT = YFMSIZExT - XFMSIZExT*bSIZExTii;

    seSIZExTsize(ii,:) = cluster_se(XFMSIZExT,eSIZExT,MSIZExT,sizecat,kSIZExT)';

    % Firm X 8 year block and Time effects
    YFMFxT8 = y - DFMFxT8*(DFMFxT8\y);
    bFxT8ii = XFMFxT8\YFMFxT8;
    bFxT8(ii,:) = bFxT8ii;
    eFxT8 = YFMFxT8 - XFMFxT8*bFxT8ii;
    
    seFxT8(ii,:) = cluster_se(XFMFxT8,eFxT8,MFxT8,tind8,kFxT8)';
    
    % Firm X 6 year block and Time effects
    YFMFxT6 = y - DFMFxT6*(DFMFxT6\y);
    bFxT6ii = XFMFxT6\YFMFxT6;
    bFxT6(ii,:) = bFxT6ii;
    eFxT6 = YFMFxT6 - XFMFxT6*bFxT6ii;
    
    seFxT6(ii,:) = cluster_se(XFMFxT6,eFxT6,MFxT6,tind6,kFxT6)';

    % Firm X 4 year block and Time effects
    YFMFxT4 = y - DFMFxT4*(DFMFxT4\y);
    bFxT4ii = XFMFxT4\YFMFxT4;
    bFxT4(ii,:) = bFxT4ii;
    eFxT4 = YFMFxT4 - XFMFxT4*bFxT4ii;
    
    seFxT4(ii,:) = cluster_se(XFMFxT4,eFxT4,MFxT4,tind4,kFxT4)';

    % Firm X 2 year block and Time effects
    YFMFxT2 = y - DFMFxT2*(DFMFxT2\y);
    bFxT2ii = XFMFxT2\YFMFxT2;
    bFxT2(ii,:) = bFxT2ii;
    eFxT2 = YFMFxT2 - XFMFxT2*bFxT2ii;

    seFxT2(ii,:) = cluster_se(XFMFxT2,eFxT2,MFxT2,tind2,kFxT2)';

    % FAMA-MACBETH

    % State groups
    [tmp1,tmp2] = FamaMacbeth(XFMSxT,YFMSxT,state);
    bFMstate(ii,:) = tmp1(1:kx);
    seFMstate(ii,:) = tmp2(1:kx);
    
    % SIC1 groups
    [tmp1,tmp2] = FamaMacbeth(XFMSIC1xT,YFMSIC1xT,sic1);
    bFMSIC1(ii,:) = tmp1(1:kx);
    seFMSIC1(ii,:) = tmp2(1:kx);
    
    % Size category groups
    [tmp1,tmp2] = FamaMacbeth(XFMSIZExT,YFMSIZExT,sizecat);
    bFMsize(ii,:) = tmp1(1:kx);
    seFMsize(ii,:) = tmp2(1:kx);

    % T8 groups
    [tmp1,tmp2] = FamaMacbeth(XFMFxT8,YFMFxT8,tind8);
    bFMT8(ii,:) = tmp1(1:kx);
    seFMT8(ii,:) = tmp2(1:kx);

    % T6 groups
    [tmp1,tmp2] = FamaMacbeth(XFMFxT6,YFMFxT6,tind6);
    bFMT6(ii,:) = tmp1(1:kx);
    seFMT6(ii,:) = tmp2(1:kx);
    
    % T4 groups
    [tmp1,tmp2] = FamaMacbeth(XFMFxT4,YFMFxT4,tind4);
    bFMT4(ii,:) = tmp1(1:kx);
    seFMT4(ii,:) = tmp2(1:kx);
    
    % T2 groups
    [tmp1,tmp2] = FamaMacbeth(XFMFxT2,YFMFxT2,tind2);
    bFMT2(ii,:) = tmp1(1:kx);
    seFMT2(ii,:) = tmp2(1:kx);
    
    % Adaptive
    if Trep(ii) == 2
        badapt(ii,:) = bFMT2(ii,:);
        seadapt(ii,:) = seFMT2(ii,:);
        dfadapt(ii) = 8;
    elseif Trep(ii) == 4
        badapt(ii,:) = bFMT4(ii,:);
        seadapt(ii,:) = seFMT4(ii,:);
        dfadapt(ii) = 4;
    elseif Trep(ii) == 6
        badapt(ii,:) = bFMT6(ii,:);
        seadapt(ii,:) = seFMT6(ii,:);
        dfadapt(ii) = 3;
    end
        
end

save FRQSim_FINAL_PART1 ;
