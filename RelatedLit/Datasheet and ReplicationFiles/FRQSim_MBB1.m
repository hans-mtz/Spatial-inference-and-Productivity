% Moving Block Bootstrap
% June 2017, Aldo Sandoval 
clc
clear
warning('off','MATLAB:rankDeficientMatrix');
yDir = 'SimulatedOutcomes' ; 
tic
% Minimum number of observations per firm
minObs=17;

% Ranges of simulations to Run
rn1=1;
rn2=1000;

%%% Load data
load XCatIndices ;      % Indices of categories used to form lagged variables
                        % [f_ind state sic2 sizecat cashcat rescat frqcat intcat]
load TimeIndices ;      % Time indices
load FirmIndices ;      % Firm indices
load xMat ;             % Design matrix from actual data
                        % [res ind frq frqint cash tq mve lev age]

load('olsparams.mat')   % True Parameters

%%% Only keep firms with minimum # of observations
F=unique(f_ind);
fms= grpstats(t_ind,f_ind,'numel');
fsel=F((fms>=minObs),1);
ind=ismember(f_ind,fsel);
I_=ind;

f_ind  = f_ind(ind,1);
t_ind  = t_ind(ind,1);
x      = x(ind,:);
indX   = indX(ind,:);
f_ind  =recode(f_ind);        % re index firms so there are no gaps

% Variables for clusters (subsample)
state = indX(:,2);
sic2  = indX(:,3);
sizec = indX(:,4);
sic1  = recode(floor(sic2/10)+1);

%%% Pre-set variables
Firms=unique(f_ind);    % Unique firm identifiers                    
J = size(Firms,1);      % Number of firms
T= max(t_ind);          % Number of periods
dataref=(1:T)';         % Complete time periods
kx = size(x,2);         % Number of variables
nSim    = 1000 ;        % Number of simulation replications
bsize   = 3;            % Bloc size
BB      = 2000;            % Number of repetitions MBB
n = size(x,1);          % Number of osbervations
rng(8)                  % Fix the seed

%%% A) Construct once for all the full matrices to be used in bloc bootsrap
xFull  =NaN*ones(J*T,kx);

r=repmat((1:J)',1,T)';
f_full=r(:);
r=repmat((1:T)',J,1);
t_full=r;

IS=[];
for i=1:J
    ind=(f_ind==i);                        % Get observations from fiirm with id==ii
    tav = t_ind(ind);                       % Time-periods available
    IS= [IS; ismember(dataref,tav)];
end
IS=(IS>0);
xFull(IS,1:kx) =x(:,1:kx);
stFull(IS,1)=state;
sic1F(IS,1) =sic1;
sic2F(IS,1) =sic2;
sizeF(IS,1) =sizec;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) Some data preprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% B) Partial time and firm effects out from design matrix for baseline results
% Firm,year 
Di = sparse(dummyvar(f_ind));
Dt = sparse(dummyvar(t_ind));
D = [Di Dt];
X = x - D*(D\x);
M = inv(X'*X);

kx = size(X,2);
k = kx+trace(D\D);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3) Predefine variables for simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A) Point estimates
b = zeros(nSim,kx);                % Baseline estimator (Firm and Time fixed effects) 

% Clustered standard errors
se_firm  = zeros(size(b));         % firm
se_state = zeros(size(b));         % state
se_sic1  = zeros(size(b));         % one-digit SIC
se_sic2  = zeros(size(b));         % two-digit SIC
se_size  = zeros(size(b));         % size


% B) Boostrap critical values
bcv_firm  = zeros(size(b));         % firm
bcv_state = zeros(size(b));         % state
bcv_sic1  = zeros(size(b));         % one-digit SIC
bcv_sic2  = zeros(size(b));         % two-digit SIC
bcv_size  = zeros(size(b));         % size

for ii = rn1:rn2

disp(ii);

% Load presimulated outcome data
filename = strcat(yDir,'/ySim',num2str(ii),'.mat');
load(filename);
y = y(I_,:);
yFull=NaN*ones(J*T,1);
yFull(IS,:)=y;

% C) Baseline estimate
Y = y - D*(D\y);
bii = X\Y  ;
b(ii,:) = bii';
e = Y-X*bii;

%%% Clustered standard errors
% Firm
[se,~,~] = cluster_se(X,e,M,f_ind,kx);
se_firm(ii,:) = (sqrt((n-kx)/(n-k))*se)';
% State
[se,~,~] = cluster_se(X,e,M,state,kx);
se_state(ii,:) = (sqrt((n-kx)/(n-k))*se)';
% SIC 1
[se,~,~] = cluster_se(X,e,M,sic1,kx);
se_sic1(ii,:) = (sqrt((n-kx)/(n-k))*se)';
% SIC 2
[se,~,~] = cluster_se(X,e,M,sic2,kx);
se_sic2(ii,:) = (sqrt((n-kx)/(n-k))*se)';
% Size
[se,~,~] = cluster_se(X,e,M,sizec,kx);
se_size(ii,:) = (sqrt((n-kx)/(n-k))*se)';

%%% D) Moving Bloc Bootsrap 
bboot = zeros(BB,kx);           % Point estimates

bbT_firm    = zeros(size(bboot));    % T-statistic cluster firm
bbT_state   = zeros(size(bboot));    % T-statistic cluster state
bbT_sic1    = zeros(size(bboot));    % T-statistic cluster sic1
bbT_sic2    = zeros(size(bboot));    % T-statistic cluster sic2
bbT_size    = zeros(size(bboot));    % T-statistic cluster size

I = randi([1 T-bsize+1],[BB,ceil(T/bsize)]) ;   % sample indexes
Isel=zeros(T*J,1);

    parfor indb=1:BB
    warning('off','MATLAB:rankDeficientMatrix');
    %%%%%% Moving Block Bootsrap
    Iplus = index_mbb( I(indb,:), bsize, T, J );
    Isel=(f_full -1)*T + Iplus;
    xb=xFull(Isel,:);               % resampled design matrix
    yb=yFull(Isel,:);               % resampled outcome
    
    % Get rid of missing observations
    inan  = (isnan(xb)==1);
    Index = (sum(inan,2)==0);
    xb    = xb(Index,:);             
    yb    = yb(Index,:);    
    fb    = f_full(Index,:);        % Firm Index
    tb    = t_full(Index,:);        % Time Index
    Sb    = stFull(Index,:);        % State Index
    sb    = sizeF(Index,:);         % Size Index  
    s1b   = sic1F(Index,:);         % SIC1 Index  
    s2b   = sic2F(Index,:);         % SIC2 index
    nb    = size(xb,1);             % Number of observations

    %%%%%% Re-construct dummy matrices
    Dbi = sparse(dummyvar(fb));
    Dbt = sparse(dummyvar(tb));
    DB  = [Dbi Dbt];
    
    % Partial out time and firm fixed effects    
    Xb = xb - DB*(DB\xb);
    Yb = yb - DB*(DB\yb);
    % Bootstrap Point Estimates
    bb= Xb\Yb;
    Mb=inv(Xb'*Xb);
    eb = Yb-Xb*bb;
    bboot(indb,:) = bb';
    %%%% Clustered Standard errors
    % Firm
    [sebb,~,~] = cluster_se(Xb,eb,Mb,fb,kx);
    sebb = (sqrt((nb-kx)/(nb-k))*sebb )';
    bbT_firm(indb,:) = (bb-bii)'./sebb; 
    % State
    [sebb,~,~] = cluster_se(Xb,eb,Mb,Sb,kx);
    sebb = (sqrt((nb-kx)/(nb-k))*sebb )';
    bbT_state(indb,:) = (bb-bii)'./sebb; 
    % SIC1
    [sebb,~,~] = cluster_se(Xb,eb,Mb,s1b,kx);
    sebb = (sqrt((nb-kx)/(nb-k))*sebb )';
    bbT_sic1(indb,:) = (bb-bii)'./sebb; 
    % SIC2
    [sebb,~,~] = cluster_se(Xb,eb,Mb,s2b,kx);
    sebb = (sqrt((nb-kx)/(nb-k))*sebb )';
    bbT_sic2(indb,:) = (bb-bii)'./sebb; 
    % Size
    [sebb,~,~] = cluster_se(Xb,eb,Mb,sb,kx);
    sebb = (sqrt((nb-kx)/(nb-k))*sebb )';
    bbT_size(indb,:) = (bb-bii)'./sebb;
    end
    
bcv_firm(ii,:)  =  prctile(abs(bbT_firm),95);         % firm
bcv_state(ii,:) =  prctile(abs(bbT_state),95);        % state
bcv_sic1(ii,:)  =  prctile(abs(bbT_sic1),95);         % one-digit SIC
bcv_sic2(ii,:)  =  prctile(abs(bbT_sic2),95);         % two-digit SIC
bcv_size(ii,:)  =  prctile(abs(bbT_size),95);         % Size
    
 fMat=strcat('Results/Pt1/SizePt1_',num2str(nSim),'X',num2str(BB),'_','bsize',num2str(bsize),'_',num2str(minObs),'Obs','_sim_',num2str(ii),'.mat');
 save(fMat,'b','se_firm','se_state','se_sic1','se_sic2','se_size','bcv_firm','bcv_state','bcv_sic1','bcv_sic2','bcv_size','bTest3')

end

toc