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

%%%% Time blocks
% 2 year block
tind2 = ceil(t_ind/2);
tind2(tind2 == 9) = 8;

% 4 year block
tind4 = ceil(t_ind/4);
tind4(tind4 == 5) = 4;

% 6 year block
tind6 = ceil(t_ind/6);

% 8 year block
tind8 = ceil(t_ind/8);
tind8(tind8 == 3) = 2;

%%% Pre-set variables
Firms=unique(f_ind);    % Unique firm identifiers                    
J = size(Firms,1);      % Number of firms
T= max(t_ind);          % Number of periods
dataref=(1:T)';         % Complete time periods
kx = size(x,2);         % Number of variables
nSim    = 1000 ;         % Number of simulation replications
bsize   = 3;            % Bloc size
BB      = 2000;           % Number of repetitions MBB
rng(8)                  % Fix the seed

%%% A) Construct once for all the full matrices to be used in bloc bootsrap
xFull  =NaN*ones(J*T,kx);

% firm and time indexes
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

n = size(X,1);
kx = size(X,2);
k = kx+trace(D\D);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3) Predefine variables for simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A) Point estimates
b = zeros(nSim,kx);                % Baseline estimator (Firm and Time fixed effects) 

% Clustered standard errors
se_t2  = zeros(size(b));          % 2 year time block
se_t4  = zeros(size(b));         % 4 year time block
se_t6  = zeros(size(b));         % 6 year time block
se_t8  = zeros(size(b));         % 8 year time block


% B) Boostrap critical values
bcv_t2 = zeros(size(b));         % 2 year time block
bcv_t4  = zeros(size(b));        % 4 year time block
bcv_t6  = zeros(size(b));        % 6 year time block
bcv_t8  = zeros(size(b));        % 8 year time block

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
% time block 2
[se,~,~] = cluster_se(X,e,M,tind2,kx);
se_t2(ii,:) = (sqrt((n-kx)/(n-k))*se)';
% time block 4
[se,~,~] = cluster_se(X,e,M,tind4,kx);
se_t4(ii,:) = (sqrt((n-kx)/(n-k))*se)';
% time block 6
[se,~,~] = cluster_se(X,e,M,tind6,kx);
se_t6(ii,:) = (sqrt((n-kx)/(n-k))*se)';
% time block 8
[se,~,~] = cluster_se(X,e,M,tind8,kx);
se_t8(ii,:) = (sqrt((n-kx)/(n-k))*se)';

%%% D) Moving Bloc Bootsrap 
bboot = zeros(BB,kx);           % Point estimates

bbT_t2    = zeros(size(bboot));     % T-statistic 2 year block
bbT_t4    = zeros(size(bboot));    % T-statistic 4 year block
bbT_t6    = zeros(size(bboot));    % T-statistic 6 year block
bbT_t8    = zeros(size(bboot));    % T-statistic 8 year block

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
    fb    = f_full(Index,:);         % Firm Index
    tb    = t_full(Index,:);         % Time Index
    nb    = size(xb,1);              % Number of observations
    
    %%%%%% Re-construct dummy matrices
    Dbi = sparse(dummyvar(fb));
    Dbt = sparse(dummyvar(tb));
    DB  = [Dbi Dbt];
            
    %%%% Time blocks for clusters
    % 2 year block
    tind2_ = ceil(tb/2);
    tind2_(tind2_ == 9) = 8;
    
    % 4 year block
    tind4_ = ceil(tb/4);
    tind4_(tind4_ == 5) = 4;
    
    % 6 year block
    tind6_ = ceil(tb/6);
    
    % 8 year block
    tind8_ = ceil(tb/8);
    tind8_(tind8_ == 3) = 2;
    
    % Partial out time and firm fixed effects    
    Xb = xb - DB*(DB\xb);
    Yb = yb - DB*(DB\yb);
    % Bootstrap Point Estimates
    bb= Xb\Yb;
    Mb=inv(Xb'*Xb);
    eb = Yb-Xb*bb;
    bboot(indb,:) = bb';
    %%%% Clustered Standard errors
    % 2 year block
    [sebb,~,~] = cluster_se(Xb,eb,Mb,tind2_,kx);
    sebb = (sqrt((nb-kx)/(nb-k))*sebb )';
    bbT_t2(indb,:) = (bb-bii)'./sebb; 
    % 4 year block
    [sebb,~,~] = cluster_se(Xb,eb,Mb,tind4_,kx);
    sebb = (sqrt((nb-kx)/(nb-k))*sebb )';
    bbT_t4(indb,:) = (bb-bii)'./sebb; 
    % 6 year block
    [sebb,~,~] = cluster_se(Xb,eb,Mb,tind6_,kx);
    sebb = (sqrt((nb-kx)/(nb-k))*sebb )';
    bbT_t6(indb,:) = (bb-bii)'./sebb; 
    % 8 year block
    [sebb,~,~] = cluster_se(Xb,eb,Mb,tind8_,kx);
    sebb = (sqrt((nb-kx)/(nb-k))*sebb )';
    bbT_t8(indb,:) = (bb-bii)'./sebb;
    end
    

bcv_t2(ii,:) =   prctile(abs(bbT_t2),95);         % 2 year block
bcv_t4(ii,:)  =  prctile(abs(bbT_t4),95);         % 4 year block
bcv_t6(ii,:)  =  prctile(abs(bbT_t6),95);         % 6 year block
bcv_t8(ii,:)  =  prctile(abs(bbT_t8),95);         % 8 year block

fMat=strcat('Results/Pt2/SizePt2_',num2str(nSim),'X',num2str(BB),'_','bsize',num2str(bsize),'_',num2str(minObs),'Obs','_sim_',num2str(ii),'.mat');
save(fMat,'b','se_t2','se_t4','se_t6','se_t8','bcv_t2','bcv_t4','bcv_t6','bcv_t8','bTest3')

end

toc
