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

%%% Only keep firms with the minimum # of observations
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
T2 = groupcross(f_ind,tind2);

% 4 year block
tind4 = ceil(t_ind/4);
tind4(tind4 == 5) = 4;
T4 = groupcross(f_ind,tind4);

% 6 year block
tind6 = ceil(t_ind/6);
T6 = groupcross(f_ind,tind6);

% 8 year block
tind8 = ceil(t_ind/8);
tind8(tind8 == 3) = 2;
T8 = groupcross(f_ind,tind8);


%%% Pre-set variables
Firms=unique(f_ind);    % Unique firm identifiers                    
J = size(Firms,1);      % Number of firms
T= max(t_ind);          % Number of periods
dataref=(1:T)';         % Complete time periods
kx = size(x,2);         % Number of variables
nSim    = 1000 ;        % Number of simulation replications
bsize   = 3;            % Bloc size
BB      = 2000;         % Number of repetitions MBB
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) Some data preprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% B) Partial time and firm effects out from design matrix for baseline results
Dt    = sparse(dummyvar(t_ind));
% Firm 2 year block, year 
DFt2    = sparse(dummyvar(T2));
ind = (sum(DFt2)==0);               % Get rid of columns of zeros
DFt2(:,ind)=[];

DFt2xT = [Dt DFt2];
XFt2xT = x - DFt2xT*(DFt2xT\x);
kFt2xT = kx+trace(DFt2xT\DFt2xT);
MFt2xT = inv(XFt2xT'*XFt2xT);

% Firm 4 year block, year
DFt4    = sparse(dummyvar(T4));
ind = (sum(DFt4)==0);
DFt4(:,ind)=[];

DFt4xT = [Dt DFt4];
XFt4xT = x - DFt4xT*(DFt4xT\x);
kFt4xT = kx+trace(DFt4xT\DFt4xT);
MFt4xT = inv(XFt4xT'*XFt4xT);

% Firm 6 year block, year 
DFt6    = sparse(dummyvar(T6));
ind = (sum(DFt6)==0);
DFt6(:,ind)=[];

DFt6xT = [Dt DFt6];
XFt6xT = x - DFt6xT*(DFt6xT\x);
kFt6xT = kx+trace(DFt6xT\DFt6xT);
MFt6xT = inv(XFt6xT'*XFt6xT);

% Firm 8 year block, year 
DFt8    = sparse(dummyvar(T8));
ind = (sum(DFt8)==0);
DFt8(:,ind)=[];

DFt8xT = [Dt DFt8];
XFt8xT = x - DFt8xT*(DFt8xT\x);
kFt8xT = kx+trace(DFt8xT\DFt8xT);
MFt8xT = inv(XFt8xT'*XFt8xT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3) Predefine variables for simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A) Point estimates
bFt2xT = zeros(nSim,kx);
bFt4xT = zeros(nSim,kx);
bFt6xT = zeros(nSim,kx);
bFt8xT = zeros(nSim,kx);

% Clustered standard errors
se_t2xT  = zeros(size(bFt2xT));         
se_t4xT  = zeros(size(bFt4xT));         
se_t6xT  = zeros(size(bFt6xT));         
se_t8xT  = zeros(size(bFt8xT));         

% B) Boostrap critical values
bcv_t2xT  = zeros(size(bFt2xT)); 
bcv_t4xT  = zeros(size(bFt4xT));
bcv_t6xT  = zeros(size(bFt6xT));
bcv_t8xT  = zeros(size(bFt8xT));


for ii = rn1:rn2
 
disp(ii);

% Load presimulated outcome data
filename = strcat(yDir,'/ySim',num2str(ii),'.mat');
load(filename);
y = y(I_,:);                % Only select firms in subsample
yFull=NaN*ones(J*T,1);
yFull(IS,:)=y;

%%%% C) Point estimates
% Firm X 2 year block, year
YFt2xT = y - DFt2xT*(DFt2xT\y);
biit2 = XFt2xT\YFt2xT  ;
bFt2xT(ii,:) = biit2';
e = YFt2xT-XFt2xT*biit2;
% standard errors
[se,~,~] = cluster_se(XFt2xT,e,MFt2xT,tind2,kx);
se_t2xT(ii,:) = (sqrt((n-kx)/(n-kFt2xT))*se)';

% Firm X 4 year block, year
YFt4xT = y - DFt4xT*(DFt4xT\y);
biit4 = XFt4xT\YFt4xT  ;
bFt4xT(ii,:) = biit4';
e = YFt4xT-XFt4xT*biit4;
% standard errors
[se,~,~] = cluster_se(XFt4xT,e,MFt4xT,tind4,kx);
se_t4xT(ii,:) = (sqrt((n-kx)/(n-kFt4xT))*se)';

% Firm X 6 year block, year
YFt6xT = y - DFt6xT*(DFt6xT\y);
biit6 = XFt6xT\YFt6xT  ;
bFt6xT(ii,:) = biit6';
e = YFt6xT-XFt6xT*biit6;
% standard errors
[se,~,~] = cluster_se(XFt6xT,e,MFt6xT,tind6,kx);
se_t6xT(ii,:) = (sqrt((n-kx)/(n-kFt6xT))*se)';

% Firm X 8 year block, year
YFt8xT = y - DFt8xT*(DFt8xT\y);
biit8 = XFt8xT\YFt8xT  ;
bFt8xT(ii,:) = biit8';
e = YFt8xT-XFt8xT*biit8;
% standard errors
[se,~,~] = cluster_se(XFt8xT,e,MFt8xT,tind8,kx);
se_t8xT(ii,:) = (sqrt((n-kx)/(n-kFt8xT))*se)';

%%%% D) Moving Bloc Bootsrap 
bb_t2 = zeros(BB,kx);          % Point estimates
bb_t4  = zeros(BB,kx);           % Point estimates
bb_t6  = zeros(BB,kx);           % Point estimates
bb_t8  = zeros(BB,kx);           % Point estimates

bbT_t2    = zeros(size(bb_t2));    % T-statistic cluster firm
bbT_t4    = zeros(size(bb_t4));    % T-statistic cluster firm
bbT_t6    = zeros(size(bb_t6));    % T-statistic cluster firm
bbT_t8    = zeros(size(bb_t8));    % T-statistic cluster firm

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
    Dbt = sparse(dummyvar(tb));
    %%%% Time blocks
    % 2 year block
    tind2_ = ceil(tb/2);
    tind2_(tind2_ == 9) = 8;
    tb2 = groupcross(fb,tind2_);
    df    = sparse(dummyvar(tb2));
    ind = (sum(df)==0);
    df(:,ind)=[];
    Dt2 = [df Dbt];
    
    % 4 year block
    tind4_= ceil(tb/4);
    tind4_(tind4_ == 5) = 4;
    tb4 = groupcross(fb,tind4_);
    df    = sparse(dummyvar(tb4));
    ind = (sum(df)==0);
    df(:,ind)=[];
    Dt4 = [df Dbt];
    
    % 6 year block
    tind6_ = ceil(tb/6);
    tb6 = groupcross(fb,tind6_);
    df    = sparse(dummyvar(tb6));
    ind = (sum(df)==0);
    df(:,ind)=[];
    Dt6 = [df Dbt];
    
    % 8 year block
    tind8_ = ceil(tb/8);
    tind8_(tind8_ == 3) = 2;
    tb8 = groupcross(fb,tind8_);
    df    = sparse(dummyvar(tb8));
    ind = (sum(df)==0);
    df(:,ind)=[];
    Dt8 = [df Dbt];
    
    % Firm X 2 year block, year    
    Xbt2 = xb - Dt2*(Dt2\xb);
    Ybt2 = yb - Dt2*(Dt2\yb);    

    bb   = Xbt2\Ybt2;
    Mbt2 = inv(Xbt2'*Xbt2);
    ebt2 = Ybt2-Xbt2*bb;
    bb_t2(indb,:) = bb';
    
    [sebb,~,~] = cluster_se(Xbt2,ebt2,Mbt2,tind2_,kx);
    kft2 = kx+trace(Dt2\Dt2);
    sebb = (sqrt((nb-kx)/(nb-kft2))*sebb )';
    bbT_t2(indb,:) = (bb-biit2)'./sebb;
    
    % Firm X 4 year block, year    
    Xbt4 = xb - Dt4*(Dt4\xb);
    Ybt4 = yb - Dt4*(Dt4\yb);    

    bb   = Xbt4\Ybt4;
    Mbt4 = inv(Xbt4'*Xbt4);
    ebt4 = Ybt4-Xbt4*bb;
    bb_t4(indb,:) = bb';
    
    [sebb,~,~] = cluster_se(Xbt4,ebt4,Mbt4,tind4_,kx);
    kft4 = kx+trace(Dt4\Dt4);
    sebb = (sqrt((nb-kx)/(nb-kft4))*sebb )';
    bbT_t4(indb,:) = (bb-biit4)'./sebb;
    
    % Firm X 6 year block, year    
    Xbt6 = xb - Dt6*(Dt6\xb);
    Ybt6 = yb - Dt6*(Dt6\yb);    

    bb   = Xbt6\Ybt6;
    Mbt6 = inv(Xbt6'*Xbt6);
    ebt6 = Ybt6-Xbt6*bb;
    bb_t6(indb,:) = bb';
    
    [sebb,~,~] = cluster_se(Xbt6,ebt6,Mbt6,tind6_,kx);
    kft6 = kx+trace(Dt6\Dt6);
    sebb = (sqrt((nb-kx)/(nb-kft6))*sebb )';
    bbT_t6(indb,:) = (bb-biit6)'./sebb;    
    
    % Firm X 8 year block, year    
    Xbt8 = xb - Dt8*(Dt8\xb);
    Ybt8 = yb - Dt8*(Dt8\yb);    

    bb   = Xbt8\Ybt8;
    Mbt8 = inv(Xbt8'*Xbt8);
    ebt8 = Ybt8-Xbt8*bb;
    bb_t8(indb,:) = bb';
    
    [sebb,~,~] = cluster_se(Xbt8,ebt8,Mbt8,tind8_,kx);
    kft8 = kx+trace(Dt8\Dt8);
    sebb = (sqrt((nb-kx)/(nb-kft8))*sebb )';
    bbT_t8(indb,:) = (bb-biit8)'./sebb;    
    
    end
    
bcv_t2xT(ii,:) =  prctile(abs(bbT_t2),95);    
bcv_t4xT(ii,:) =  prctile(abs(bbT_t4),95); 
bcv_t6xT(ii,:) =  prctile(abs(bbT_t6),95); 
bcv_t8xT(ii,:) =  prctile(abs(bbT_t8),95);

fMat=strcat('Results/Pt4/SizePt4_',num2str(nSim),'X',num2str(BB),'_','bsize',num2str(bsize),'_',num2str(minObs),'Obs','_sim_',num2str(ii),'.mat');
save(fMat,'bFt2xT','bFt4xT','bFt6xT','bFt8xT','se_t2xT','se_t4xT','se_t6xT','se_t8xT','bcv_t2xT','bcv_t4xT','bcv_t6xT','bcv_t8xT','bTest3');


end

% Bsize_t2    = mean(abs(bFt2xT-ones(nSim,1)*bTest3')./se_t2xT > bcv_t2xT)
% Bsize_t4    = mean(abs(bFt4xT-ones(nSim,1)*bTest3')./se_t4xT > bcv_t4xT)
% Bsize_t6    = mean(abs(bFt6xT-ones(nSim,1)*bTest3')./se_t6xT > bcv_t6xT)
% Bsize_t8    = mean(abs(bFt8xT-ones(nSim,1)*bTest3')./se_t8xT > bcv_t8xT)
toc
