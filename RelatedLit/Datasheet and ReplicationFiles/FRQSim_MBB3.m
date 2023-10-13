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

% Variables for fixed effects & clusters
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
nSim    = 1000;         % Number of simulation replications
bsize   = 3;            % Bloc size
BB      = 2000;         % Number of repetitions MBB
n = size(x,1);
rng(8)                  % Fix the seed

%%% A) Construct once for all the full matrices to be used in bloc bootsrap
xFull  =NaN*ones(J*T,kx);

% firm and time indexes
r=repmat((1:J)',1,T)';
f_full=r(:);
r=repmat((1:T)',J,1);
t_full=r;

% Create full arrays for clusters and fixed effects
stFull =NaN*ones(J*T,1);
sic1F  =NaN*ones(J*T,1);
sic2F  =NaN*ones(J*T,1); 
sizeF  =NaN*ones(J*T,1);

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

Di  = sparse(dummyvar(f_ind));

% Firm, state x year
SxT = groupcross(state,t_ind);
DSxT = sparse(dummyvar(SxT));
DFMSxT = [Di DSxT];
XFMSxT = x - DFMSxT*(DFMSxT\x);
kSxT = kx+trace(DFMSxT\DFMSxT);
MSxT = inv(XFMSxT'*XFMSxT);

% Firm, one-digit x year
sc1xT = groupcross(sic1,t_ind);

Dsc1xT = sparse(dummyvar(sc1xT));
DFsc1xT = [Di Dsc1xT];
XFsc1xT = x - DFsc1xT*(DFsc1xT\x);
ksc1xT = kx+trace(DFsc1xT\DFsc1xT);
Msc1xT = inv(XFsc1xT'*XFsc1xT);

% Firm, two-digit x year
sc2xT = groupcross(sic2,t_ind);

Dsc2xT = sparse(dummyvar(sc2xT));
DFsc2xT = [Di Dsc2xT];
XFsc2xT = x - DFsc2xT*(DFsc2xT\x);
ksc2xT = kx+trace(DFsc2xT\DFsc2xT);
Msc2xT = inv(XFsc2xT'*XFsc2xT);

% Firm, size x year
szexT = groupcross(sizec,t_ind);

DszexT = sparse(dummyvar(szexT));
DFszexT = [Di DszexT];
XFszexT = x - DFszexT*(DFszexT\x);
kszexT = kx+trace(DFszexT\DFszexT);
MszexT = inv(XFszexT'*XFszexT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3) Predefine variables for simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A) Point estimates
bFMSxT = zeros(nSim,kx);           % Baseline estimator (Firm and StateXTime fixed effects)
bFsc1xT = zeros(nSim,kx);           % Baseline estimator (Firm and one-digitXTime fixed effects)
bFsc2xT = zeros(nSim,kx);           % Baseline estimator (Firm and two-digitXTime fixed effects)
bFszexT = zeros(nSim,kx);           % Baseline estimator (Firm and SizeXTime fixed effects)

% Clustered standard errors
se_state = zeros(size(bFMSxT));         % state
se_sic1  = zeros(size(bFsc1xT));         % one-digit SIC
se_sic2  = zeros(size(bFsc2xT));         % two-digit SIC
se_size  = zeros(size(bFszexT));         % size

% B) Boostrap critical values
bcv_state = zeros(size(bFMSxT));         % state
bcv_sic1  = zeros(size(bFsc1xT));         % one-digit SIC
bcv_sic2  = zeros(size(bFsc2xT));         % two-digit SIC
bcv_size  = zeros(size(bFszexT));         % size

for ii = rn1:rn2
disp(ii);

% Load presimulated outcome data
filename = strcat(yDir,'/ySim',num2str(ii),'.mat');
load(filename);
y = y(I_,:);            % Only select firms with minimum # of obs.
yFull=NaN*ones(J*T,1);
yFull(IS,:)=y;

% C) Point estimates

% Firm, State X Year
YFMSxT = y - DFMSxT*(DFMSxT\y);
biist = XFMSxT\YFMSxT  ;
bFMSxT(ii,:) = biist';
e = YFMSxT-XFMSxT*biist;
% State
[se,~,~] = cluster_se(XFMSxT,e,MSxT,state,kx);
se_state(ii,:) = (sqrt((n-kx)/(n-kSxT))*se)';

% Firm, sic1 X Year
YFsc1xT = y - DFsc1xT*(DFsc1xT\y);
biis1 = XFsc1xT\YFsc1xT  ;
bFsc1xT(ii,:) = biis1';
e = YFsc1xT-XFsc1xT*biis1;
% one-digit SIC
[se,~,~] = cluster_se(XFsc1xT,e,Msc1xT,sic1,kx);
se_sic1(ii,:) = (sqrt((n-kx)/(n-ksc1xT))*se)';

% Firm, sic2 X Year
YFsc2xT = y - DFsc2xT*(DFsc2xT\y);
biis2 = XFsc2xT\YFsc2xT  ;
bFsc2xT(ii,:) = biis2';
e = YFsc2xT-XFsc2xT*biis2;
% one-digit SIC
[se,~,~] = cluster_se(XFsc2xT,e,Msc2xT,sic2,kx);
se_sic2(ii,:) = (sqrt((n-kx)/(n-ksc2xT))*se)';

% Firm, size X Year
YFszexT = y - DFszexT*(DFszexT\y);
biisz = XFszexT\YFszexT  ;
bFszexT(ii,:) = biisz';
e = YFszexT-XFszexT*biisz;
% one-digit SIC
[se,~,~] = cluster_se(XFszexT,e,MszexT,sizec,kx);
se_size(ii,:) = (sqrt((n-kx)/(n-kszexT))*se)';

%%% D) Moving Bloc Bootsrap 
bb_state = zeros(BB,kx);          % Point estimates
bb_sic1  = zeros(BB,kx);           % Point estimates
bb_sic2  = zeros(BB,kx);           % Point estimates
bb_size  = zeros(BB,kx);           % Point estimates

bbT_state   = zeros(size(bb_state));    % T-statistic cluster firm
bbT_sic1    = zeros(size(bb_sic1));    % T-statistic cluster firm
bbT_sic2    = zeros(size(bb_sic2));    % T-statistic cluster firm
bbT_size    = zeros(size(bb_size));    % T-statistic cluster firm

I = randi([1 T-bsize+1],[BB,ceil(T/bsize)]) ;   % sample indexes
Isel=zeros(T*J,1);

    parfor indb=1:BB
    warning('off','MATLAB:rankDeficientMatrix');
    %%%%%% Moving Block Bootsrap
    Iplus = index_mbb( I(indb,:), bsize, T, J );
    Isel=(f_full -1)*T + Iplus;
    xb=xFull(Isel,:);
    yb=yFull(Isel,:);

    % Get rid of missing observations
    inan  = (isnan(xb)==1);
    Index = (sum(inan,2)==0);
    xb    = xb(Index,:);             
    yb    = yb(Index,:);             
    fb    = f_full(Index,:);      % Firm Index
    tb    = t_full(Index,:);      % Time Index
    Sb    = stFull(Index,:);      % State Index
    sb    = sizeF(Index,:);       % Size Index  
    s1b   = sic1F(Index,:);       % SIC1 Index  
    s2b   = sic2F(Index,:);       % SIC2 index
    nb    = size(xb,1);           % Number of observations 

    %%%%%% Re-construct dummy matrices
    Dbi  = sparse(dummyvar(fb));
    % Firm, State X Year
    Iaux = groupcross(Sb,tb);
    Daux = sparse(dummyvar(Iaux));
    Dst  = [Dbi Daux];
    
    % Firm, one-digit x year
    Iaux = groupcross(s1b,tb);
    Daux = sparse(dummyvar(Iaux));
    Ds1  = [Dbi Daux];
    
    % Firm, two-digit x year
    Iaux = groupcross(s2b,tb);
    Daux = sparse(dummyvar(Iaux));
    Ds2  = [Dbi Daux];    

    % Firm, size x year
    Iaux = groupcross(sb,tb);
    Daux = sparse(dummyvar(Iaux));
    Dsz  = [Dbi Daux];    
    
    %%%%%% point estimates and t-statistics
    % Firm, State X Year  
    Xbst = xb - Dst*(Dst\xb);
    Ybst = yb - Dst*(Dst\yb);    

    bb   = Xbst\Ybst;
    Mbst = inv(Xbst'*Xbst);
    ebst = Ybst-Xbst*bb;
    bb_state(indb,:) = bb';
    
    [sebb,~,~] = cluster_se(Xbst,ebst,Mbst,Sb,kx);
    sebb = (sqrt((nb-kx)/(nb-kSxT))*sebb )';
    bbT_state(indb,:) = (bb-biist)'./sebb;
    
    % Firm, sic1 X Year    
    Xbs1 = xb - Ds1*(Ds1\xb);
    Ybs1 = yb - Ds1*(Ds1\yb);    

    bb   = Xbs1\Ybs1;
    Mbs1 = inv(Xbs1'*Xbs1);
    ebs1 = Ybs1-Xbs1*bb;
    bb_sic1(indb,:) = bb';
    
    [sebb,~,~] = cluster_se(Xbs1,ebs1,Mbs1,s1b,kx);
    sebb = (sqrt((nb-kx)/(nb-ksc1xT))*sebb )';
    bbT_sic1(indb,:) = (bb-biis1)'./sebb;     
    
    % Firm, sic2 X Year    
    Xbs2 = xb - Ds2*(Ds2\xb);
    Ybs2 = yb - Ds2*(Ds2\yb);    

    bb   = Xbs2\Ybs2;
    Mbs2 = inv(Xbs2'*Xbs2);
    ebs2 = Ybs2-Xbs2*bb;
    bb_sic2(indb,:) = bb';
    
    [sebb,~,~] = cluster_se(Xbs2,ebs2,Mbs2,s2b,kx);
    sebb = (sqrt((nb-kx)/(nb-ksc2xT))*sebb )';
    bbT_sic2(indb,:) = (bb-biis2)'./sebb;    
    
    % Firm, sic2 X Year    
    Xbsz = xb - Dsz*(Dsz\xb);
    Ybsz = yb - Dsz*(Dsz\yb);    

    bb   = Xbsz\Ybsz;
    Mbsz = inv(Xbsz'*Xbsz);
    ebsz = Ybsz-Xbsz*bb;
    bb_size(indb,:) = bb';
    
    [sebb,~,~] = cluster_se(Xbsz,ebsz,Mbsz,sb,kx);
    sebb = (sqrt((nb-kx)/(nb-kszexT))*sebb )';
    bbT_size(indb,:) = (bb-biisz)'./sebb;
    
    end

bcv_state(ii,:) =  prctile(abs(bbT_state),95);        % state
bcv_sic1(ii,:)  =  prctile(abs(bbT_sic1),95);         % one-digit SIC
bcv_sic2(ii,:)  =  prctile(abs(bbT_sic2),95);         % two-digit SIC
bcv_size(ii,:)  =  prctile(abs(bbT_size),95);         % Size

fMat=strcat('Results/Pt3/SizePt3_',num2str(nSim),'X',num2str(BB),'_','bsize',num2str(bsize),'_',num2str(minObs),'Obs','_sim_',num2str(ii),'.mat');
save(fMat,'bFMSxT','bFsc1xT','bFsc2xT','bFszexT','se_state','se_sic1','se_sic2','se_size','bcv_state','bcv_sic1','bcv_sic2','bcv_size','bTest3');
end


toc