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


% Firm, sic1 X Year
YFsc1xT = y - DFsc1xT*(DFsc1xT\y);
biis1 = XFsc1xT\YFsc1xT  ;
bFsc1xT(ii,:) = biis1';


% Firm, sic2 X Year
YFsc2xT = y - DFsc2xT*(DFsc2xT\y);
biis2 = XFsc2xT\YFsc2xT  ;
bFsc2xT(ii,:) = biis2';


% Firm, size X Year
YFszexT = y - DFszexT*(DFszexT\y);
biisz = XFszexT\YFszexT  ;
bFszexT(ii,:) = biisz';


end


rmse_state = sqrt(mean((bFMSxT - ones(nSim,1)*bTest3').^2));
rmse_sic1  = sqrt(mean((bFsc1xT - ones(nSim,1)*bTest3').^2));
rmse_sic2  = sqrt(mean((bFsc2xT - ones(nSim,1)*bTest3').^2));
rmse_size  = sqrt(mean((bFszexT - ones(nSim,1)*bTest3').^2));

num2str(rmse_state,'%12.3f')
num2str(rmse_sic1,'%12.3f')
num2str(rmse_sic2,'%12.3f')
num2str(rmse_size,'%12.3f')

disp([ 'firm, state x year                     & ' num2str(rmse_state(1,1),'%12.3f') ' & '  num2str(rmse_state(1,2),'%12.3f') ' & ' ...
     num2str(rmse_state(1,3),'%12.3f') ' & '  num2str(rmse_state(1,4),'%12.3f') ' & '  num2str(rmse_state(1,5),'%12.3f') ' & ' ...
     num2str(rmse_state(1,6),'%12.3f') ' & '  num2str(rmse_state(1,7),'%12.3f') ' & '  num2str(rmse_state(1,8),'%12.3f') ' & '  num2str(rmse_state(1,9),'%12.3f') '\\ '])
 
disp([ 'firm, one-digit SIC x year                     & ' num2str(rmse_sic1(1,1),'%12.3f') ' & '  num2str(rmse_sic1(1,2),'%12.3f') ' & ' ...
     num2str(rmse_sic1(1,3),'%12.3f') ' & '  num2str(rmse_sic1(1,4),'%12.3f') ' & '  num2str(rmse_sic1(1,5),'%12.3f') ' & ' ...
     num2str(rmse_sic1(1,6),'%12.3f') ' & '  num2str(rmse_sic1(1,7),'%12.3f') ' & '  num2str(rmse_sic1(1,8),'%12.3f') ' & '  num2str(rmse_sic1(1,9),'%12.3f') '\\ '])

 disp([ 'firm, two-digit SIC x year                     & ' num2str(rmse_sic2(1,1),'%12.3f') ' & '  num2str(rmse_sic2(1,2),'%12.3f') ' & ' ...
     num2str(rmse_sic2(1,3),'%12.3f') ' & '  num2str(rmse_sic2(1,4),'%12.3f') ' & '  num2str(rmse_sic2(1,5),'%12.3f') ' & ' ...
     num2str(rmse_sic2(1,6),'%12.3f') ' & '  num2str(rmse_sic2(1,7),'%12.3f') ' & '  num2str(rmse_sic2(1,8),'%12.3f') ' & '  num2str(rmse_sic2(1,9),'%12.3f') '\\ '])
 
 disp([ 'firm, size category x year                     & ' num2str(rmse_size(1,1),'%12.3f') ' & '  num2str(rmse_size(1,2),'%12.3f') ' & ' ...
     num2str(rmse_size(1,3),'%12.3f') ' & '  num2str(rmse_size(1,4),'%12.3f') ' & '  num2str(rmse_size(1,5),'%12.3f') ' & ' ...
     num2str(rmse_size(1,6),'%12.3f') ' & '  num2str(rmse_size(1,7),'%12.3f') ' & '  num2str(rmse_size(1,8),'%12.3f') ' & '  num2str(rmse_size(1,9),'%12.3f') '\\ '])
 

toc