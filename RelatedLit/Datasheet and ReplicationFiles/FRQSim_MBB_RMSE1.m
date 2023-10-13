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


end

rmse= sqrt(mean((b - ones(nSim,1)*bTest3').^2));

num2str(rmse,'%12.3f')

disp([ 'firm, year                     & ' num2str(rmse(1,1),'%12.3f') ' & '  num2str(rmse(1,2),'%12.3f') ' & ' ...
     num2str(rmse(1,3),'%12.3f') ' & '  num2str(rmse(1,4),'%12.3f') ' & '  num2str(rmse(1,5),'%12.3f') ' & ' ...
     num2str(rmse(1,6),'%12.3f') ' & '  num2str(rmse(1,7),'%12.3f') ' & '  num2str(rmse(1,8),'%12.3f') ' & '  num2str(rmse(1,9),'%12.3f') '\\ '])







toc