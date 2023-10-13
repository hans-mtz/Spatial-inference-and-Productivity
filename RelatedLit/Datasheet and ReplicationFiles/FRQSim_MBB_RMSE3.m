% Moving Block Bootstrap
% June 2017, Aldo Sandoval 
clc
clear
warning('off','MATLAB:rankDeficientMatrix');
yDir = 'SimulatedOutcomes' ; 
tic
% Minimum number of observations per firm
minObs=17;

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
nSim    = 1000 ;           % Number of simulation replications
bsize   = 2;            % Bloc size
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
ind = (sum(DFt2)==0);
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


for ii = 1:nSim
 
disp(ii);

% Load presimulated outcome data
filename = strcat(yDir,'/ySim',num2str(ii),'.mat');
load(filename);
y = y(I_,:);
yFull=NaN*ones(J*T,1);
yFull(IS,:)=y;

%%%% C) Point estimates
% Firm X 2 year block, year
YFt2xT = y - DFt2xT*(DFt2xT\y);
biit2 = XFt2xT\YFt2xT  ;
bFt2xT(ii,:) = biit2';

% Firm X 4 year block, year
YFt4xT = y - DFt4xT*(DFt4xT\y);
biit4 = XFt4xT\YFt4xT  ;
bFt4xT(ii,:) = biit4';

% Firm X 6 year block, year
YFt6xT = y - DFt6xT*(DFt6xT\y);
biit6 = XFt6xT\YFt6xT  ;
bFt6xT(ii,:) = biit6';

% Firm X 8 year block, year
YFt8xT = y - DFt8xT*(DFt8xT\y);
biit8 = XFt8xT\YFt8xT  ;
bFt8xT(ii,:) = biit8';


end

rmse_t2 = sqrt(mean((bFt2xT - ones(nSim,1)*bTest3').^2));
rmse_t4 = sqrt(mean((bFt4xT - ones(nSim,1)*bTest3').^2));
rmse_t6 = sqrt(mean((bFt6xT - ones(nSim,1)*bTest3').^2));
rmse_t8 = sqrt(mean((bFt8xT - ones(nSim,1)*bTest3').^2));
num2str(rmse_t8,'%12.3f')
num2str(rmse_t6,'%12.3f')
num2str(rmse_t4,'%12.3f')
num2str(rmse_t2,'%12.3f')

 disp([ 'firm x 8 year time block, year  & ' num2str(rmse_t8(1,1),'%12.3f') ' & '  num2str(rmse_t8(1,2),'%12.3f') ' & ' ...
     num2str(rmse_t8(1,3),'%12.3f') ' & '  num2str(rmse_t8(1,4),'%12.3f') ' & '  num2str(rmse_t8(1,5),'%12.3f') ' & ' ...
     num2str(rmse_t8(1,6),'%12.3f') ' & '  num2str(rmse_t8(1,7),'%12.3f') ' & '  num2str(rmse_t8(1,8),'%12.3f') ' & '  num2str(rmse_t8(1,9),'%12.3f') '\\ '])
 
  disp([ 'firm x 6 year time block, year  & ' num2str(rmse_t6(1,1),'%12.3f') ' & '  num2str(rmse_t6(1,2),'%12.3f') ' & ' ...
     num2str(rmse_t6(1,3),'%12.3f') ' & '  num2str(rmse_t6(1,4),'%12.3f') ' & '  num2str(rmse_t6(1,5),'%12.3f') ' & ' ...
     num2str(rmse_t6(1,6),'%12.3f') ' & '  num2str(rmse_t6(1,7),'%12.3f') ' & '  num2str(rmse_t6(1,8),'%12.3f') ' & '  num2str(rmse_t6(1,9),'%12.3f') '\\ '])
 
  disp([ 'firm x 4 year time block, year  & ' num2str(rmse_t4(1,1),'%12.3f') ' & '  num2str(rmse_t4(1,2),'%12.3f') ' & ' ...
     num2str(rmse_t4(1,3),'%12.3f') ' & '  num2str(rmse_t4(1,4),'%12.3f') ' & '  num2str(rmse_t4(1,5),'%12.3f') ' & ' ...
     num2str(rmse_t4(1,6),'%12.3f') ' & '  num2str(rmse_t4(1,7),'%12.3f') ' & '  num2str(rmse_t4(1,8),'%12.3f') ' & '  num2str(rmse_t4(1,9),'%12.3f') '\\ '])
 
disp([ 'firm x 2 year time block, year  & ' num2str(rmse_t2(1,1),'%12.3f') ' & '  num2str(rmse_t2(1,2),'%12.3f') ' & ' ...
     num2str(rmse_t2(1,3),'%12.3f') ' & '  num2str(rmse_t2(1,4),'%12.3f') ' & '  num2str(rmse_t2(1,5),'%12.3f') ' & ' ...
     num2str(rmse_t2(1,6),'%12.3f') ' & '  num2str(rmse_t2(1,7),'%12.3f') ' & '  num2str(rmse_t2(1,8),'%12.3f') ' & '  num2str(rmse_t2(1,9),'%12.3f') '\\ '])
 


toc
