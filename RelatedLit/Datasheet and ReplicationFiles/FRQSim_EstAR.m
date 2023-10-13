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
                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Some data preprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Firm dummies, time dummies
Di = sparse(dummyvar(f_ind));
Dt = sparse(dummyvar(t_ind));
D = [Di Dt];

% Create quintile in x category variables
xCat = zeros(size(x));

for jj = 1:size(x,2)
    gmx = Di*(Di\x(:,jj));
    
    % Make 5 groups based on quantiles
    px = prctile(gmx,20:20:80);

    for pp = 1:5
        if pp == 1
            xCat(gmx <= px(pp),jj) = pp;
        elseif pp > 1 && pp < 5
            xCat(gmx <= px(pp) & gmx > px(pp-1),jj) = pp;
        else
            xCat(gmx > px(pp-1),jj) = pp;
        end
    end
    
end

% Partial out fixed effects for various strategies and form dummy matrices
X = x - D*(D\x);
M = inv(X'*X);

n = size(X,1);
kx = size(X,2);                     

% Other categorical variables that may be useful
state = indX(:,2);
SxT = groupcross(state,t_ind);

sic2 = indX(:,3);
sic1 = recode(floor(sic2/10)+1);
SIC1xT = groupcross(sic1,t_ind);

sizecat = indX(:,4);
SIZExT = groupcross(sizecat,t_ind);

tind4 = ceil(t_ind/4);
tind4(tind4 == 5) = 4;
FxT4 = groupcross(f_ind,tind4);

tind2 = ceil(t_ind/2);
tind2(tind2 == 9) = 8;

% Index variables for constructing cross-sectional averages
indvars = [xCat state sic1 sic2];
ki = size(indvars,2);

% Generate functions of some lagged x variable to get indices for use in
% constructing distance matrices - 3 lags?
tempL1x = zeros(n,ki);
tempL2x = zeros(n,ki);
tempL3x = zeros(n,ki);
for jj = 1:ki
    tempL1x(:,jj) = FRQLagGroupMean(X(:,1),indvars(:,jj),t_ind,1);
    tempL2x(:,jj) = FRQLagGroupMean(X(:,1),indvars(:,jj),t_ind,2);
    tempL3x(:,jj) = FRQLagGroupMean(X(:,1),indvars(:,jj),t_ind,3);
end
tempL1x = [FRQlag(X(:,1),f_ind,t_ind,1),tempL1x];
tempL2x = [FRQlag(X(:,1),f_ind,t_ind,2),tempL2x];
tempL3x = [FRQlag(X(:,1),f_ind,t_ind,3),tempL3x];
tempLx = [tempL1x,tempL2x,tempL3x];

% indices of observations for which three lags in the variables can
% be constructed
use = ~isnan(sum(tempLx,2));

% Get subsets of the data for relevant variables with three lags available
XU = X(use,:);
tU = t_ind(use,:);
stateU = state(use,:);
sic1U = sic1(use,:);
sic2U = sic2(use,:);

% Form distance vectors
for tt = min(tU):max(tU)
    tempInd = tU == tt;
    tempN = sum(tempInd);
    tempV = stateU(tempInd);
    tempSic1 = sic1U(tempInd);
    tempSic2 = sic2U(tempInd);
    tempX = XU(tempInd,:);

    WState = zeros(tempN);
    WSic1 = zeros(tempN);
    WSic2 = zeros(tempN);
    for jj = 1:tempN
        WState(:,jj) = tempV == tempV(jj);
        WSic1(:,jj) = tempSic1 == tempSic1(jj);
        WSic2(:,jj) = tempSic2 == tempSic2(jj);
    end
    VState = WState(logical(tril(ones(tempN),-1)));
    VSic1 = WSic1(logical(tril(ones(tempN),-1)));
    VSic2 = WSic2(logical(tril(ones(tempN),-1)));
    
    WX1 = zeros(tempN);
    WX2 = zeros(tempN);
    WX3 = zeros(tempN);
    WX4 = zeros(tempN);
    WX5 = zeros(tempN);
    WX6 = zeros(tempN);
    WX7 = zeros(tempN);
    WX8 = zeros(tempN);
    WX9 = zeros(tempN);
    for jj = 1:tempN
        WX1(:,jj) = abs(tempX(:,1)-tempX(jj,1));
        WX2(:,jj) = abs(tempX(:,2)-tempX(jj,2));
        WX3(:,jj) = abs(tempX(:,3)-tempX(jj,3));
        WX4(:,jj) = abs(tempX(:,4)-tempX(jj,4));
        WX5(:,jj) = abs(tempX(:,5)-tempX(jj,5));
        WX6(:,jj) = abs(tempX(:,6)-tempX(jj,6));
        WX7(:,jj) = abs(tempX(:,7)-tempX(jj,7));
        WX8(:,jj) = abs(tempX(:,8)-tempX(jj,8));
        WX9(:,jj) = abs(tempX(:,9)-tempX(jj,9));
    end
    VX1 = WX1(logical(tril(ones(tempN),-1)));
    VX2 = WX2(logical(tril(ones(tempN),-1)));
    VX3 = WX3(logical(tril(ones(tempN),-1)));
    VX4 = WX4(logical(tril(ones(tempN),-1)));
    VX5 = WX5(logical(tril(ones(tempN),-1)));
    VX6 = WX6(logical(tril(ones(tempN),-1)));
    VX7 = WX7(logical(tril(ones(tempN),-1)));
    VX8 = WX8(logical(tril(ones(tempN),-1)));
    VX9 = WX9(logical(tril(ones(tempN),-1)));
    
    assignin('base',strcat('VState_',num2str(tt)),VState);
    assignin('base',strcat('VSic1_',num2str(tt)),VSic1);
    assignin('base',strcat('VSic2_',num2str(tt)),VSic2);
    assignin('base',strcat('VX1_',num2str(tt)),VX1);
    assignin('base',strcat('VX2_',num2str(tt)),VX2);
    assignin('base',strcat('VX3_',num2str(tt)),VX3);
    assignin('base',strcat('VX4_',num2str(tt)),VX4);
    assignin('base',strcat('VX5_',num2str(tt)),VX5);
    assignin('base',strcat('VX6_',num2str(tt)),VX6);
    assignin('base',strcat('VX7_',num2str(tt)),VX7);
    assignin('base',strcat('VX8_',num2str(tt)),VX8);
    assignin('base',strcat('VX9_',num2str(tt)),VX9);
end
VSic1 = sparse([VSic1_4 ; VSic1_5 ; VSic1_6 ; VSic1_7 ; VSic1_8 ; VSic1_9 ; ... 
    VSic1_10 ; VSic1_11 ; VSic1_12 ; VSic1_13 ; VSic1_14 ; VSic1_15 ; ... 
    VSic1_16 ; VSic1_17]);    
VSic2 = sparse([VSic2_4 ; VSic2_5 ; VSic2_6 ; VSic2_7 ; VSic2_8 ; VSic2_9 ; ... 
    VSic2_10 ; VSic2_11 ; VSic2_12 ; VSic2_13 ; VSic2_14 ; VSic2_15 ; ... 
    VSic2_16 ; VSic2_17]);    
VState = sparse([VState_4 ; VState_5 ; VState_6 ; VState_7 ; VState_8 ; VState_9 ; ... 
    VState_10 ; VState_11 ; VState_12 ; VState_13 ; VState_14 ; VState_15 ; ... 
    VState_16 ; VState_17]);    
VX1 = [VX1_4 ; VX1_5 ; VX1_6 ; VX1_7 ; VX1_8 ; VX1_9 ; ... 
    VX1_10 ; VX1_11 ; VX1_12 ; VX1_13 ; VX1_14 ; VX1_15 ; ... 
    VX1_16 ; VX1_17];    
VX2 = [VX2_4 ; VX2_5 ; VX2_6 ; VX2_7 ; VX2_8 ; VX2_9 ; ... 
    VX2_10 ; VX2_11 ; VX2_12 ; VX2_13 ; VX2_14 ; VX2_15 ; ... 
    VX2_16 ; VX2_17];    
VX3 = [VX3_4 ; VX3_5 ; VX3_6 ; VX3_7 ; VX3_8 ; VX3_9 ; ... 
    VX3_10 ; VX3_11 ; VX3_12 ; VX3_13 ; VX3_14 ; VX3_15 ; ... 
    VX3_16 ; VX3_17];    
VX4 = [VX4_4 ; VX4_5 ; VX4_6 ; VX4_7 ; VX4_8 ; VX4_9 ; ... 
    VX4_10 ; VX4_11 ; VX4_12 ; VX4_13 ; VX4_14 ; VX4_15 ; ... 
    VX4_16 ; VX4_17];    
VX5 = [VX5_4 ; VX5_5 ; VX5_6 ; VX5_7 ; VX5_8 ; VX5_9 ; ... 
    VX5_10 ; VX5_11 ; VX5_12 ; VX5_13 ; VX5_14 ; VX5_15 ; ... 
    VX5_16 ; VX5_17];    
VX6 = [VX6_4 ; VX6_5 ; VX6_6 ; VX6_7 ; VX6_8 ; VX6_9 ; ... 
    VX6_10 ; VX6_11 ; VX6_12 ; VX6_13 ; VX6_14 ; VX6_15 ; ... 
    VX6_16 ; VX6_17];    
VX7 = [VX7_4 ; VX7_5 ; VX7_6 ; VX7_7 ; VX7_8 ; VX7_9 ; ... 
    VX7_10 ; VX7_11 ; VX7_12 ; VX7_13 ; VX7_14 ; VX7_15 ; ... 
    VX7_16 ; VX7_17];    
VX8 = [VX8_4 ; VX8_5 ; VX8_6 ; VX8_7 ; VX8_8 ; VX8_9 ; ... 
    VX8_10 ; VX8_11 ; VX8_12 ; VX8_13 ; VX8_14 ; VX8_15 ; ... 
    VX8_16 ; VX8_17];    
VX9 = [VX9_4 ; VX9_5 ; VX9_6 ; VX9_7 ; VX9_8 ; VX9_9 ; ... 
    VX9_10 ; VX9_11 ; VX9_12 ; VX9_13 ; VX9_14 ; VX9_15 ; ... 
    VX9_16 ; VX9_17];    

% Create quintile in x-distance category variables
distCat = zeros(size(VX1,1),kx);

for jj = 1:kx
    currvar = strcat('VX',num2str(jj));
    
    % Make 5 groups based on quantiles
    px = prctile(eval(currvar),20:20:80);

    for pp = 1:5
        if pp == 1
            distCat(eval(currvar) <= px(pp),jj) = pp;
        elseif pp > 1 && pp < 5
            distCat(eval(currvar) <= px(pp) & eval(currvar) > px(pp-1),jj) = pp;
        else
            distCat(eval(currvar) > px(pp-1),jj) = pp;
        end
    end
    
end

distMat = [ones(size(VSic1)) VSic1 VSic2 VState];
for jj = 1:kx
    temp = sparse(dummyvar(distCat(:,jj)));
%     distMat = [distMat sparse(dummy(distCat(:,jj)))]; %#ok<AGROW>
    distMat = [distMat temp(:,1:end-1)]; %#ok<AGROW>
end
QdM = distMat'*distMat;
iQdM = inv(distMat'*distMat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Predefine variables for simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of simulation replications
nSim = 1000 ;
b = zeros(nSim,kx);

pac1 = zeros(3*(ki+1),nSim);
spac1 = zeros(3*(ki+1),nSim);

acf1 = zeros(ki+1,nSim);
acf2 = zeros(ki+1,nSim);
acf3 = zeros(ki+1,nSim);
sacf1 = zeros(ki+1,nSim);
sacf2 = zeros(ki+1,nSim);
sacf3 = zeros(ki+1,nSim);

spcoef = zeros(size(QdM,1),nSim);
se_spcoef = zeros(size(QdM,1),nSim);

sptest = zeros(kx+3,nSim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main simulation step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:nSim
    if (ii-1)/20 == floor((ii-1)/20)
        disp(ii);
    end
    
    % Load presimulated outcome data
    filename = strcat(yDir,'\ySim',num2str(ii),'.mat');
    load(filename);

    
    % Baseline estimate
    Y = y - D*(D\y);
    bii = X\Y;
    b(ii,:) = bii';
    e = Y-X*bii;
     
    % Generate functions of lagged residuals
    L11 = zeros(n,ki);
    L12 = zeros(n,ki);
    L13 = zeros(n,ki);
    for jj = 1:ki
        L11(:,jj) = FRQLagGroupMean(e,indvars(:,jj),t_ind,1);
        L12(:,jj) = FRQLagGroupMean(e,indvars(:,jj),t_ind,2);
        L13(:,jj) = FRQLagGroupMean(e,indvars(:,jj),t_ind,3);
    end
    
    L1 = [FRQlag(e,f_ind,t_ind,1),L11,FRQlag(e,f_ind,t_ind,2),L12,...
        FRQlag(e,f_ind,t_ind,3),L13];
    
    % Estimated "autoregressive" parameters
    pac1(:,ii) = L1(use,:)\e(use,1);
    
    v = e(use,1)-L1(use,:)*pac1(:,ii);
    
    spac1(:,ii) = cluster_se(L1(use,:),v,pinv(L1(use,:)'*L1(use,:)),t_ind(use,:));
    
    % Estimated "autocorrelation" function
    acf1(:,ii) = L1(use,1:13)\e(use,1);
    acf2(:,ii) = L1(use,14:26)\e(use,1);
    acf3(:,ii) = L1(use,27:39)\e(use,1);

    sacf1(:,ii) = cluster_se(L1(use,1:13),e(use,1)-L1(use,1:13)*acf1(:,ii),pinv(L1(use,1:13)'*L1(use,1:13)),t_ind(use,:));    
    sacf2(:,ii) = cluster_se(L1(use,14:26),e(use,1)-L1(use,14:26)*acf2(:,ii),pinv(L1(use,14:26)'*L1(use,14:26)),t_ind(use,:));    
    sacf3(:,ii) = cluster_se(L1(use,27:39),e(use,1)-L1(use,27:39)*acf3(:,ii),pinv(L1(use,27:39)'*L1(use,27:39)),t_ind(use,:));    

    % Note that under iid errors and assuming a balanced panel, own
    % autocorrelations should all be -1/(T-1).  Using this as a 
    % back-of-the-envelope benchmark replacing T with \bar{T} because of 
    % unbalanced panel gives -1/(9.86-1) \approx -.1129.  
  
%     tslabel = {'e_1','xbar1_1','xbar2_1','xbar3_1','xbar4_1','xbar5_1',...
%         'xbar6_1','xbar7_1','xbar8_1','xbar9_1','statebar_1','sic1bar_1','sic2bar_1',...
%         'e_2','xbar1_2','xbar2_2','xbar3_2','xbar4_2','xbar5_2',...
%         'xbar6_2','xbar7_2','xbar8_2','xbar9_2','statebar_2','sic1bar_2','sic2bar_2',...
%         'e_3','xbar1_3','xbar2_3','xbar3_3','xbar4_3','xbar5_3',...
%         'xbar6_3','xbar7_3','xbar8_3','xbar9_3','statebar_3','sic1bar_3','sic2bar_3'};
%     space = '     ';
%     space = repmat(space,39,1);
    
    % Construct cross-products of residuals from time series model
    for tt = min(tU):max(tU)
        tempInd = tU == tt;
        tempN = sum(tempInd);
        
        tempV = v(tempInd);
        Wv = tempV*tempV';
        Vv = Wv(logical(tril(ones(tempN),-1)));
        
        assignin('base',strcat('Vv_',num2str(tt)),Vv);
    end
    Vv = [Vv_4 ; Vv_5 ; Vv_6 ; Vv_7 ; Vv_8 ; Vv_9 ; ...
        Vv_10 ; Vv_11 ; Vv_12 ; Vv_13 ; Vv_14 ; Vv_15 ; ...
        Vv_16 ; Vv_17];
    dMV = distMat'*Vv;
    spcoef(:,ii) = iQdM*dMV; %#ok<MINV>
    r = Vv - distMat*spcoef(:,ii);
    V_spcoef = var(r)*iQdM; %#ok<MINV>
    se_spcoef(:,ii) = sqrt(diag(V_spcoef));
    Wx = zeros(kx,1);
    for jj = 1:kx
        Wx(jj) = spcoef(4*(jj-1)+4:4*jj+4,ii)'*inv(V_spcoef(4*(jj-1)+4:4*jj+4,4*(jj-1)+4:4*jj+4))*spcoef(4*(jj-1)+4:4*jj+4,ii); %#ok<MINV>
    end
    
    sptest(:,ii) = [(spcoef(2:4,ii)./se_spcoef(2:4,ii)).^2 ; Wx];
    
%     label = {'c','sic1','sic2','state','1','1','1','1','2','2','2','2',...
%         '3','3','3','3','4','4','4','4','5','5','5','5','6','6','6','6',...
%         '7','7','7','7','8','8','8','8','9','9','9','9'};
%     space = '     ';
%     space = repmat(space,40,1);

end

save BCVResidualRegressions pac1 spac1 acf1 acf2 acf3 sacf1 sacf2 sacf3 spcoef se_spcoef sptest ;

