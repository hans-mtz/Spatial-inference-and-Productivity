clear ;

MainDir = 'C:\Users\chansen1\Dropbox\JARReview\Example_VerdiInvestment' ;
cd(MainDir);

%% Load preestimated parameters for error DGP

load bcv_optim6 ;    % Spatial covariance parameters (not actually used here
                     % but used to generate Sigma
load SigmaEsts ;     % Time by time period spatial covariance matrices
load distparams ;    % Parameters of the iid innovation distribution
load XCatIndices ;   % Indices of categories used to form lagged variables
load TimeIndices ;   % Time indices
load lagparams ;     % Estimated "autoregressive" parameters 
load fixedeffects ;  % Estimates of \alpha_i + \delta_t
load FirmIndices ;   % Firm indices
load xMat ;          % Design matrix from actual data
load olsparams ;     % OLS estimates of parameters on x variables

%% Generate y's for new data

rng(1262015) ;  % Setting used to generate first 1000 datasets is 1262015.  
                % Change to generate new data sets that are not perfectly
                % correlated with the first.
nSim = 1000 ;

yDir = 'C:\Users\chansen1\Dropbox\JARReview\Example_VerdiInvestment\SimulatedOutcomes' ;  
                % Directory where new y's will be stored
                
for ii = 1:nSim
    disp(ii)
    err = FRQSim_GenErr(egb2hatC,SIGMA,aTestM,indX,t_ind);
    y = x*bTest3 + FE + err ;
    filename = strcat(yDir,'\ySim',num2str(ii),'.mat');
    save(filename,'y');
end
