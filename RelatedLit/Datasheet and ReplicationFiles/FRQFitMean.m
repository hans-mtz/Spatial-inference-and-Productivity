clear ;

cd C:\Users\chansen1\Dropbox\JARReview\Example_VerdiInvestment ;

%% Read in data and define variables

DATA = importdata('FRQMainMat.csv',',',1);

firmid = (DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'gvkey'), DATA.colheaders, 'UniformOutput', false))));
year = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'fyear'), DATA.colheaders, 'UniformOutput', false)));
state = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'statenum'), DATA.colheaders, 'UniformOutput', false)));
sic1 = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'sic1'), DATA.colheaders, 'UniformOutput', false)));
sic2 = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'sic2'), DATA.colheaders, 'UniformOutput', false)));
sic3 = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'sic3'), DATA.colheaders, 'UniformOutput', false)));
invest = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'inv'), DATA.colheaders, 'UniformOutput', false)));
res = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'REAL_ESTATE_STATE'), DATA.colheaders, 'UniformOutput', false)));
ind = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'STATE_INDEX'), DATA.colheaders, 'UniformOutput', false)));
frq = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'std_frqind'), DATA.colheaders, 'UniformOutput', false)));
frqint = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'int_frqind'), DATA.colheaders, 'UniformOutput', false)));
cash = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'cash1'), DATA.colheaders, 'UniformOutput', false)));
tq = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'lag_q'), DATA.colheaders, 'UniformOutput', false)));
mve = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'ln_lag_mve'), DATA.colheaders, 'UniformOutput', false)));
lev = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'lag_leverage'), DATA.colheaders, 'UniformOutput', false)));
age = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'ln_lag_age'), DATA.colheaders, 'UniformOutput', false)));
resid0 = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'resid'), DATA.colheaders, 'UniformOutput', false))); 
mstate = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'mstate'), DATA.colheaders, 'UniformOutput', false))); 
nstate = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'nstate'), DATA.colheaders, 'UniformOutput', false))); 
msic2 = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'msic'), DATA.colheaders, 'UniformOutput', false))); 
nsic2 = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'nsic'), DATA.colheaders, 'UniformOutput', false))); 
msize = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'msize'), DATA.colheaders, 'UniformOutput', false))); 
nsize = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'nsize'), DATA.colheaders, 'UniformOutput', false))); 
mcash = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'mcash'), DATA.colheaders, 'UniformOutput', false))); 
ncash = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'ncash'), DATA.colheaders, 'UniformOutput', false))); 
mres = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'mres'), DATA.colheaders, 'UniformOutput', false))); 
nres = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'nres'), DATA.colheaders, 'UniformOutput', false))); 
mfrq = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'mfrq'), DATA.colheaders, 'UniformOutput', false))); 
nfrq = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'nfrq'), DATA.colheaders, 'UniformOutput', false))); 
mint = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'mint'), DATA.colheaders, 'UniformOutput', false))); 
nint = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'nint'), DATA.colheaders, 'UniformOutput', false))); 
flresid = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'flresid'), DATA.colheaders, 'UniformOutput', false))); 
flmstate = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'flmstate'), DATA.colheaders, 'UniformOutput', false))); 
flmsic2 = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'flmsic'), DATA.colheaders, 'UniformOutput', false))); 
flmsize = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'flmsize'), DATA.colheaders, 'UniformOutput', false))); 
flmcash = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'flmcash'), DATA.colheaders, 'UniformOutput', false))); 
flmres = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'flmres'), DATA.colheaders, 'UniformOutput', false))); 
flmfrq = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'flmfrq'), DATA.colheaders, 'UniformOutput', false))); 
flmint = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'flmint'), DATA.colheaders, 'UniformOutput', false))); 
vhat = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'vhat'), DATA.colheaders, 'UniformOutput', false))); 
d_inv = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'d_inv'), DATA.colheaders, 'UniformOutput', false)));
d_res = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'d_REAL_ESTATE_STATE'), DATA.colheaders, 'UniformOutput', false)));
d_ind = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'d_STATE_INDEX'), DATA.colheaders, 'UniformOutput', false)));
d_frq = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'d_std_frqind'), DATA.colheaders, 'UniformOutput', false)));
d_int = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'d_int_frqind'), DATA.colheaders, 'UniformOutput', false)));
d_cash = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'d_cash1'), DATA.colheaders, 'UniformOutput', false)));
d_tq = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'d_lag_q'), DATA.colheaders, 'UniformOutput', false)));
d_mve = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'d_ln_lag_mve'), DATA.colheaders, 'UniformOutput', false)));
d_lev = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'d_lag_leverage'), DATA.colheaders, 'UniformOutput', false)));
d_age = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'d_ln_lag_age'), DATA.colheaders, 'UniformOutput', false)));
d_yr2 = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'d_yr2'), DATA.colheaders, 'UniformOutput', false)));
d_yr3 = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'d_yr3'), DATA.colheaders, 'UniformOutput', false)));
d_yr4 = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'d_yr4'), DATA.colheaders, 'UniformOutput', false)));
d_yr5 = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'d_yr5'), DATA.colheaders, 'UniformOutput', false)));
d_yr6 = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'d_yr6'), DATA.colheaders, 'UniformOutput', false)));
d_yr7 = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'d_yr7'), DATA.colheaders, 'UniformOutput', false)));
d_yr8 = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'d_yr8'), DATA.colheaders, 'UniformOutput', false)));
d_yr9 = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'d_yr9'), DATA.colheaders, 'UniformOutput', false)));
d_yr10 = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'d_yr10'), DATA.colheaders, 'UniformOutput', false)));
d_yr11 = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'d_yr11'), DATA.colheaders, 'UniformOutput', false)));
d_yr12 = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'d_yr12'), DATA.colheaders, 'UniformOutput', false)));
d_yr13 = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'d_yr13'), DATA.colheaders, 'UniformOutput', false)));
d_yr14 = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'d_yr14'), DATA.colheaders, 'UniformOutput', false)));
d_yr15 = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'d_yr15'), DATA.colheaders, 'UniformOutput', false)));
d_yr16 = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'d_yr16'), DATA.colheaders, 'UniformOutput', false)));
d_yr17 = DATA.data(:,...
    cell2mat(cellfun(@(x) strcmp(x,'d_yr17'), DATA.colheaders, 'UniformOutput', false)));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make sure we get the same results as we got in Stata
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% First/Main regression
bStata = [2.199 ; -.1304 ; .0541 ; -.2478 ; .0903 ; 1.1072 ; .3223 ; -4.704 ; -.7741];

% Demeaned data imported from Stata
dY = d_inv;
dX = [d_res d_ind d_frq d_int d_cash d_tq d_mve d_lev d_age ...
    d_yr2 d_yr3 d_yr4 d_yr5 d_yr6 d_yr7 d_yr8 d_yr9 d_yr10 ...
    d_yr11 d_yr12 d_yr13 d_yr14 d_yr15 d_yr16 d_yr17];
bTest1 = dX\dY;

% Data matrices
y = invest;
x = [res ind frq frqint cash tq mve lev age];

% Partial out firm effects by within firm demeaning
f_ind = recode(firmid);
t_ind = recode(year);

Dt = dummyvar(t_ind);
X = [x Dt(:,2:end)];
Y = y;
for ii = 1:max(f_ind)
    fii = f_ind == ii;
    Y(fii) = Y(fii) - mean(Y(fii));
    X(fii,:) = X(fii,:) - ones(sum(fii),1)*mean(X(fii,:));
end
bTest2 = X\Y;

% Brute force including all the dummies
Di = dummyvar(f_ind);

Mx = x - [Di Dt(:,2:end)]*([Di Dt(:,2:end)]\x);
My = y - [Di Dt(:,2:end)]*([Di Dt(:,2:end)]\y);

bTest3 = Mx\My;
corr([bStata,bTest1(1:size(Mx,2)),bTest2(1:size(Mx,2)),bTest3])

resid1 = My-Mx*bTest3;
corr([resid0 resid1])

save FirmIndices f_ind ;
save xMat x ;
save olsparams bTest3 ;

FE = zeros(size(y)); 
for ii = 1:max(f_ind)
    fii = f_ind == ii;
    FE(fii) = ones(sum(fii),1)*(mean(y(fii)) - mean([x(fii,:) Dt(fii,2:end)])*bTest2);
end
FE = FE+[zeros(size(x)) Dt(:,2:end)]*bTest2;

save fixedeffects FE ;


%% Auxiliary regression using "filled" variables
aStata = [.2784 ; .0309 ; .2027 ; .2514 ; .305 ; .2350 ; .2781 ; -.2445];

flx = [flresid flmstate flmsic2 flmsize flmcash flmres flmfrq flmint];
aTest = flx\resid1;

corr([aStata aTest])

vhat1 = resid1-flx*aTest;
corr([vhat vhat1])

%% Form variables for auxiliary regression in MATLAB
% Firm level averages of a few variables to form x-dependent categories
gmsize = Di*(Di\mve);
gmcash = Di*(Di\cash);
gmres = Di*(Di\res);
gmfrq = Di*(Di\frq);
gmint = Di*(Di\frqint);

% Make 5 groups based on quantiles
psize = prctile(gmsize,20:20:80);
pcash = prctile(gmcash,20:20:80);
pres = prctile(gmres,20:20:80);
pfrq = prctile(gmfrq,20:20:80);
pint = prctile(gmint,20:20:80);

% Create categorical variables
sizecat = zeros(size(gmsize));
cashcat = zeros(size(gmsize));
rescat = zeros(size(gmsize));
frqcat = zeros(size(gmsize));
intcat = zeros(size(gmsize));
for pp = 1:5
    if pp == 1
        sizecat(gmsize <= psize(pp))= pp;
        cashcat(gmcash <= pcash(pp)) = pp;
        rescat(gmres <= pres(pp)) = pp;
        frqcat(gmfrq <= pfrq(pp)) = pp;
        intcat(gmint <= pint(pp)) = pp;
    elseif pp > 1 && pp < 5
        sizecat(gmsize <= psize(pp) & gmsize > psize(pp-1)) = pp;
        cashcat(gmcash <= pcash(pp) & gmcash > pcash(pp-1)) = pp;
        rescat(gmres <= pres(pp) & gmres > pres(pp-1)) = pp;
        frqcat(gmfrq <= pfrq(pp) & gmfrq > pfrq(pp-1)) = pp;
        intcat(gmint <= pint(pp) & gmint > pint(pp-1)) = pp;
    else
        sizecat(gmsize > psize(pp-1)) = pp;
        cashcat(gmcash > pcash(pp-1)) = pp;
        rescat(gmres > pres(pp-1)) = pp;
        frqcat(gmfrq > pfrq(pp-1)) = pp;
        intcat(gmint > pint(pp-1)) = pp;
    end
end

% Note:  This breaks "ties" for discreteness in a different way than Stata
% so answers won't be identical.

% Want means of lagged residuals within category x year groups. Categorical
% variables considered: firm, state, sic2, sizecat, cashcat, rescat, frqcat,
% intcat

flresidM = FRQFillLagGroupMean(resid1,f_ind,t_ind); 
flmstateM = FRQFillLagGroupMean(resid1,state,t_ind); 
flmsic2M = FRQFillLagGroupMean(resid1,sic2,t_ind); 
flmsizeM = FRQFillLagGroupMean(resid1,sizecat,t_ind); 
flmcashM = FRQFillLagGroupMean(resid1,cashcat,t_ind); 
flmresM = FRQFillLagGroupMean(resid1,rescat,t_ind); 
flmfrqM = FRQFillLagGroupMean(resid1,frqcat,t_ind); 
flmintM = FRQFillLagGroupMean(resid1,intcat,t_ind); 

indX = [f_ind state sic2 sizecat cashcat rescat frqcat intcat];

save XCatIndices indX ;
save TimeIndices t_ind ;

%% Auxiliary regression using variables constructed in MATLAB

flxM = [flresidM flmstateM flmsic2M flmsizeM flmcashM flmresM flmfrqM flmintM];
aTestM = flxM\resid1;

corr([aStata aTest aTestM])

vhat1M = resid1-flxM*aTestM;
corr([vhat vhat1 vhat1M])

save lagparams aTestM ;

% Close enough for me.

