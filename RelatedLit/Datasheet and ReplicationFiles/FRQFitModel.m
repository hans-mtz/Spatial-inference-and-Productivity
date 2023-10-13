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


%% Auxiliary regression using "filled" variables
aStata = [.2784 ; .0309 ; .2027 ; .2514 ; .305 ; .2350 ; .2781 ; -.2445];

flx = [flresid flmstate flmsic2 flmsize flmcash flmres flmfrq flmint];
aTest = flx\resid1;

corr([aStata aTest])

vhat1 = resid1-flx*aTest;
corr([vhat vhat1])

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Same results as in Stata.  Now the "fun" part.  Need to build a 
% time period by time period spatial model for the residuals from
% the auxiliary regression.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start by forming Nt x Nt spatial weight matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Will want some dummy variables
Dst = dummyvar(recode(state));
Ds1 = dummyvar(recode(sic1));
Ds2 = dummyvar(recode(sic2));
Ds3 = dummyvar(recode(sic3));

%% Time period 1
tt = t_ind == 1;
N1 = sum(tt);

% Data for time period 1
Dst_1 = Dst(tt,:);
Ds1_1 = Ds1(tt,:);
Ds2_1 = Ds2(tt,:);
Ds3_1 = Ds3(tt,:);

Dt_1 = Dt(tt,:);
x_1 = x(tt,:);

mve_1 = mve(tt,1);
cash_1 = cash(tt,1);
res_1 = res(tt,1);
frq_1 = frq(tt,1);
frqint_1 = frqint(tt,1);

v_1 = vhat1(tt,1);

% State level "random effect"
WSt_1 = Dst_1*Dst_1';

% SIC1 level "random effect"
WS1_1 = Ds1_1*Ds1_1';

% SIC2 level "random effect"
WS2_1 = Ds2_1*Ds2_1';

% SIC3 level "random effect"
WS3_1 = Ds3_1*Ds3_1';

% Distance matrices
mveDist_1 = zeros(N1);
cashDist_1 = zeros(N1);
frqintDist_1 = zeros(N1);
for ii = 1:N1
    mveDist_1(:,ii) = abs(mve_1 - mve_1(ii));
    cashDist_1(:,ii) = abs(cash_1 - cash_1(ii));
    frqintDist_1(:,ii) = sqrt((res_1 - res_1(ii)).^2 + (frq_1 - frq_1(ii)).^2 ...
        + (frqint_1 - frqint_1(ii)).^2);
end

% Characteristic "neighbors"
Wmve_1 = (1-.5*mveDist_1).*(.5*mveDist_1 < 1) - eye(N1);
Wcash_1 = (1-cashDist_1).*(cashDist_1 < 1) - eye(N1);
Wfrqint_1 = (1-.5*frqintDist_1).*(.5*frqintDist_1 < 1) - eye(N1);


%% Time period 2
tt = t_ind == 2;
N2 = sum(tt);

% Data for time period 2
Dst_2 = Dst(tt,:);
Ds1_2 = Ds1(tt,:);
Ds2_2 = Ds2(tt,:);
Ds3_2 = Ds3(tt,:);

Dt_2 = Dt(tt,:);
x_2 = x(tt,:);

mve_2 = mve(tt,1);
cash_2 = cash(tt,1);
res_2 = res(tt,1);
frq_2 = frq(tt,1);
frqint_2 = frqint(tt,1);

v_2 = vhat1(tt,1);

% State level "random effect"
WSt_2 = Dst_2*Dst_2';

% SIC1 level "random effect"
WS1_2 = Ds1_2*Ds1_2';

% SIC2 level "random effect"
WS2_2 = Ds2_2*Ds2_2';

% SIC3 level "random effect"
WS3_2 = Ds3_2*Ds3_2';

% Distance matrices
mveDist_2 = zeros(N2);
cashDist_2 = zeros(N2);
frqintDist_2 = zeros(N2);
for ii = 1:N2
    mveDist_2(:,ii) = abs(mve_2 - mve_2(ii));
    cashDist_2(:,ii) = abs(cash_2 - cash_2(ii));
    frqintDist_2(:,ii) = sqrt((res_2 - res_2(ii)).^2 + (frq_2 - frq_2(ii)).^2 ...
        + (frqint_2 - frqint_2(ii)).^2);
end

% Characteristic "neighbors"
Wmve_2 = (1-.5*mveDist_2).*(.5*mveDist_2 < 1) - eye(N2);
Wcash_2 = (1-cashDist_2).*(cashDist_2 < 1) - eye(N2);
Wfrqint_2 = (1-.5*frqintDist_2).*(.5*frqintDist_2 < 1) - eye(N2);


%% Time period 3
tt = t_ind == 3;
N3 = sum(tt);

% Data for time period 3
Dst_3 = Dst(tt,:);
Ds1_3 = Ds1(tt,:);
Ds2_3 = Ds2(tt,:);
Ds3_3 = Ds3(tt,:);

Dt_3 = Dt(tt,:);
x_3 = x(tt,:);

mve_3 = mve(tt,1);
cash_3 = cash(tt,1);
res_3 = res(tt,1);
frq_3 = frq(tt,1);
frqint_3 = frqint(tt,1);

v_3 = vhat1(tt,1);

% State level "random effect"
WSt_3 = Dst_3*Dst_3';

% SIC1 level "random effect"
WS1_3 = Ds1_3*Ds1_3';

% SIC2 level "random effect"
WS2_3 = Ds2_3*Ds2_3';

% SIC3 level "random effect"
WS3_3 = Ds3_3*Ds3_3';

% Distance matrices
mveDist_3 = zeros(N3);
cashDist_3 = zeros(N3);
frqintDist_3 = zeros(N3);
for ii = 1:N3
    mveDist_3(:,ii) = abs(mve_3 - mve_3(ii));
    cashDist_3(:,ii) = abs(cash_3 - cash_3(ii));
    frqintDist_3(:,ii) = sqrt((res_3 - res_3(ii)).^2 + (frq_3 - frq_3(ii)).^2 ...
        + (frqint_3 - frqint_3(ii)).^2);
end

% Characteristic "neighbors"
Wmve_3 = (1-.5*mveDist_3).*(.5*mveDist_3 < 1) - eye(N3);
Wcash_3 = (1-cashDist_3).*(cashDist_3 < 1) - eye(N3);
Wfrqint_3 = (1-.5*frqintDist_3).*(.5*frqintDist_3 < 1) - eye(N3);


%% Time period 4
tt = t_ind == 4;
N4 = sum(tt);

% Data for time period 4
Dst_4 = Dst(tt,:);
Ds1_4 = Ds1(tt,:);
Ds2_4 = Ds2(tt,:);
Ds3_4 = Ds3(tt,:);

Dt_4 = Dt(tt,:);
x_4 = x(tt,:);

mve_4 = mve(tt,1);
cash_4 = cash(tt,1);
res_4 = res(tt,1);
frq_4 = frq(tt,1);
frqint_4 = frqint(tt,1);

v_4 = vhat1(tt,1);

% State level "random effect"
WSt_4 = Dst_4*Dst_4';

% SIC1 level "random effect"
WS1_4 = Ds1_4*Ds1_4';

% SIC2 level "random effect"
WS2_4 = Ds2_4*Ds2_4';

% SIC3 level "random effect"
WS3_4 = Ds3_4*Ds3_4';

% Distance matrices
mveDist_4 = zeros(N4);
cashDist_4 = zeros(N4);
frqintDist_4 = zeros(N4);
for ii = 1:N4
    mveDist_4(:,ii) = abs(mve_4 - mve_4(ii));
    cashDist_4(:,ii) = abs(cash_4 - cash_4(ii));
    frqintDist_4(:,ii) = sqrt((res_4 - res_4(ii)).^2 + (frq_4 - frq_4(ii)).^2 ...
        + (frqint_4 - frqint_4(ii)).^2);
end

% Characteristic "neighbors"
Wmve_4 = (1-.5*mveDist_4).*(.5*mveDist_4 < 1) - eye(N4);
Wcash_4 = (1-cashDist_4).*(cashDist_4 < 1) - eye(N4);
Wfrqint_4 = (1-.5*frqintDist_4).*(.5*frqintDist_4 < 1) - eye(N4);


%% Time period 5
tt = t_ind == 5;
N5 = sum(tt);

% Data for time period 5
Dst_5 = Dst(tt,:);
Ds1_5 = Ds1(tt,:);
Ds2_5 = Ds2(tt,:);
Ds3_5 = Ds3(tt,:);

Dt_5 = Dt(tt,:);
x_5 = x(tt,:);

mve_5 = mve(tt,1);
cash_5 = cash(tt,1);
res_5 = res(tt,1);
frq_5 = frq(tt,1);
frqint_5 = frqint(tt,1);

v_5 = vhat1(tt,1);

% State level "random effect"
WSt_5 = Dst_5*Dst_5';

% SIC1 level "random effect"
WS1_5 = Ds1_5*Ds1_5';

% SIC2 level "random effect"
WS2_5 = Ds2_5*Ds2_5';

% SIC3 level "random effect"
WS3_5 = Ds3_5*Ds3_5';

% Distance matrices
mveDist_5 = zeros(N5);
cashDist_5 = zeros(N5);
frqintDist_5 = zeros(N5);
for ii = 1:N5
    mveDist_5(:,ii) = abs(mve_5 - mve_5(ii));
    cashDist_5(:,ii) = abs(cash_5 - cash_5(ii));
    frqintDist_5(:,ii) = sqrt((res_5 - res_5(ii)).^2 + (frq_5 - frq_5(ii)).^2 ...
        + (frqint_5 - frqint_5(ii)).^2);
end

% Characteristic "neighbors"
Wmve_5 = (1-.5*mveDist_5).*(.5*mveDist_5 < 1) - eye(N5);
Wcash_5 = (1-cashDist_5).*(cashDist_5 < 1) - eye(N5);
Wfrqint_5 = (1-.5*frqintDist_5).*(.5*frqintDist_5 < 1) - eye(N5);


%% Time period 6
tt = t_ind == 6;
N6 = sum(tt);

% Data for time period 6
Dst_6 = Dst(tt,:);
Ds1_6 = Ds1(tt,:);
Ds2_6 = Ds2(tt,:);
Ds3_6 = Ds3(tt,:);

Dt_6 = Dt(tt,:);
x_6 = x(tt,:);

mve_6 = mve(tt,1);
cash_6 = cash(tt,1);
res_6 = res(tt,1);
frq_6 = frq(tt,1);
frqint_6 = frqint(tt,1);

v_6 = vhat1(tt,1);

% State level "random effect"
WSt_6 = Dst_6*Dst_6';

% SIC1 level "random effect"
WS1_6 = Ds1_6*Ds1_6';

% SIC2 level "random effect"
WS2_6 = Ds2_6*Ds2_6';

% SIC3 level "random effect"
WS3_6 = Ds3_6*Ds3_6';

% Distance matrices
mveDist_6 = zeros(N6);
cashDist_6 = zeros(N6);
frqintDist_6 = zeros(N6);
for ii = 1:N6
    mveDist_6(:,ii) = abs(mve_6 - mve_6(ii));
    cashDist_6(:,ii) = abs(cash_6 - cash_6(ii));
    frqintDist_6(:,ii) = sqrt((res_6 - res_6(ii)).^2 + (frq_6 - frq_6(ii)).^2 ...
        + (frqint_6 - frqint_6(ii)).^2);
end

% Characteristic "neighbors"
Wmve_6 = (1-.6*mveDist_6).*(.6*mveDist_6 < 1) - eye(N6);
Wcash_6 = (1-cashDist_6).*(cashDist_6 < 1) - eye(N6);
Wfrqint_6 = (1-.6*frqintDist_6).*(.6*frqintDist_6 < 1) - eye(N6);



%% Time period 7
tt = t_ind == 7;
N7 = sum(tt);

% Data for time period 7
Dst_7 = Dst(tt,:);
Ds1_7 = Ds1(tt,:);
Ds2_7 = Ds2(tt,:);
Ds3_7 = Ds3(tt,:);

Dt_7 = Dt(tt,:);
x_7 = x(tt,:);

mve_7 = mve(tt,1);
cash_7 = cash(tt,1);
res_7 = res(tt,1);
frq_7 = frq(tt,1);
frqint_7 = frqint(tt,1);

v_7 = vhat1(tt,1);

% State level "random effect"
WSt_7 = Dst_7*Dst_7';

% SIC1 level "random effect"
WS1_7 = Ds1_7*Ds1_7';

% SIC2 level "random effect"
WS2_7 = Ds2_7*Ds2_7';

% SIC3 level "random effect"
WS3_7 = Ds3_7*Ds3_7';

% Distance matrices
mveDist_7 = zeros(N7);
cashDist_7 = zeros(N7);
frqintDist_7 = zeros(N7);
for ii = 1:N7
    mveDist_7(:,ii) = abs(mve_7 - mve_7(ii));
    cashDist_7(:,ii) = abs(cash_7 - cash_7(ii));
    frqintDist_7(:,ii) = sqrt((res_7 - res_7(ii)).^2 + (frq_7 - frq_7(ii)).^2 ...
        + (frqint_7 - frqint_7(ii)).^2);
end

% Characteristic "neighbors"
Wmve_7 = (1-.7*mveDist_7).*(.7*mveDist_7 < 1) - eye(N7);
Wcash_7 = (1-cashDist_7).*(cashDist_7 < 1) - eye(N7);
Wfrqint_7 = (1-.7*frqintDist_7).*(.7*frqintDist_7 < 1) - eye(N7);


%% Time period 8
tt = t_ind == 8;
N8 = sum(tt);

% Data for time period 8
Dst_8 = Dst(tt,:);
Ds1_8 = Ds1(tt,:);
Ds2_8 = Ds2(tt,:);
Ds3_8 = Ds3(tt,:);

Dt_8 = Dt(tt,:);
x_8 = x(tt,:);

mve_8 = mve(tt,1);
cash_8 = cash(tt,1);
res_8 = res(tt,1);
frq_8 = frq(tt,1);
frqint_8 = frqint(tt,1);

v_8 = vhat1(tt,1);

% State level "random effect"
WSt_8 = Dst_8*Dst_8';

% SIC1 level "random effect"
WS1_8 = Ds1_8*Ds1_8';

% SIC2 level "random effect"
WS2_8 = Ds2_8*Ds2_8';

% SIC3 level "random effect"
WS3_8 = Ds3_8*Ds3_8';

% Distance matrices
mveDist_8 = zeros(N8);
cashDist_8 = zeros(N8);
frqintDist_8 = zeros(N8);
for ii = 1:N8
    mveDist_8(:,ii) = abs(mve_8 - mve_8(ii));
    cashDist_8(:,ii) = abs(cash_8 - cash_8(ii));
    frqintDist_8(:,ii) = sqrt((res_8 - res_8(ii)).^2 + (frq_8 - frq_8(ii)).^2 ...
        + (frqint_8 - frqint_8(ii)).^2);
end

% Characteristic "neighbors"
Wmve_8 = (1-.8*mveDist_8).*(.8*mveDist_8 < 1) - eye(N8);
Wcash_8 = (1-cashDist_8).*(cashDist_8 < 1) - eye(N8);
Wfrqint_8 = (1-.8*frqintDist_8).*(.8*frqintDist_8 < 1) - eye(N8);

%% Time period 9
tt = t_ind == 9;
N9 = sum(tt);

% Data for time period 9
Dst_9 = Dst(tt,:);
Ds1_9 = Ds1(tt,:);
Ds2_9 = Ds2(tt,:);
Ds3_9 = Ds3(tt,:);

Dt_9 = Dt(tt,:);
x_9 = x(tt,:);

mve_9 = mve(tt,1);
cash_9 = cash(tt,1);
res_9 = res(tt,1);
frq_9 = frq(tt,1);
frqint_9 = frqint(tt,1);

v_9 = vhat1(tt,1);

% State level "random effect"
WSt_9 = Dst_9*Dst_9';

% SIC1 level "random effect"
WS1_9 = Ds1_9*Ds1_9';

% SIC2 level "random effect"
WS2_9 = Ds2_9*Ds2_9';

% SIC3 level "random effect"
WS3_9 = Ds3_9*Ds3_9';

% Distance matrices
mveDist_9 = zeros(N9);
cashDist_9 = zeros(N9);
frqintDist_9 = zeros(N9);
for ii = 1:N9
    mveDist_9(:,ii) = abs(mve_9 - mve_9(ii));
    cashDist_9(:,ii) = abs(cash_9 - cash_9(ii));
    frqintDist_9(:,ii) = sqrt((res_9 - res_9(ii)).^2 + (frq_9 - frq_9(ii)).^2 ...
        + (frqint_9 - frqint_9(ii)).^2);
end

% Characteristic "neighbors"
Wmve_9 = (1-.9*mveDist_9).*(.9*mveDist_9 < 1) - eye(N9);
Wcash_9 = (1-cashDist_9).*(cashDist_9 < 1) - eye(N9);
Wfrqint_9 = (1-.9*frqintDist_9).*(.9*frqintDist_9 < 1) - eye(N9);


%% Time period 10
tt = t_ind == 10;
N10 = sum(tt);

% Data for time period 10
Dst_10 = Dst(tt,:);
Ds1_10 = Ds1(tt,:);
Ds2_10 = Ds2(tt,:);
Ds3_10 = Ds3(tt,:);

Dt_10 = Dt(tt,:);
x_10 = x(tt,:);

mve_10 = mve(tt,1);
cash_10 = cash(tt,1);
res_10 = res(tt,1);
frq_10 = frq(tt,1);
frqint_10 = frqint(tt,1);

v_10 = vhat1(tt,1);

% State level "random effect"
WSt_10 = Dst_10*Dst_10';

% SIC1 level "random effect"
WS1_10 = Ds1_10*Ds1_10';

% SIC2 level "random effect"
WS2_10 = Ds2_10*Ds2_10';

% SIC3 level "random effect"
WS3_10 = Ds3_10*Ds3_10';

% Distance matrices
mveDist_10 = zeros(N10);
cashDist_10 = zeros(N10);
frqintDist_10 = zeros(N10);
for ii = 1:N10
    mveDist_10(:,ii) = abs(mve_10 - mve_10(ii));
    cashDist_10(:,ii) = abs(cash_10 - cash_10(ii));
    frqintDist_10(:,ii) = sqrt((res_10 - res_10(ii)).^2 + (frq_10 - frq_10(ii)).^2 ...
        + (frqint_10 - frqint_10(ii)).^2);
end

% Characteristic "neighbors"
Wmve_10 = (1-.10*mveDist_10).*(.10*mveDist_10 < 1) - eye(N10);
Wcash_10 = (1-cashDist_10).*(cashDist_10 < 1) - eye(N10);
Wfrqint_10 = (1-.10*frqintDist_10).*(.10*frqintDist_10 < 1) - eye(N10);


%% Time period 11
tt = t_ind == 11;
N11 = sum(tt);

% Data for time period 11
Dst_11 = Dst(tt,:);
Ds1_11 = Ds1(tt,:);
Ds2_11 = Ds2(tt,:);
Ds3_11 = Ds3(tt,:);

Dt_11 = Dt(tt,:);
x_11 = x(tt,:);

mve_11 = mve(tt,1);
cash_11 = cash(tt,1);
res_11 = res(tt,1);
frq_11 = frq(tt,1);
frqint_11 = frqint(tt,1);

v_11 = vhat1(tt,1);

% State level "random effect"
WSt_11 = Dst_11*Dst_11';

% SIC1 level "random effect"
WS1_11 = Ds1_11*Ds1_11';

% SIC2 level "random effect"
WS2_11 = Ds2_11*Ds2_11';

% SIC3 level "random effect"
WS3_11 = Ds3_11*Ds3_11';

% Distance matrices
mveDist_11 = zeros(N11);
cashDist_11 = zeros(N11);
frqintDist_11 = zeros(N11);
for ii = 1:N11
    mveDist_11(:,ii) = abs(mve_11 - mve_11(ii));
    cashDist_11(:,ii) = abs(cash_11 - cash_11(ii));
    frqintDist_11(:,ii) = sqrt((res_11 - res_11(ii)).^2 + (frq_11 - frq_11(ii)).^2 ...
        + (frqint_11 - frqint_11(ii)).^2);
end

% Characteristic "neighbors"
Wmve_11 = (1-.11*mveDist_11).*(.11*mveDist_11 < 1) - eye(N11);
Wcash_11 = (1-cashDist_11).*(cashDist_11 < 1) - eye(N11);
Wfrqint_11 = (1-.11*frqintDist_11).*(.11*frqintDist_11 < 1) - eye(N11);



%% Time period 12
tt = t_ind == 12;
N12 = sum(tt);

% Data for time period 12
Dst_12 = Dst(tt,:);
Ds1_12 = Ds1(tt,:);
Ds2_12 = Ds2(tt,:);
Ds3_12 = Ds3(tt,:);

Dt_12 = Dt(tt,:);
x_12 = x(tt,:);

mve_12 = mve(tt,1);
cash_12 = cash(tt,1);
res_12 = res(tt,1);
frq_12 = frq(tt,1);
frqint_12 = frqint(tt,1);

v_12 = vhat1(tt,1);

% State level "random effect"
WSt_12 = Dst_12*Dst_12';

% SIC1 level "random effect"
WS1_12 = Ds1_12*Ds1_12';

% SIC2 level "random effect"
WS2_12 = Ds2_12*Ds2_12';

% SIC3 level "random effect"
WS3_12 = Ds3_12*Ds3_12';

% Distance matrices
mveDist_12 = zeros(N12);
cashDist_12 = zeros(N12);
frqintDist_12 = zeros(N12);
for ii = 1:N12
    mveDist_12(:,ii) = abs(mve_12 - mve_12(ii));
    cashDist_12(:,ii) = abs(cash_12 - cash_12(ii));
    frqintDist_12(:,ii) = sqrt((res_12 - res_12(ii)).^2 + (frq_12 - frq_12(ii)).^2 ...
        + (frqint_12 - frqint_12(ii)).^2);
end

% Characteristic "neighbors"
Wmve_12 = (1-.12*mveDist_12).*(.12*mveDist_12 < 1) - eye(N12);
Wcash_12 = (1-cashDist_12).*(cashDist_12 < 1) - eye(N12);
Wfrqint_12 = (1-.12*frqintDist_12).*(.12*frqintDist_12 < 1) - eye(N12);



%% Time period 13
tt = t_ind == 13;
N13 = sum(tt);

% Data for time period 13
Dst_13 = Dst(tt,:);
Ds1_13 = Ds1(tt,:);
Ds2_13 = Ds2(tt,:);
Ds3_13 = Ds3(tt,:);

Dt_13 = Dt(tt,:);
x_13 = x(tt,:);

mve_13 = mve(tt,1);
cash_13 = cash(tt,1);
res_13 = res(tt,1);
frq_13 = frq(tt,1);
frqint_13 = frqint(tt,1);

v_13 = vhat1(tt,1);

% State level "random effect"
WSt_13 = Dst_13*Dst_13';

% SIC1 level "random effect"
WS1_13 = Ds1_13*Ds1_13';

% SIC2 level "random effect"
WS2_13 = Ds2_13*Ds2_13';

% SIC3 level "random effect"
WS3_13 = Ds3_13*Ds3_13';

% Distance matrices
mveDist_13 = zeros(N13);
cashDist_13 = zeros(N13);
frqintDist_13 = zeros(N13);
for ii = 1:N13
    mveDist_13(:,ii) = abs(mve_13 - mve_13(ii));
    cashDist_13(:,ii) = abs(cash_13 - cash_13(ii));
    frqintDist_13(:,ii) = sqrt((res_13 - res_13(ii)).^2 + (frq_13 - frq_13(ii)).^2 ...
        + (frqint_13 - frqint_13(ii)).^2);
end

% Characteristic "neighbors"
Wmve_13 = (1-.13*mveDist_13).*(.13*mveDist_13 < 1) - eye(N13);
Wcash_13 = (1-cashDist_13).*(cashDist_13 < 1) - eye(N13);
Wfrqint_13 = (1-.13*frqintDist_13).*(.13*frqintDist_13 < 1) - eye(N13);



%% Time period 14
tt = t_ind == 14;
N14 = sum(tt);

% Data for time period 14
Dst_14 = Dst(tt,:);
Ds1_14 = Ds1(tt,:);
Ds2_14 = Ds2(tt,:);
Ds3_14 = Ds3(tt,:);

Dt_14 = Dt(tt,:);
x_14 = x(tt,:);

mve_14 = mve(tt,1);
cash_14 = cash(tt,1);
res_14 = res(tt,1);
frq_14 = frq(tt,1);
frqint_14 = frqint(tt,1);

v_14 = vhat1(tt,1);

% State level "random effect"
WSt_14 = Dst_14*Dst_14';

% SIC1 level "random effect"
WS1_14 = Ds1_14*Ds1_14';

% SIC2 level "random effect"
WS2_14 = Ds2_14*Ds2_14';

% SIC3 level "random effect"
WS3_14 = Ds3_14*Ds3_14';

% Distance matrices
mveDist_14 = zeros(N14);
cashDist_14 = zeros(N14);
frqintDist_14 = zeros(N14);
for ii = 1:N14
    mveDist_14(:,ii) = abs(mve_14 - mve_14(ii));
    cashDist_14(:,ii) = abs(cash_14 - cash_14(ii));
    frqintDist_14(:,ii) = sqrt((res_14 - res_14(ii)).^2 + (frq_14 - frq_14(ii)).^2 ...
        + (frqint_14 - frqint_14(ii)).^2);
end

% Characteristic "neighbors"
Wmve_14 = (1-.14*mveDist_14).*(.14*mveDist_14 < 1) - eye(N14);
Wcash_14 = (1-cashDist_14).*(cashDist_14 < 1) - eye(N14);
Wfrqint_14 = (1-.14*frqintDist_14).*(.14*frqintDist_14 < 1) - eye(N14);



%% Time period 15
tt = t_ind == 15;
N15 = sum(tt);

% Data for time period 15
Dst_15 = Dst(tt,:);
Ds1_15 = Ds1(tt,:);
Ds2_15 = Ds2(tt,:);
Ds3_15 = Ds3(tt,:);

Dt_15 = Dt(tt,:);
x_15 = x(tt,:);

mve_15 = mve(tt,1);
cash_15 = cash(tt,1);
res_15 = res(tt,1);
frq_15 = frq(tt,1);
frqint_15 = frqint(tt,1);

v_15 = vhat1(tt,1);

% State level "random effect"
WSt_15 = Dst_15*Dst_15';

% SIC1 level "random effect"
WS1_15 = Ds1_15*Ds1_15';

% SIC2 level "random effect"
WS2_15 = Ds2_15*Ds2_15';

% SIC3 level "random effect"
WS3_15 = Ds3_15*Ds3_15';

% Distance matrices
mveDist_15 = zeros(N15);
cashDist_15 = zeros(N15);
frqintDist_15 = zeros(N15);
for ii = 1:N15
    mveDist_15(:,ii) = abs(mve_15 - mve_15(ii));
    cashDist_15(:,ii) = abs(cash_15 - cash_15(ii));
    frqintDist_15(:,ii) = sqrt((res_15 - res_15(ii)).^2 + (frq_15 - frq_15(ii)).^2 ...
        + (frqint_15 - frqint_15(ii)).^2);
end

% Characteristic "neighbors"
Wmve_15 = (1-.15*mveDist_15).*(.15*mveDist_15 < 1) - eye(N15);
Wcash_15 = (1-cashDist_15).*(cashDist_15 < 1) - eye(N15);
Wfrqint_15 = (1-.15*frqintDist_15).*(.15*frqintDist_15 < 1) - eye(N15);



%% Time period 16
tt = t_ind == 16;
N16 = sum(tt);

% Data for time period 16
Dst_16 = Dst(tt,:);
Ds1_16 = Ds1(tt,:);
Ds2_16 = Ds2(tt,:);
Ds3_16 = Ds3(tt,:);

Dt_16 = Dt(tt,:);
x_16 = x(tt,:);

mve_16 = mve(tt,1);
cash_16 = cash(tt,1);
res_16 = res(tt,1);
frq_16 = frq(tt,1);
frqint_16 = frqint(tt,1);

v_16 = vhat1(tt,1);

% State level "random effect"
WSt_16 = Dst_16*Dst_16';

% SIC1 level "random effect"
WS1_16 = Ds1_16*Ds1_16';

% SIC2 level "random effect"
WS2_16 = Ds2_16*Ds2_16';

% SIC3 level "random effect"
WS3_16 = Ds3_16*Ds3_16';

% Distance matrices
mveDist_16 = zeros(N16);
cashDist_16 = zeros(N16);
frqintDist_16 = zeros(N16);
for ii = 1:N16
    mveDist_16(:,ii) = abs(mve_16 - mve_16(ii));
    cashDist_16(:,ii) = abs(cash_16 - cash_16(ii));
    frqintDist_16(:,ii) = sqrt((res_16 - res_16(ii)).^2 + (frq_16 - frq_16(ii)).^2 ...
        + (frqint_16 - frqint_16(ii)).^2);
end

% Characteristic "neighbors"
Wmve_16 = (1-.16*mveDist_16).*(.16*mveDist_16 < 1) - eye(N16);
Wcash_16 = (1-cashDist_16).*(cashDist_16 < 1) - eye(N16);
Wfrqint_16 = (1-.16*frqintDist_16).*(.16*frqintDist_16 < 1) - eye(N16);



%% Time period 17
tt = t_ind == 17;
N17 = sum(tt);

% Data for time period 17
Dst_17 = Dst(tt,:);
Ds1_17 = Ds1(tt,:);
Ds2_17 = Ds2(tt,:);
Ds3_17 = Ds3(tt,:);

Dt_17 = Dt(tt,:);
x_17 = x(tt,:);

mve_17 = mve(tt,1);
cash_17 = cash(tt,1);
res_17 = res(tt,1);
frq_17 = frq(tt,1);
frqint_17 = frqint(tt,1);

v_17 = vhat1(tt,1);

% State level "random effect"
WSt_17 = Dst_17*Dst_17';

% SIC1 level "random effect"
WS1_17 = Ds1_17*Ds1_17';

% SIC2 level "random effect"
WS2_17 = Ds2_17*Ds2_17';

% SIC3 level "random effect"
WS3_17 = Ds3_17*Ds3_17';

% Distance matrices
mveDist_17 = zeros(N17);
cashDist_17 = zeros(N17);
frqintDist_17 = zeros(N17);
for ii = 1:N17
    mveDist_17(:,ii) = abs(mve_17 - mve_17(ii));
    cashDist_17(:,ii) = abs(cash_17 - cash_17(ii));
    frqintDist_17(:,ii) = sqrt((res_17 - res_17(ii)).^2 + (frq_17 - frq_17(ii)).^2 ...
        + (frqint_17 - frqint_17(ii)).^2);
end

% Characteristic "neighbors"
Wmve_17 = (1-.17*mveDist_17).*(.17*mveDist_17 < 1) - eye(N17);
Wcash_17 = (1-cashDist_17).*(cashDist_17 < 1) - eye(N17);
Wfrqint_17 = (1-.17*frqintDist_17).*(.17*frqintDist_17 < 1) - eye(N17);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define functions for covariance matrices
% "Model":
% v_t = (a1*Wmve_t + a2*Wcash_t + a3*Wfrqint_t)*v_t + u_t
%     = sum(ak*Wk_t)*v_t + u_t
% =>
% (I_Nt - sum(ak*Wk_t))*v_t = u_t
% => 
% v_t = (I_Nt - sum(ak*Wk_t))^(-1) u_t = Omega_t u_t
%
% u_{it} = delta_{state(i,t)} + delta_{sic1(i,t)} + delta_{sic2(i,t)}
%              + delta_{sic3(i,t)} + sigma_{x_{it}'phi} eta_{it}
% eta_{it} iid (0,1) r.v.
%
% u_t = DSt_t*delta_state + Dsic1_t*delta_sic1 + Dsic2_t*delta_sic2
%              + Dsic3_t*delta_sic3 + diag(sigma_t)*eta_t
% E[u_t u_t'] = s2_state*WSt_t + s2_sic1*WS1_t + s2_sic2*WS2_t
%              + s2_sic3*WS3_t + diag(sigma_t^2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OmegaFun = @(coef,W1,W2,W3,Nt) inv(eye(Nt) - coef(1)*W1 - coef(2)*W2 - coef(3)*W3);
REFun = @(coef,W1,W2,W3,W4) coef(1)*W1 + coef(2)*W2 + coef(3)*W3 + coef(4)*W4;
VarFun = @(coef,Z) exp(Z*coef);  
        % Want this to depend on the same nine variables as the mean,
        % time effects, and 1 digit industry maybe?        
Sigma_t = @(coef,WO1,WO2,WO3,Nt,WR1,WR2,WR3,WR4,Z) ...
    OmegaFun(coef(1:3),WO1,WO2,WO3,Nt)*(REFun(coef(4:7),WR1,WR2,WR3,WR4) ...
    + diag(VarFun(coef(8:end),Z)))*OmegaFun(coef(1:3),WO1,WO2,WO3,Nt)';
        
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Need to complete the likelihood by tying everything together.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z_1 = [Dt_1 x_1 Ds1_1(:,2:end)];
z_2 = [Dt_2 x_2 Ds1_2(:,2:end)];
z_3 = [Dt_3 x_3 Ds1_3(:,2:end)];
z_4 = [Dt_4 x_4 Ds1_4(:,2:end)];
z_5 = [Dt_5 x_5 Ds1_5(:,2:end)];
z_6 = [Dt_6 x_6 Ds1_6(:,2:end)];
z_7 = [Dt_7 x_7 Ds1_7(:,2:end)];
z_8 = [Dt_8 x_8 Ds1_8(:,2:end)];
z_9 = [Dt_9 x_9 Ds1_9(:,2:end)];
z_10 = [Dt_10 x_10 Ds1_10(:,2:end)];
z_11 = [Dt_11 x_11 Ds1_11(:,2:end)];
z_12 = [Dt_12 x_12 Ds1_12(:,2:end)];
z_13 = [Dt_13 x_13 Ds1_13(:,2:end)];
z_14 = [Dt_14 x_14 Ds1_14(:,2:end)];
z_15 = [Dt_15 x_15 Ds1_15(:,2:end)];
z_16 = [Dt_16 x_16 Ds1_16(:,2:end)];
z_17 = [Dt_17 x_17 Ds1_17(:,2:end)];

pvar = size(z_1,2);

% log likelihood contribution for time period t
ltt = @(coef,WO1,WO2,WO3,Nt,WR1,WR2,WR3,WR4,Z,V) ...
    -(.5*Nt)*log(2*pi)+(.5*Nt)*log(var(V))...
    -.5*log(det(Sigma_t(coef,WO1,WO2,WO3,Nt,WR1,WR2,WR3,WR4,Z)/var(V))) ...
    -.5*V'*(Sigma_t(coef,WO1,WO2,WO3,Nt,WR1,WR2,WR3,WR4,Z)\V);

% log likelihood
ll = @(coef,WO1_1, WO2_1 ,WO3_1 ,Nt_1 ,WR1_1 ,WR2_1 ,WR3_1 ,WR4_1 ,Z_1 ,V_1,...
            WO1_2, WO2_2 ,WO3_2 ,Nt_2 ,WR1_2 ,WR2_2 ,WR3_2 ,WR4_2 ,Z_2 ,V_2,...
            WO1_3, WO2_3 ,WO3_3 ,Nt_3 ,WR1_3 ,WR2_3 ,WR3_3 ,WR4_3 ,Z_3 ,V_3,...
            WO1_4, WO2_4 ,WO3_4 ,Nt_4 ,WR1_4 ,WR2_4 ,WR3_4 ,WR4_4 ,Z_4 ,V_4,...
            WO1_5, WO2_5 ,WO3_5 ,Nt_5 ,WR1_5 ,WR2_5 ,WR3_5 ,WR4_5 ,Z_5 ,V_5,...
            WO1_6, WO2_6 ,WO3_6 ,Nt_6 ,WR1_6 ,WR2_6 ,WR3_6 ,WR4_6 ,Z_6 ,V_6,...
            WO1_7, WO2_7 ,WO3_7 ,Nt_7 ,WR1_7 ,WR2_7 ,WR3_7 ,WR4_7 ,Z_7 ,V_7,...
            WO1_8, WO2_8 ,WO3_8 ,Nt_8 ,WR1_8 ,WR2_8 ,WR3_8 ,WR4_8 ,Z_8 ,V_8,...
            WO1_9, WO2_9 ,WO3_9 ,Nt_9 ,WR1_9 ,WR2_9 ,WR3_9 ,WR4_9 ,Z_9 ,V_9,...
            WO1_10,WO2_10,WO3_10,Nt_10,WR1_10,WR2_10,WR3_10,WR4_10,Z_10,V_10,...
            WO1_11,WO2_11,WO3_11,Nt_11,WR1_11,WR2_11,WR3_11,WR4_11,Z_11,V_11,...
            WO1_12,WO2_12,WO3_12,Nt_12,WR1_12,WR2_12,WR3_12,WR4_12,Z_12,V_12,...
            WO1_13,WO2_13,WO3_13,Nt_13,WR1_13,WR2_13,WR3_13,WR4_13,Z_13,V_13,...
            WO1_14,WO2_14,WO3_14,Nt_14,WR1_14,WR2_14,WR3_14,WR4_14,Z_14,V_14,...
            WO1_15,WO2_15,WO3_15,Nt_15,WR1_15,WR2_15,WR3_15,WR4_15,Z_15,V_15,...
            WO1_16,WO2_16,WO3_16,Nt_16,WR1_16,WR2_16,WR3_16,WR4_16,Z_16,V_16,...
            WO1_17,WO2_17,WO3_17,Nt_17,WR1_17,WR2_17,WR3_17,WR4_17,Z_17,V_17) ...
    ltt(coef,WO1_1, WO2_1 ,WO3_1 ,Nt_1 ,WR1_1 ,WR2_1 ,WR3_1 ,WR4_1 ,Z_1 ,V_1) ...
    + ltt(coef,WO1_2, WO2_2 ,WO3_2 ,Nt_2 ,WR1_2 ,WR2_2 ,WR3_2 ,WR4_2 ,Z_2 ,V_2) ...
    + ltt(coef,WO1_3, WO2_3 ,WO3_3 ,Nt_3 ,WR1_3 ,WR2_3 ,WR3_3 ,WR4_3 ,Z_3 ,V_3) ...
    + ltt(coef,WO1_4, WO2_4 ,WO3_4 ,Nt_4 ,WR1_4 ,WR2_4 ,WR3_4 ,WR4_4 ,Z_4 ,V_4) ...
    + ltt(coef,WO1_5, WO2_5 ,WO3_5 ,Nt_5 ,WR1_5 ,WR2_5 ,WR3_5 ,WR4_5 ,Z_5 ,V_5) ...
    + ltt(coef,WO1_6, WO2_6 ,WO3_6 ,Nt_6 ,WR1_6 ,WR2_6 ,WR3_6 ,WR4_6 ,Z_6 ,V_6) ...
    + ltt(coef,WO1_7, WO2_7 ,WO3_7 ,Nt_7 ,WR1_7 ,WR2_7 ,WR3_7 ,WR4_7 ,Z_7 ,V_7) ...
    + ltt(coef,WO1_8, WO2_8 ,WO3_8 ,Nt_8 ,WR1_8 ,WR2_8 ,WR3_8 ,WR4_8 ,Z_8 ,V_8) ...
    + ltt(coef,WO1_9, WO2_9 ,WO3_9 ,Nt_9 ,WR1_9 ,WR2_9 ,WR3_9 ,WR4_9 ,Z_9 ,V_9) ...
    + ltt(coef,WO1_10,WO2_10,WO3_10,Nt_10,WR1_10,WR2_10,WR3_10,WR4_10,Z_10,V_10) ...
    + ltt(coef,WO1_11,WO2_11,WO3_11,Nt_11,WR1_11,WR2_11,WR3_11,WR4_11,Z_11,V_11) ...
    + ltt(coef,WO1_12,WO2_12,WO3_12,Nt_12,WR1_12,WR2_12,WR3_12,WR4_12,Z_12,V_12) ...
    + ltt(coef,WO1_13,WO2_13,WO3_13,Nt_13,WR1_13,WR2_13,WR3_13,WR4_13,Z_13,V_13) ...
    + ltt(coef,WO1_14,WO2_14,WO3_14,Nt_14,WR1_14,WR2_14,WR3_14,WR4_14,Z_14,V_14) ...
    + ltt(coef,WO1_15,WO2_15,WO3_15,Nt_15,WR1_15,WR2_15,WR3_15,WR4_15,Z_15,V_15) ...
    + ltt(coef,WO1_16,WO2_16,WO3_16,Nt_16,WR1_16,WR2_16,WR3_16,WR4_16,Z_16,V_16) ...
    + ltt(coef,WO1_17,WO2_17,WO3_17,Nt_17,WR1_17,WR2_17,WR3_17,WR4_17,Z_17,V_17);
    
% %% Try to maximize the likelihood
% 
% init = [.01;.01;.01;.001;.001;.001;.001;...
%     log(var(v_1));log(var(v_2));log(var(v_3));log(var(v_4));log(var(v_5));...
%     log(var(v_6));log(var(v_7));log(var(v_8));log(var(v_9));log(var(v_10));...
%     log(var(v_11));log(var(v_12));log(var(v_13));log(var(v_14));log(var(v_15));...
%     log(var(v_16));log(var(v_17));
%     zeros(pvar-17,1)];
% 
% cvp_est = fminsearch(@(theta) -ll(theta,Wmve_1, Wcash_1 ,Wfrqint_1 ,N1 ,WSt_1 ,WS1_1 ,WS2_1 ,WS3_1 ,z_1 ,v_1,...
%             Wmve_2, Wcash_2 ,Wfrqint_2 ,N2 ,WSt_2 ,WS1_2 ,WS2_2 ,WS3_2 ,z_2 ,v_2,...
%             Wmve_3, Wcash_3 ,Wfrqint_3 ,N3 ,WSt_3 ,WS1_3 ,WS2_3 ,WS3_3 ,z_3 ,v_3,...
%             Wmve_4, Wcash_4 ,Wfrqint_4 ,N4 ,WSt_4 ,WS1_4 ,WS2_4 ,WS3_4 ,z_4 ,v_4,...
%             Wmve_5, Wcash_5 ,Wfrqint_5 ,N5 ,WSt_5 ,WS1_5 ,WS2_5 ,WS3_5 ,z_5 ,v_5,...
%             Wmve_6, Wcash_6 ,Wfrqint_6 ,N6 ,WSt_6 ,WS1_6 ,WS2_6 ,WS3_6 ,z_6 ,v_6,...
%             Wmve_7, Wcash_7 ,Wfrqint_7 ,N7 ,WSt_7 ,WS1_7 ,WS2_7 ,WS3_7 ,z_7 ,v_7,...
%             Wmve_8, Wcash_8 ,Wfrqint_8 ,N8 ,WSt_8 ,WS1_8 ,WS2_8 ,WS3_8 ,z_8 ,v_8,...
%             Wmve_9, Wcash_9 ,Wfrqint_9 ,N9 ,WSt_9 ,WS1_9 ,WS2_9 ,WS3_9 ,z_9 ,v_9,...
%             Wmve_10,Wcash_10,Wfrqint_10,N10,WSt_10,WS1_10,WS2_10,WS3_10,z_10,v_10,...
%             Wmve_11,Wcash_11,Wfrqint_11,N11,WSt_11,WS1_11,WS2_11,WS3_11,z_11,v_11,...
%             Wmve_12,Wcash_12,Wfrqint_12,N12,WSt_12,WS1_12,WS2_12,WS3_12,z_12,v_12,...
%             Wmve_13,Wcash_13,Wfrqint_13,N13,WSt_13,WS1_13,WS2_13,WS3_13,z_13,v_13,...
%             Wmve_14,Wcash_14,Wfrqint_14,N14,WSt_14,WS1_14,WS2_14,WS3_14,z_14,v_14,...
%             Wmve_15,Wcash_15,Wfrqint_15,N15,WSt_15,WS1_15,WS2_15,WS3_15,z_15,v_15,...
%             Wmve_16,Wcash_16,Wfrqint_16,N16,WSt_16,WS1_16,WS2_16,WS3_16,z_16,v_16,...
%             Wmve_17,Wcash_17,Wfrqint_17,N17,WSt_17,WS1_17,WS2_17,WS3_17,z_17,v_17),...
%             init,optimset('disp','iter'));
% 
% % Note that estimation doesn't converge here.  Save current values.
% save bcv_optim1 cvp_est ;
% 
% %% Load in estimation results from above.  Comment either this block or the block above out.
% load bcv_optim1 ;
%         
% %% See whether these estimates make any sense        
% Sigma_1 = Sigma_t(cvp_est,Wmve_1, Wcash_1 ,Wfrqint_1 ,N1 ,WSt_1 ,WS1_1 ,WS2_1 ,WS3_1 ,z_1);
% Sigma_2 = Sigma_t(cvp_est,Wmve_2, Wcash_2 ,Wfrqint_2 ,N2 ,WSt_2 ,WS1_2 ,WS2_2 ,WS3_2 ,z_2);
% Sigma_3 = Sigma_t(cvp_est,Wmve_3, Wcash_3 ,Wfrqint_3 ,N3 ,WSt_3 ,WS1_3 ,WS2_3 ,WS3_3 ,z_3);
% Sigma_4 = Sigma_t(cvp_est,Wmve_4, Wcash_4 ,Wfrqint_4 ,N4 ,WSt_4 ,WS1_4 ,WS2_4 ,WS3_4 ,z_4);
% Sigma_5 = Sigma_t(cvp_est,Wmve_5, Wcash_5 ,Wfrqint_5 ,N5 ,WSt_5 ,WS1_5 ,WS2_5 ,WS3_5 ,z_5);
% Sigma_6 = Sigma_t(cvp_est,Wmve_6, Wcash_6 ,Wfrqint_6 ,N6 ,WSt_6 ,WS1_6 ,WS2_6 ,WS3_6 ,z_6);
% Sigma_7 = Sigma_t(cvp_est,Wmve_7, Wcash_7 ,Wfrqint_7 ,N7 ,WSt_7 ,WS1_7 ,WS2_7 ,WS3_7 ,z_7);
% Sigma_8 = Sigma_t(cvp_est,Wmve_8, Wcash_8 ,Wfrqint_8 ,N8 ,WSt_8 ,WS1_8 ,WS2_8 ,WS3_8 ,z_8);
% Sigma_9 = Sigma_t(cvp_est,Wmve_9, Wcash_9 ,Wfrqint_9 ,N9 ,WSt_9 ,WS1_9 ,WS2_9 ,WS3_9 ,z_9);
% Sigma_10 = Sigma_t(cvp_est,Wmve_10, Wcash_10 ,Wfrqint_10 ,N10 ,WSt_10 ,WS1_10 ,WS2_10 ,WS3_10 ,z_10);
% Sigma_11 = Sigma_t(cvp_est,Wmve_11, Wcash_11 ,Wfrqint_11 ,N11 ,WSt_11 ,WS1_11 ,WS2_11 ,WS3_11 ,z_11);
% Sigma_12 = Sigma_t(cvp_est,Wmve_12, Wcash_12 ,Wfrqint_12 ,N12 ,WSt_12 ,WS1_12 ,WS2_12 ,WS3_12 ,z_12);
% Sigma_13 = Sigma_t(cvp_est,Wmve_13, Wcash_13 ,Wfrqint_13 ,N13 ,WSt_13 ,WS1_13 ,WS2_13 ,WS3_13 ,z_13);
% Sigma_14 = Sigma_t(cvp_est,Wmve_14, Wcash_14 ,Wfrqint_14 ,N14 ,WSt_14 ,WS1_14 ,WS2_14 ,WS3_14 ,z_14);
% Sigma_15 = Sigma_t(cvp_est,Wmve_15, Wcash_15 ,Wfrqint_15 ,N15 ,WSt_15 ,WS1_15 ,WS2_15 ,WS3_15 ,z_15);
% Sigma_16 = Sigma_t(cvp_est,Wmve_16, Wcash_16 ,Wfrqint_16 ,N16 ,WSt_16 ,WS1_16 ,WS2_16 ,WS3_16 ,z_16);
% Sigma_17 = Sigma_t(cvp_est,Wmve_17, Wcash_17 ,Wfrqint_17 ,N17 ,WSt_17 ,WS1_17 ,WS2_17 ,WS3_17 ,z_17);
% 
% % p.d.?
% disp(min(eig(Sigma_1)))
% disp(min(eig(Sigma_2)))
% disp(min(eig(Sigma_3)))
% disp(min(eig(Sigma_4)))
% disp(min(eig(Sigma_5)))
% disp(min(eig(Sigma_6)))
% disp(min(eig(Sigma_7)))
% disp(min(eig(Sigma_8)))
% disp(min(eig(Sigma_9)))
% disp(min(eig(Sigma_10)))
% disp(min(eig(Sigma_11)))
% disp(min(eig(Sigma_12)))
% disp(min(eig(Sigma_13)))
% disp(min(eig(Sigma_14)))
% disp(min(eig(Sigma_15)))
% disp(min(eig(Sigma_16)))
% disp(min(eig(Sigma_17)))
% 
% % Variances sort of reasonable?
% disp(min(diag(Sigma_1)))
% disp(min(diag(Sigma_2)))
% disp(min(diag(Sigma_3)))
% disp(min(diag(Sigma_4)))
% disp(min(diag(Sigma_5)))
% disp(min(diag(Sigma_6)))
% disp(min(diag(Sigma_7)))
% disp(min(diag(Sigma_8)))
% disp(min(diag(Sigma_9)))
% disp(min(diag(Sigma_10)))
% disp(min(diag(Sigma_11)))
% disp(min(diag(Sigma_12)))
% disp(min(diag(Sigma_13)))
% disp(min(diag(Sigma_14)))
% disp(min(diag(Sigma_15)))
% disp(min(diag(Sigma_16)))
% disp(min(diag(Sigma_17)))
% 
% % Any big correlations?
% disp(max(max(abs(diag(sqrt(1./diag(Sigma_1)))*Sigma_1*diag(sqrt(1./diag(Sigma_1))) - eye(N1)))))
% disp(max(max(abs(diag(sqrt(1./diag(Sigma_2)))*Sigma_2*diag(sqrt(1./diag(Sigma_2))) - eye(N2)))))
% disp(max(max(abs(diag(sqrt(1./diag(Sigma_3)))*Sigma_3*diag(sqrt(1./diag(Sigma_3))) - eye(N3)))))
% disp(max(max(abs(diag(sqrt(1./diag(Sigma_4)))*Sigma_4*diag(sqrt(1./diag(Sigma_4))) - eye(N4)))))
% disp(max(max(abs(diag(sqrt(1./diag(Sigma_5)))*Sigma_5*diag(sqrt(1./diag(Sigma_5))) - eye(N5)))))
% disp(max(max(abs(diag(sqrt(1./diag(Sigma_6)))*Sigma_6*diag(sqrt(1./diag(Sigma_6))) - eye(N6)))))
% disp(max(max(abs(diag(sqrt(1./diag(Sigma_7)))*Sigma_7*diag(sqrt(1./diag(Sigma_7))) - eye(N7)))))
% disp(max(max(abs(diag(sqrt(1./diag(Sigma_8)))*Sigma_8*diag(sqrt(1./diag(Sigma_8))) - eye(N8)))))
% disp(max(max(abs(diag(sqrt(1./diag(Sigma_9)))*Sigma_9*diag(sqrt(1./diag(Sigma_9))) - eye(N9)))))
% disp(max(max(abs(diag(sqrt(1./diag(Sigma_10)))*Sigma_10*diag(sqrt(1./diag(Sigma_10))) - eye(N10)))))
% disp(max(max(abs(diag(sqrt(1./diag(Sigma_11)))*Sigma_11*diag(sqrt(1./diag(Sigma_11))) - eye(N11)))))
% disp(max(max(abs(diag(sqrt(1./diag(Sigma_12)))*Sigma_12*diag(sqrt(1./diag(Sigma_12))) - eye(N12)))))
% disp(max(max(abs(diag(sqrt(1./diag(Sigma_13)))*Sigma_13*diag(sqrt(1./diag(Sigma_13))) - eye(N13)))))
% disp(max(max(abs(diag(sqrt(1./diag(Sigma_14)))*Sigma_14*diag(sqrt(1./diag(Sigma_14))) - eye(N14)))))
% disp(max(max(abs(diag(sqrt(1./diag(Sigma_15)))*Sigma_15*diag(sqrt(1./diag(Sigma_15))) - eye(N15)))))
% disp(max(max(abs(diag(sqrt(1./diag(Sigma_16)))*Sigma_16*diag(sqrt(1./diag(Sigma_16))) - eye(N16)))))
% disp(max(max(abs(diag(sqrt(1./diag(Sigma_17)))*Sigma_17*diag(sqrt(1./diag(Sigma_17))) - eye(N17)))))
% 
% %% Try another optimization step
% 
% [cvp_est2,~,~,~,~,hess2] = fminunc(@(theta) -ll(theta,Wmve_1, Wcash_1 ,Wfrqint_1 ,N1 ,WSt_1 ,WS1_1 ,WS2_1 ,WS3_1 ,z_1 ,v_1,...
%             Wmve_2, Wcash_2 ,Wfrqint_2 ,N2 ,WSt_2 ,WS1_2 ,WS2_2 ,WS3_2 ,z_2 ,v_2,...
%             Wmve_3, Wcash_3 ,Wfrqint_3 ,N3 ,WSt_3 ,WS1_3 ,WS2_3 ,WS3_3 ,z_3 ,v_3,...
%             Wmve_4, Wcash_4 ,Wfrqint_4 ,N4 ,WSt_4 ,WS1_4 ,WS2_4 ,WS3_4 ,z_4 ,v_4,...
%             Wmve_5, Wcash_5 ,Wfrqint_5 ,N5 ,WSt_5 ,WS1_5 ,WS2_5 ,WS3_5 ,z_5 ,v_5,...
%             Wmve_6, Wcash_6 ,Wfrqint_6 ,N6 ,WSt_6 ,WS1_6 ,WS2_6 ,WS3_6 ,z_6 ,v_6,...
%             Wmve_7, Wcash_7 ,Wfrqint_7 ,N7 ,WSt_7 ,WS1_7 ,WS2_7 ,WS3_7 ,z_7 ,v_7,...
%             Wmve_8, Wcash_8 ,Wfrqint_8 ,N8 ,WSt_8 ,WS1_8 ,WS2_8 ,WS3_8 ,z_8 ,v_8,...
%             Wmve_9, Wcash_9 ,Wfrqint_9 ,N9 ,WSt_9 ,WS1_9 ,WS2_9 ,WS3_9 ,z_9 ,v_9,...
%             Wmve_10,Wcash_10,Wfrqint_10,N10,WSt_10,WS1_10,WS2_10,WS3_10,z_10,v_10,...
%             Wmve_11,Wcash_11,Wfrqint_11,N11,WSt_11,WS1_11,WS2_11,WS3_11,z_11,v_11,...
%             Wmve_12,Wcash_12,Wfrqint_12,N12,WSt_12,WS1_12,WS2_12,WS3_12,z_12,v_12,...
%             Wmve_13,Wcash_13,Wfrqint_13,N13,WSt_13,WS1_13,WS2_13,WS3_13,z_13,v_13,...
%             Wmve_14,Wcash_14,Wfrqint_14,N14,WSt_14,WS1_14,WS2_14,WS3_14,z_14,v_14,...
%             Wmve_15,Wcash_15,Wfrqint_15,N15,WSt_15,WS1_15,WS2_15,WS3_15,z_15,v_15,...
%             Wmve_16,Wcash_16,Wfrqint_16,N16,WSt_16,WS1_16,WS2_16,WS3_16,z_16,v_16,...
%             Wmve_17,Wcash_17,Wfrqint_17,N17,WSt_17,WS1_17,WS2_17,WS3_17,z_17,v_17),...
%             cvp_est,optimset('disp','iter'));
%         
% save bcv_optim2 cvp_est2 hess2 ;      
% 
% %% Try optimization step 3
% 
% [cvp_est3,~,~,~,~,hess3] = fminunc(@(theta) -ll(theta,Wmve_1, Wcash_1 ,Wfrqint_1 ,N1 ,WSt_1 ,WS1_1 ,WS2_1 ,WS3_1 ,z_1 ,v_1,...
%             Wmve_2, Wcash_2 ,Wfrqint_2 ,N2 ,WSt_2 ,WS1_2 ,WS2_2 ,WS3_2 ,z_2 ,v_2,...
%             Wmve_3, Wcash_3 ,Wfrqint_3 ,N3 ,WSt_3 ,WS1_3 ,WS2_3 ,WS3_3 ,z_3 ,v_3,...
%             Wmve_4, Wcash_4 ,Wfrqint_4 ,N4 ,WSt_4 ,WS1_4 ,WS2_4 ,WS3_4 ,z_4 ,v_4,...
%             Wmve_5, Wcash_5 ,Wfrqint_5 ,N5 ,WSt_5 ,WS1_5 ,WS2_5 ,WS3_5 ,z_5 ,v_5,...
%             Wmve_6, Wcash_6 ,Wfrqint_6 ,N6 ,WSt_6 ,WS1_6 ,WS2_6 ,WS3_6 ,z_6 ,v_6,...
%             Wmve_7, Wcash_7 ,Wfrqint_7 ,N7 ,WSt_7 ,WS1_7 ,WS2_7 ,WS3_7 ,z_7 ,v_7,...
%             Wmve_8, Wcash_8 ,Wfrqint_8 ,N8 ,WSt_8 ,WS1_8 ,WS2_8 ,WS3_8 ,z_8 ,v_8,...
%             Wmve_9, Wcash_9 ,Wfrqint_9 ,N9 ,WSt_9 ,WS1_9 ,WS2_9 ,WS3_9 ,z_9 ,v_9,...
%             Wmve_10,Wcash_10,Wfrqint_10,N10,WSt_10,WS1_10,WS2_10,WS3_10,z_10,v_10,...
%             Wmve_11,Wcash_11,Wfrqint_11,N11,WSt_11,WS1_11,WS2_11,WS3_11,z_11,v_11,...
%             Wmve_12,Wcash_12,Wfrqint_12,N12,WSt_12,WS1_12,WS2_12,WS3_12,z_12,v_12,...
%             Wmve_13,Wcash_13,Wfrqint_13,N13,WSt_13,WS1_13,WS2_13,WS3_13,z_13,v_13,...
%             Wmve_14,Wcash_14,Wfrqint_14,N14,WSt_14,WS1_14,WS2_14,WS3_14,z_14,v_14,...
%             Wmve_15,Wcash_15,Wfrqint_15,N15,WSt_15,WS1_15,WS2_15,WS3_15,z_15,v_15,...
%             Wmve_16,Wcash_16,Wfrqint_16,N16,WSt_16,WS1_16,WS2_16,WS3_16,z_16,v_16,...
%             Wmve_17,Wcash_17,Wfrqint_17,N17,WSt_17,WS1_17,WS2_17,WS3_17,z_17,v_17),...
%             cvp_est2,optimset('disp','iter'));
%         
% save bcv_optim3 cvp_est3 hess3 ;      
% 
% %% Try optimization step 4
% 
% [cvp_est4,~,~,~,~,hess4] = fminunc(@(theta) -ll(theta,Wmve_1, Wcash_1 ,Wfrqint_1 ,N1 ,WSt_1 ,WS1_1 ,WS2_1 ,WS3_1 ,z_1 ,v_1,...
%             Wmve_2, Wcash_2 ,Wfrqint_2 ,N2 ,WSt_2 ,WS1_2 ,WS2_2 ,WS3_2 ,z_2 ,v_2,...
%             Wmve_3, Wcash_3 ,Wfrqint_3 ,N3 ,WSt_3 ,WS1_3 ,WS2_3 ,WS3_3 ,z_3 ,v_3,...
%             Wmve_4, Wcash_4 ,Wfrqint_4 ,N4 ,WSt_4 ,WS1_4 ,WS2_4 ,WS3_4 ,z_4 ,v_4,...
%             Wmve_5, Wcash_5 ,Wfrqint_5 ,N5 ,WSt_5 ,WS1_5 ,WS2_5 ,WS3_5 ,z_5 ,v_5,...
%             Wmve_6, Wcash_6 ,Wfrqint_6 ,N6 ,WSt_6 ,WS1_6 ,WS2_6 ,WS3_6 ,z_6 ,v_6,...
%             Wmve_7, Wcash_7 ,Wfrqint_7 ,N7 ,WSt_7 ,WS1_7 ,WS2_7 ,WS3_7 ,z_7 ,v_7,...
%             Wmve_8, Wcash_8 ,Wfrqint_8 ,N8 ,WSt_8 ,WS1_8 ,WS2_8 ,WS3_8 ,z_8 ,v_8,...
%             Wmve_9, Wcash_9 ,Wfrqint_9 ,N9 ,WSt_9 ,WS1_9 ,WS2_9 ,WS3_9 ,z_9 ,v_9,...
%             Wmve_10,Wcash_10,Wfrqint_10,N10,WSt_10,WS1_10,WS2_10,WS3_10,z_10,v_10,...
%             Wmve_11,Wcash_11,Wfrqint_11,N11,WSt_11,WS1_11,WS2_11,WS3_11,z_11,v_11,...
%             Wmve_12,Wcash_12,Wfrqint_12,N12,WSt_12,WS1_12,WS2_12,WS3_12,z_12,v_12,...
%             Wmve_13,Wcash_13,Wfrqint_13,N13,WSt_13,WS1_13,WS2_13,WS3_13,z_13,v_13,...
%             Wmve_14,Wcash_14,Wfrqint_14,N14,WSt_14,WS1_14,WS2_14,WS3_14,z_14,v_14,...
%             Wmve_15,Wcash_15,Wfrqint_15,N15,WSt_15,WS1_15,WS2_15,WS3_15,z_15,v_15,...
%             Wmve_16,Wcash_16,Wfrqint_16,N16,WSt_16,WS1_16,WS2_16,WS3_16,z_16,v_16,...
%             Wmve_17,Wcash_17,Wfrqint_17,N17,WSt_17,WS1_17,WS2_17,WS3_17,z_17,v_17),...
%             cvp_est3,optimset('disp','iter','MaxFunEvals',8000));
%         
% save bcv_optim3 cvp_est4 hess4 ;      % Oops.  Wrote over previous :(
% save bcv_optim4 cvp_est4 hess4 ;
% save bcv_optim3 cvp_est3 hess3 ;
% 
% %% Load in estimation results from above.  Comment either this block or the block above out.
% load bcv_optim4 ;
% 
% %% See whether these estimates make any sense        
% Sigma_1 = Sigma_t(cvp_est4,Wmve_1, Wcash_1 ,Wfrqint_1 ,N1 ,WSt_1 ,WS1_1 ,WS2_1 ,WS3_1 ,z_1);
% Sigma_2 = Sigma_t(cvp_est4,Wmve_2, Wcash_2 ,Wfrqint_2 ,N2 ,WSt_2 ,WS1_2 ,WS2_2 ,WS3_2 ,z_2);
% Sigma_3 = Sigma_t(cvp_est4,Wmve_3, Wcash_3 ,Wfrqint_3 ,N3 ,WSt_3 ,WS1_3 ,WS2_3 ,WS3_3 ,z_3);
% Sigma_4 = Sigma_t(cvp_est4,Wmve_4, Wcash_4 ,Wfrqint_4 ,N4 ,WSt_4 ,WS1_4 ,WS2_4 ,WS3_4 ,z_4);
% Sigma_5 = Sigma_t(cvp_est4,Wmve_5, Wcash_5 ,Wfrqint_5 ,N5 ,WSt_5 ,WS1_5 ,WS2_5 ,WS3_5 ,z_5);
% Sigma_6 = Sigma_t(cvp_est4,Wmve_6, Wcash_6 ,Wfrqint_6 ,N6 ,WSt_6 ,WS1_6 ,WS2_6 ,WS3_6 ,z_6);
% Sigma_7 = Sigma_t(cvp_est4,Wmve_7, Wcash_7 ,Wfrqint_7 ,N7 ,WSt_7 ,WS1_7 ,WS2_7 ,WS3_7 ,z_7);
% Sigma_8 = Sigma_t(cvp_est4,Wmve_8, Wcash_8 ,Wfrqint_8 ,N8 ,WSt_8 ,WS1_8 ,WS2_8 ,WS3_8 ,z_8);
% Sigma_9 = Sigma_t(cvp_est4,Wmve_9, Wcash_9 ,Wfrqint_9 ,N9 ,WSt_9 ,WS1_9 ,WS2_9 ,WS3_9 ,z_9);
% Sigma_10 = Sigma_t(cvp_est4,Wmve_10, Wcash_10 ,Wfrqint_10 ,N10 ,WSt_10 ,WS1_10 ,WS2_10 ,WS3_10 ,z_10);
% Sigma_11 = Sigma_t(cvp_est4,Wmve_11, Wcash_11 ,Wfrqint_11 ,N11 ,WSt_11 ,WS1_11 ,WS2_11 ,WS3_11 ,z_11);
% Sigma_12 = Sigma_t(cvp_est4,Wmve_12, Wcash_12 ,Wfrqint_12 ,N12 ,WSt_12 ,WS1_12 ,WS2_12 ,WS3_12 ,z_12);
% Sigma_13 = Sigma_t(cvp_est4,Wmve_13, Wcash_13 ,Wfrqint_13 ,N13 ,WSt_13 ,WS1_13 ,WS2_13 ,WS3_13 ,z_13);
% Sigma_14 = Sigma_t(cvp_est4,Wmve_14, Wcash_14 ,Wfrqint_14 ,N14 ,WSt_14 ,WS1_14 ,WS2_14 ,WS3_14 ,z_14);
% Sigma_15 = Sigma_t(cvp_est4,Wmve_15, Wcash_15 ,Wfrqint_15 ,N15 ,WSt_15 ,WS1_15 ,WS2_15 ,WS3_15 ,z_15);
% Sigma_16 = Sigma_t(cvp_est4,Wmve_16, Wcash_16 ,Wfrqint_16 ,N16 ,WSt_16 ,WS1_16 ,WS2_16 ,WS3_16 ,z_16);
% Sigma_17 = Sigma_t(cvp_est4,Wmve_17, Wcash_17 ,Wfrqint_17 ,N17 ,WSt_17 ,WS1_17 ,WS2_17 ,WS3_17 ,z_17);
% 
% % p.d.?
% disp(min(eig(Sigma_1)))
% disp(min(eig(Sigma_2)))
% disp(min(eig(Sigma_3)))
% disp(min(eig(Sigma_4)))
% disp(min(eig(Sigma_5)))
% disp(min(eig(Sigma_6)))
% disp(min(eig(Sigma_7)))
% disp(min(eig(Sigma_8)))
% disp(min(eig(Sigma_9)))
% disp(min(eig(Sigma_10)))
% disp(min(eig(Sigma_11)))
% disp(min(eig(Sigma_12)))
% disp(min(eig(Sigma_13)))
% disp(min(eig(Sigma_14)))
% disp(min(eig(Sigma_15)))
% disp(min(eig(Sigma_16)))
% disp(min(eig(Sigma_17)))
% 
% % Variances sort of reasonable?
% disp(min(diag(Sigma_1)))
% disp(min(diag(Sigma_2)))
% disp(min(diag(Sigma_3)))
% disp(min(diag(Sigma_4)))
% disp(min(diag(Sigma_5)))
% disp(min(diag(Sigma_6)))
% disp(min(diag(Sigma_7)))
% disp(min(diag(Sigma_8)))
% disp(min(diag(Sigma_9)))
% disp(min(diag(Sigma_10)))
% disp(min(diag(Sigma_11)))
% disp(min(diag(Sigma_12)))
% disp(min(diag(Sigma_13)))
% disp(min(diag(Sigma_14)))
% disp(min(diag(Sigma_15)))
% disp(min(diag(Sigma_16)))
% disp(min(diag(Sigma_17)))
% 
% % Any big correlations?
% disp(max(max(abs(diag(sqrt(1./diag(Sigma_1)))*Sigma_1*diag(sqrt(1./diag(Sigma_1))) - eye(N1)))))
% disp(max(max(abs(diag(sqrt(1./diag(Sigma_2)))*Sigma_2*diag(sqrt(1./diag(Sigma_2))) - eye(N2)))))
% disp(max(max(abs(diag(sqrt(1./diag(Sigma_3)))*Sigma_3*diag(sqrt(1./diag(Sigma_3))) - eye(N3)))))
% disp(max(max(abs(diag(sqrt(1./diag(Sigma_4)))*Sigma_4*diag(sqrt(1./diag(Sigma_4))) - eye(N4)))))
% disp(max(max(abs(diag(sqrt(1./diag(Sigma_5)))*Sigma_5*diag(sqrt(1./diag(Sigma_5))) - eye(N5)))))
% disp(max(max(abs(diag(sqrt(1./diag(Sigma_6)))*Sigma_6*diag(sqrt(1./diag(Sigma_6))) - eye(N6)))))
% disp(max(max(abs(diag(sqrt(1./diag(Sigma_7)))*Sigma_7*diag(sqrt(1./diag(Sigma_7))) - eye(N7)))))
% disp(max(max(abs(diag(sqrt(1./diag(Sigma_8)))*Sigma_8*diag(sqrt(1./diag(Sigma_8))) - eye(N8)))))
% disp(max(max(abs(diag(sqrt(1./diag(Sigma_9)))*Sigma_9*diag(sqrt(1./diag(Sigma_9))) - eye(N9)))))
% disp(max(max(abs(diag(sqrt(1./diag(Sigma_10)))*Sigma_10*diag(sqrt(1./diag(Sigma_10))) - eye(N10)))))
% disp(max(max(abs(diag(sqrt(1./diag(Sigma_11)))*Sigma_11*diag(sqrt(1./diag(Sigma_11))) - eye(N11)))))
% disp(max(max(abs(diag(sqrt(1./diag(Sigma_12)))*Sigma_12*diag(sqrt(1./diag(Sigma_12))) - eye(N12)))))
% disp(max(max(abs(diag(sqrt(1./diag(Sigma_13)))*Sigma_13*diag(sqrt(1./diag(Sigma_13))) - eye(N13)))))
% disp(max(max(abs(diag(sqrt(1./diag(Sigma_14)))*Sigma_14*diag(sqrt(1./diag(Sigma_14))) - eye(N14)))))
% disp(max(max(abs(diag(sqrt(1./diag(Sigma_15)))*Sigma_15*diag(sqrt(1./diag(Sigma_15))) - eye(N15)))))
% disp(max(max(abs(diag(sqrt(1./diag(Sigma_16)))*Sigma_16*diag(sqrt(1./diag(Sigma_16))) - eye(N16)))))
% disp(max(max(abs(diag(sqrt(1./diag(Sigma_17)))*Sigma_17*diag(sqrt(1./diag(Sigma_17))) - eye(N17)))))
% 
% %% Try optimization step 5
% 
% [cvp_est5,~,~,~,~,hess5] = fminunc(@(theta) -ll(theta,Wmve_1, Wcash_1 ,Wfrqint_1 ,N1 ,WSt_1 ,WS1_1 ,WS2_1 ,WS3_1 ,z_1 ,v_1,...
%             Wmve_2, Wcash_2 ,Wfrqint_2 ,N2 ,WSt_2 ,WS1_2 ,WS2_2 ,WS3_2 ,z_2 ,v_2,...
%             Wmve_3, Wcash_3 ,Wfrqint_3 ,N3 ,WSt_3 ,WS1_3 ,WS2_3 ,WS3_3 ,z_3 ,v_3,...
%             Wmve_4, Wcash_4 ,Wfrqint_4 ,N4 ,WSt_4 ,WS1_4 ,WS2_4 ,WS3_4 ,z_4 ,v_4,...
%             Wmve_5, Wcash_5 ,Wfrqint_5 ,N5 ,WSt_5 ,WS1_5 ,WS2_5 ,WS3_5 ,z_5 ,v_5,...
%             Wmve_6, Wcash_6 ,Wfrqint_6 ,N6 ,WSt_6 ,WS1_6 ,WS2_6 ,WS3_6 ,z_6 ,v_6,...
%             Wmve_7, Wcash_7 ,Wfrqint_7 ,N7 ,WSt_7 ,WS1_7 ,WS2_7 ,WS3_7 ,z_7 ,v_7,...
%             Wmve_8, Wcash_8 ,Wfrqint_8 ,N8 ,WSt_8 ,WS1_8 ,WS2_8 ,WS3_8 ,z_8 ,v_8,...
%             Wmve_9, Wcash_9 ,Wfrqint_9 ,N9 ,WSt_9 ,WS1_9 ,WS2_9 ,WS3_9 ,z_9 ,v_9,...
%             Wmve_10,Wcash_10,Wfrqint_10,N10,WSt_10,WS1_10,WS2_10,WS3_10,z_10,v_10,...
%             Wmve_11,Wcash_11,Wfrqint_11,N11,WSt_11,WS1_11,WS2_11,WS3_11,z_11,v_11,...
%             Wmve_12,Wcash_12,Wfrqint_12,N12,WSt_12,WS1_12,WS2_12,WS3_12,z_12,v_12,...
%             Wmve_13,Wcash_13,Wfrqint_13,N13,WSt_13,WS1_13,WS2_13,WS3_13,z_13,v_13,...
%             Wmve_14,Wcash_14,Wfrqint_14,N14,WSt_14,WS1_14,WS2_14,WS3_14,z_14,v_14,...
%             Wmve_15,Wcash_15,Wfrqint_15,N15,WSt_15,WS1_15,WS2_15,WS3_15,z_15,v_15,...
%             Wmve_16,Wcash_16,Wfrqint_16,N16,WSt_16,WS1_16,WS2_16,WS3_16,z_16,v_16,...
%             Wmve_17,Wcash_17,Wfrqint_17,N17,WSt_17,WS1_17,WS2_17,WS3_17,z_17,v_17),...
%             cvp_est4,optimset('disp','iter','MaxFunEvals',24000));
%         
% save bcv_optim5 cvp_est5 hess5 ;      
% 
% % Still hasn't converged.  Must be close?
% 
% %% See whether these estimates make any sense        
% Sigma_1 = Sigma_t(cvp_est5,Wmve_1, Wcash_1 ,Wfrqint_1 ,N1 ,WSt_1 ,WS1_1 ,WS2_1 ,WS3_1 ,z_1);
% Sigma_2 = Sigma_t(cvp_est5,Wmve_2, Wcash_2 ,Wfrqint_2 ,N2 ,WSt_2 ,WS1_2 ,WS2_2 ,WS3_2 ,z_2);
% Sigma_3 = Sigma_t(cvp_est5,Wmve_3, Wcash_3 ,Wfrqint_3 ,N3 ,WSt_3 ,WS1_3 ,WS2_3 ,WS3_3 ,z_3);
% Sigma_4 = Sigma_t(cvp_est5,Wmve_4, Wcash_4 ,Wfrqint_4 ,N4 ,WSt_4 ,WS1_4 ,WS2_4 ,WS3_4 ,z_4);
% Sigma_5 = Sigma_t(cvp_est5,Wmve_5, Wcash_5 ,Wfrqint_5 ,N5 ,WSt_5 ,WS1_5 ,WS2_5 ,WS3_5 ,z_5);
% Sigma_6 = Sigma_t(cvp_est5,Wmve_6, Wcash_6 ,Wfrqint_6 ,N6 ,WSt_6 ,WS1_6 ,WS2_6 ,WS3_6 ,z_6);
% Sigma_7 = Sigma_t(cvp_est5,Wmve_7, Wcash_7 ,Wfrqint_7 ,N7 ,WSt_7 ,WS1_7 ,WS2_7 ,WS3_7 ,z_7);
% Sigma_8 = Sigma_t(cvp_est5,Wmve_8, Wcash_8 ,Wfrqint_8 ,N8 ,WSt_8 ,WS1_8 ,WS2_8 ,WS3_8 ,z_8);
% Sigma_9 = Sigma_t(cvp_est5,Wmve_9, Wcash_9 ,Wfrqint_9 ,N9 ,WSt_9 ,WS1_9 ,WS2_9 ,WS3_9 ,z_9);
% Sigma_10 = Sigma_t(cvp_est5,Wmve_10, Wcash_10 ,Wfrqint_10 ,N10 ,WSt_10 ,WS1_10 ,WS2_10 ,WS3_10 ,z_10);
% Sigma_11 = Sigma_t(cvp_est5,Wmve_11, Wcash_11 ,Wfrqint_11 ,N11 ,WSt_11 ,WS1_11 ,WS2_11 ,WS3_11 ,z_11);
% Sigma_12 = Sigma_t(cvp_est5,Wmve_12, Wcash_12 ,Wfrqint_12 ,N12 ,WSt_12 ,WS1_12 ,WS2_12 ,WS3_12 ,z_12);
% Sigma_13 = Sigma_t(cvp_est5,Wmve_13, Wcash_13 ,Wfrqint_13 ,N13 ,WSt_13 ,WS1_13 ,WS2_13 ,WS3_13 ,z_13);
% Sigma_14 = Sigma_t(cvp_est5,Wmve_14, Wcash_14 ,Wfrqint_14 ,N14 ,WSt_14 ,WS1_14 ,WS2_14 ,WS3_14 ,z_14);
% Sigma_15 = Sigma_t(cvp_est5,Wmve_15, Wcash_15 ,Wfrqint_15 ,N15 ,WSt_15 ,WS1_15 ,WS2_15 ,WS3_15 ,z_15);
% Sigma_16 = Sigma_t(cvp_est5,Wmve_16, Wcash_16 ,Wfrqint_16 ,N16 ,WSt_16 ,WS1_16 ,WS2_16 ,WS3_16 ,z_16);
% Sigma_17 = Sigma_t(cvp_est5,Wmve_17, Wcash_17 ,Wfrqint_17 ,N17 ,WSt_17 ,WS1_17 ,WS2_17 ,WS3_17 ,z_17);
% 
% % p.d.?
% disp(min(eig(Sigma_1)))
% disp(min(eig(Sigma_2)))
% disp(min(eig(Sigma_3)))
% disp(min(eig(Sigma_4)))
% disp(min(eig(Sigma_5)))
% disp(min(eig(Sigma_6)))
% disp(min(eig(Sigma_7)))
% disp(min(eig(Sigma_8)))
% disp(min(eig(Sigma_9)))
% disp(min(eig(Sigma_10)))
% disp(min(eig(Sigma_11)))
% disp(min(eig(Sigma_12)))
% disp(min(eig(Sigma_13)))
% disp(min(eig(Sigma_14)))
% disp(min(eig(Sigma_15)))
% disp(min(eig(Sigma_16)))
% disp(min(eig(Sigma_17)))
% 
% % Variances sort of reasonable?
% disp(min(diag(Sigma_1)))
% disp(min(diag(Sigma_2)))
% disp(min(diag(Sigma_3)))
% disp(min(diag(Sigma_4)))
% disp(min(diag(Sigma_5)))
% disp(min(diag(Sigma_6)))
% disp(min(diag(Sigma_7)))
% disp(min(diag(Sigma_8)))
% disp(min(diag(Sigma_9)))
% disp(min(diag(Sigma_10)))
% disp(min(diag(Sigma_11)))
% disp(min(diag(Sigma_12)))
% disp(min(diag(Sigma_13)))
% disp(min(diag(Sigma_14)))
% disp(min(diag(Sigma_15)))
% disp(min(diag(Sigma_16)))
% disp(min(diag(Sigma_17)))
% 
% disp(mean(diag(Sigma_1)))
% disp(mean(diag(Sigma_2)))
% disp(mean(diag(Sigma_3)))
% disp(mean(diag(Sigma_4)))
% disp(mean(diag(Sigma_5)))
% disp(mean(diag(Sigma_6)))
% disp(mean(diag(Sigma_7)))
% disp(mean(diag(Sigma_8)))
% disp(mean(diag(Sigma_9)))
% disp(mean(diag(Sigma_10)))
% disp(mean(diag(Sigma_11)))
% disp(mean(diag(Sigma_12)))
% disp(mean(diag(Sigma_13)))
% disp(mean(diag(Sigma_14)))
% disp(mean(diag(Sigma_15)))
% disp(mean(diag(Sigma_16)))
% disp(mean(diag(Sigma_17)))
% 
% % Any big correlations?
% disp(max(max(abs(diag(sqrt(1./diag(Sigma_1)))*Sigma_1*diag(sqrt(1./diag(Sigma_1))) - eye(N1)))))
% disp(max(max(abs(diag(sqrt(1./diag(Sigma_2)))*Sigma_2*diag(sqrt(1./diag(Sigma_2))) - eye(N2)))))
% disp(max(max(abs(diag(sqrt(1./diag(Sigma_3)))*Sigma_3*diag(sqrt(1./diag(Sigma_3))) - eye(N3)))))
% disp(max(max(abs(diag(sqrt(1./diag(Sigma_4)))*Sigma_4*diag(sqrt(1./diag(Sigma_4))) - eye(N4)))))
% disp(max(max(abs(diag(sqrt(1./diag(Sigma_5)))*Sigma_5*diag(sqrt(1./diag(Sigma_5))) - eye(N5)))))
% disp(max(max(abs(diag(sqrt(1./diag(Sigma_6)))*Sigma_6*diag(sqrt(1./diag(Sigma_6))) - eye(N6)))))
% disp(max(max(abs(diag(sqrt(1./diag(Sigma_7)))*Sigma_7*diag(sqrt(1./diag(Sigma_7))) - eye(N7)))))
% disp(max(max(abs(diag(sqrt(1./diag(Sigma_8)))*Sigma_8*diag(sqrt(1./diag(Sigma_8))) - eye(N8)))))
% disp(max(max(abs(diag(sqrt(1./diag(Sigma_9)))*Sigma_9*diag(sqrt(1./diag(Sigma_9))) - eye(N9)))))
% disp(max(max(abs(diag(sqrt(1./diag(Sigma_10)))*Sigma_10*diag(sqrt(1./diag(Sigma_10))) - eye(N10)))))
% disp(max(max(abs(diag(sqrt(1./diag(Sigma_11)))*Sigma_11*diag(sqrt(1./diag(Sigma_11))) - eye(N11)))))
% disp(max(max(abs(diag(sqrt(1./diag(Sigma_12)))*Sigma_12*diag(sqrt(1./diag(Sigma_12))) - eye(N12)))))
% disp(max(max(abs(diag(sqrt(1./diag(Sigma_13)))*Sigma_13*diag(sqrt(1./diag(Sigma_13))) - eye(N13)))))
% disp(max(max(abs(diag(sqrt(1./diag(Sigma_14)))*Sigma_14*diag(sqrt(1./diag(Sigma_14))) - eye(N14)))))
% disp(max(max(abs(diag(sqrt(1./diag(Sigma_15)))*Sigma_15*diag(sqrt(1./diag(Sigma_15))) - eye(N15)))))
% disp(max(max(abs(diag(sqrt(1./diag(Sigma_16)))*Sigma_16*diag(sqrt(1./diag(Sigma_16))) - eye(N16)))))
% disp(max(max(abs(diag(sqrt(1./diag(Sigma_17)))*Sigma_17*diag(sqrt(1./diag(Sigma_17))) - eye(N17)))))
% 
% %% Try optimization step 6
% 
% [cvp_est6,~,~,~,~,hess6] = fminunc(@(theta) -ll(theta,Wmve_1, Wcash_1 ,Wfrqint_1 ,N1 ,WSt_1 ,WS1_1 ,WS2_1 ,WS3_1 ,z_1 ,v_1,...
%             Wmve_2, Wcash_2 ,Wfrqint_2 ,N2 ,WSt_2 ,WS1_2 ,WS2_2 ,WS3_2 ,z_2 ,v_2,...
%             Wmve_3, Wcash_3 ,Wfrqint_3 ,N3 ,WSt_3 ,WS1_3 ,WS2_3 ,WS3_3 ,z_3 ,v_3,...
%             Wmve_4, Wcash_4 ,Wfrqint_4 ,N4 ,WSt_4 ,WS1_4 ,WS2_4 ,WS3_4 ,z_4 ,v_4,...
%             Wmve_5, Wcash_5 ,Wfrqint_5 ,N5 ,WSt_5 ,WS1_5 ,WS2_5 ,WS3_5 ,z_5 ,v_5,...
%             Wmve_6, Wcash_6 ,Wfrqint_6 ,N6 ,WSt_6 ,WS1_6 ,WS2_6 ,WS3_6 ,z_6 ,v_6,...
%             Wmve_7, Wcash_7 ,Wfrqint_7 ,N7 ,WSt_7 ,WS1_7 ,WS2_7 ,WS3_7 ,z_7 ,v_7,...
%             Wmve_8, Wcash_8 ,Wfrqint_8 ,N8 ,WSt_8 ,WS1_8 ,WS2_8 ,WS3_8 ,z_8 ,v_8,...
%             Wmve_9, Wcash_9 ,Wfrqint_9 ,N9 ,WSt_9 ,WS1_9 ,WS2_9 ,WS3_9 ,z_9 ,v_9,...
%             Wmve_10,Wcash_10,Wfrqint_10,N10,WSt_10,WS1_10,WS2_10,WS3_10,z_10,v_10,...
%             Wmve_11,Wcash_11,Wfrqint_11,N11,WSt_11,WS1_11,WS2_11,WS3_11,z_11,v_11,...
%             Wmve_12,Wcash_12,Wfrqint_12,N12,WSt_12,WS1_12,WS2_12,WS3_12,z_12,v_12,...
%             Wmve_13,Wcash_13,Wfrqint_13,N13,WSt_13,WS1_13,WS2_13,WS3_13,z_13,v_13,...
%             Wmve_14,Wcash_14,Wfrqint_14,N14,WSt_14,WS1_14,WS2_14,WS3_14,z_14,v_14,...
%             Wmve_15,Wcash_15,Wfrqint_15,N15,WSt_15,WS1_15,WS2_15,WS3_15,z_15,v_15,...
%             Wmve_16,Wcash_16,Wfrqint_16,N16,WSt_16,WS1_16,WS2_16,WS3_16,z_16,v_16,...
%             Wmve_17,Wcash_17,Wfrqint_17,N17,WSt_17,WS1_17,WS2_17,WS3_17,z_17,v_17),...
%             cvp_est5,optimset('disp','iter','MaxFunEvals',24000,'MaxIter',1000));
%         
% save bcv_optim6 cvp_est6 hess6 ;      
% 
% % Stopped right away with "local minimum possible".  Good enough for me for
% % this exercise.

%% Load in estimation results from above.  Comment either this block or the block above out.
load bcv_optim6 ;


%% See whether these estimates make any sense        
Sigma_1 = Sigma_t(cvp_est6,Wmve_1, Wcash_1 ,Wfrqint_1 ,N1 ,WSt_1 ,WS1_1 ,WS2_1 ,WS3_1 ,z_1);
Sigma_2 = Sigma_t(cvp_est6,Wmve_2, Wcash_2 ,Wfrqint_2 ,N2 ,WSt_2 ,WS1_2 ,WS2_2 ,WS3_2 ,z_2);
Sigma_3 = Sigma_t(cvp_est6,Wmve_3, Wcash_3 ,Wfrqint_3 ,N3 ,WSt_3 ,WS1_3 ,WS2_3 ,WS3_3 ,z_3);
Sigma_4 = Sigma_t(cvp_est6,Wmve_4, Wcash_4 ,Wfrqint_4 ,N4 ,WSt_4 ,WS1_4 ,WS2_4 ,WS3_4 ,z_4);
Sigma_5 = Sigma_t(cvp_est6,Wmve_5, Wcash_5 ,Wfrqint_5 ,N5 ,WSt_5 ,WS1_5 ,WS2_5 ,WS3_5 ,z_5);
Sigma_6 = Sigma_t(cvp_est6,Wmve_6, Wcash_6 ,Wfrqint_6 ,N6 ,WSt_6 ,WS1_6 ,WS2_6 ,WS3_6 ,z_6);
Sigma_7 = Sigma_t(cvp_est6,Wmve_7, Wcash_7 ,Wfrqint_7 ,N7 ,WSt_7 ,WS1_7 ,WS2_7 ,WS3_7 ,z_7);
Sigma_8 = Sigma_t(cvp_est6,Wmve_8, Wcash_8 ,Wfrqint_8 ,N8 ,WSt_8 ,WS1_8 ,WS2_8 ,WS3_8 ,z_8);
Sigma_9 = Sigma_t(cvp_est6,Wmve_9, Wcash_9 ,Wfrqint_9 ,N9 ,WSt_9 ,WS1_9 ,WS2_9 ,WS3_9 ,z_9);
Sigma_10 = Sigma_t(cvp_est6,Wmve_10, Wcash_10 ,Wfrqint_10 ,N10 ,WSt_10 ,WS1_10 ,WS2_10 ,WS3_10 ,z_10);
Sigma_11 = Sigma_t(cvp_est6,Wmve_11, Wcash_11 ,Wfrqint_11 ,N11 ,WSt_11 ,WS1_11 ,WS2_11 ,WS3_11 ,z_11);
Sigma_12 = Sigma_t(cvp_est6,Wmve_12, Wcash_12 ,Wfrqint_12 ,N12 ,WSt_12 ,WS1_12 ,WS2_12 ,WS3_12 ,z_12);
Sigma_13 = Sigma_t(cvp_est6,Wmve_13, Wcash_13 ,Wfrqint_13 ,N13 ,WSt_13 ,WS1_13 ,WS2_13 ,WS3_13 ,z_13);
Sigma_14 = Sigma_t(cvp_est6,Wmve_14, Wcash_14 ,Wfrqint_14 ,N14 ,WSt_14 ,WS1_14 ,WS2_14 ,WS3_14 ,z_14);
Sigma_15 = Sigma_t(cvp_est6,Wmve_15, Wcash_15 ,Wfrqint_15 ,N15 ,WSt_15 ,WS1_15 ,WS2_15 ,WS3_15 ,z_15);
Sigma_16 = Sigma_t(cvp_est6,Wmve_16, Wcash_16 ,Wfrqint_16 ,N16 ,WSt_16 ,WS1_16 ,WS2_16 ,WS3_16 ,z_16);
Sigma_17 = Sigma_t(cvp_est6,Wmve_17, Wcash_17 ,Wfrqint_17 ,N17 ,WSt_17 ,WS1_17 ,WS2_17 ,WS3_17 ,z_17);

% p.d.?
disp(min(eig(Sigma_1)))
disp(min(eig(Sigma_2)))
disp(min(eig(Sigma_3)))
disp(min(eig(Sigma_4)))
disp(min(eig(Sigma_5)))
disp(min(eig(Sigma_6)))
disp(min(eig(Sigma_7)))
disp(min(eig(Sigma_8)))
disp(min(eig(Sigma_9)))
disp(min(eig(Sigma_10)))
disp(min(eig(Sigma_11)))
disp(min(eig(Sigma_12)))
disp(min(eig(Sigma_13)))
disp(min(eig(Sigma_14)))
disp(min(eig(Sigma_15)))
disp(min(eig(Sigma_16)))
disp(min(eig(Sigma_17)))

% Variances sort of reasonable?
disp(min(diag(Sigma_1)))
disp(min(diag(Sigma_2)))
disp(min(diag(Sigma_3)))
disp(min(diag(Sigma_4)))
disp(min(diag(Sigma_5)))
disp(min(diag(Sigma_6)))
disp(min(diag(Sigma_7)))
disp(min(diag(Sigma_8)))
disp(min(diag(Sigma_9)))
disp(min(diag(Sigma_10)))
disp(min(diag(Sigma_11)))
disp(min(diag(Sigma_12)))
disp(min(diag(Sigma_13)))
disp(min(diag(Sigma_14)))
disp(min(diag(Sigma_15)))
disp(min(diag(Sigma_16)))
disp(min(diag(Sigma_17)))

disp(mean(diag(Sigma_1)))
disp(mean(diag(Sigma_2)))
disp(mean(diag(Sigma_3)))
disp(mean(diag(Sigma_4)))
disp(mean(diag(Sigma_5)))
disp(mean(diag(Sigma_6)))
disp(mean(diag(Sigma_7)))
disp(mean(diag(Sigma_8)))
disp(mean(diag(Sigma_9)))
disp(mean(diag(Sigma_10)))
disp(mean(diag(Sigma_11)))
disp(mean(diag(Sigma_12)))
disp(mean(diag(Sigma_13)))
disp(mean(diag(Sigma_14)))
disp(mean(diag(Sigma_15)))
disp(mean(diag(Sigma_16)))
disp(mean(diag(Sigma_17)))

% Any big correlations?
disp(max(max(abs(diag(sqrt(1./diag(Sigma_1)))*Sigma_1*diag(sqrt(1./diag(Sigma_1))) - eye(N1)))))
disp(max(max(abs(diag(sqrt(1./diag(Sigma_2)))*Sigma_2*diag(sqrt(1./diag(Sigma_2))) - eye(N2)))))
disp(max(max(abs(diag(sqrt(1./diag(Sigma_3)))*Sigma_3*diag(sqrt(1./diag(Sigma_3))) - eye(N3)))))
disp(max(max(abs(diag(sqrt(1./diag(Sigma_4)))*Sigma_4*diag(sqrt(1./diag(Sigma_4))) - eye(N4)))))
disp(max(max(abs(diag(sqrt(1./diag(Sigma_5)))*Sigma_5*diag(sqrt(1./diag(Sigma_5))) - eye(N5)))))
disp(max(max(abs(diag(sqrt(1./diag(Sigma_6)))*Sigma_6*diag(sqrt(1./diag(Sigma_6))) - eye(N6)))))
disp(max(max(abs(diag(sqrt(1./diag(Sigma_7)))*Sigma_7*diag(sqrt(1./diag(Sigma_7))) - eye(N7)))))
disp(max(max(abs(diag(sqrt(1./diag(Sigma_8)))*Sigma_8*diag(sqrt(1./diag(Sigma_8))) - eye(N8)))))
disp(max(max(abs(diag(sqrt(1./diag(Sigma_9)))*Sigma_9*diag(sqrt(1./diag(Sigma_9))) - eye(N9)))))
disp(max(max(abs(diag(sqrt(1./diag(Sigma_10)))*Sigma_10*diag(sqrt(1./diag(Sigma_10))) - eye(N10)))))
disp(max(max(abs(diag(sqrt(1./diag(Sigma_11)))*Sigma_11*diag(sqrt(1./diag(Sigma_11))) - eye(N11)))))
disp(max(max(abs(diag(sqrt(1./diag(Sigma_12)))*Sigma_12*diag(sqrt(1./diag(Sigma_12))) - eye(N12)))))
disp(max(max(abs(diag(sqrt(1./diag(Sigma_13)))*Sigma_13*diag(sqrt(1./diag(Sigma_13))) - eye(N13)))))
disp(max(max(abs(diag(sqrt(1./diag(Sigma_14)))*Sigma_14*diag(sqrt(1./diag(Sigma_14))) - eye(N14)))))
disp(max(max(abs(diag(sqrt(1./diag(Sigma_15)))*Sigma_15*diag(sqrt(1./diag(Sigma_15))) - eye(N15)))))
disp(max(max(abs(diag(sqrt(1./diag(Sigma_16)))*Sigma_16*diag(sqrt(1./diag(Sigma_16))) - eye(N16)))))
disp(max(max(abs(diag(sqrt(1./diag(Sigma_17)))*Sigma_17*diag(sqrt(1./diag(Sigma_17))) - eye(N17)))))

%% Need to fit a distribution for the error terms now.
% Assume Sigma_t^(-1/2) v_t \sim iid (0,I_{N_t})
eps1 = chol(inv(Sigma_1))*v_1;
eps2 = chol(inv(Sigma_2))*v_2;
eps3 = chol(inv(Sigma_3))*v_3;
eps4 = chol(inv(Sigma_4))*v_4;
eps5 = chol(inv(Sigma_5))*v_5;
eps6 = chol(inv(Sigma_6))*v_6;
eps7 = chol(inv(Sigma_7))*v_7;
eps8 = chol(inv(Sigma_8))*v_8;
eps9 = chol(inv(Sigma_9))*v_9;
eps10 = chol(inv(Sigma_10))*v_10;
eps11 = chol(inv(Sigma_11))*v_11;
eps12 = chol(inv(Sigma_12))*v_12;
eps13 = chol(inv(Sigma_13))*v_13;
eps14 = chol(inv(Sigma_14))*v_14;
eps15 = chol(inv(Sigma_15))*v_15;
eps16 = chol(inv(Sigma_16))*v_16;
eps17 = chol(inv(Sigma_17))*v_17;

eps = [eps1;eps2;eps3;eps4;eps5;eps6;eps7;eps8;eps9;eps10;eps11;eps12;eps13;eps14;eps15;eps16;eps17];

% for ii = 1:17
%     name = strcat('eps',num2str(ii));
%     figure; hist(eval(name),20);
% end
histfit(eps);

% Looks like a little asymmetry and lots of kurtosis.  T?  EGB2?
% Need to constrain mean and variance to 0 and 1.

loglikt = @(df,data) -sum(log(sqrt(df/(df-2))*tpdf(sqrt(df/(df-2))*data,df)));

[dfhat,~,~,~,~,~,dfhess] = fmincon(@(param) loglikt(param,eps),8,[],[],[],[],3,Inf,[],optimset('disp','iter'));

loglikegb2 = @(param,data) -((param(3)/param(2))*sum(data-param(1))-size(data,1)*log(param(2))-size(data,1)*betaln(param(3),param(4))-(param(3)+param(4))*sum(log(1+exp((data-param(1))./param(2)))));

[egb2hat,~,~,~,~,~,egb2hess] = fmincon(@(param) loglikegb2(param,eps),[0;1;2;2],[],[],[],[],[-Inf;1e-8;1e-8;1e-8],Inf*ones(4,1),[],optimset('disp','iter'));
[egb2hatC,~,~,~,~,~,egb2hessC] = fmincon(@(param) loglikegb2(param,eps),egb2hat,[],[],[],[],[-Inf;1e-8;1e-8;1e-8],Inf*ones(4,1),@(param) egb2con(param),optimset('disp','iter'));

pd = fitdist(eps,'kernel');

plot(-6:.01:6,[pdf(pd,-6:.01:6)' , sqrt(dfhat/(dfhat-2))*tpdf(sqrt(dfhat/(dfhat-2))*(-6:.01:6),dfhat)' , egb2_pdf(-6:.01:6,egb2hatC)']);

% EGB2 looks pretty good.

save distparams egb2hatC ;

SIGMA = {Sigma_1, Sigma_2, Sigma_3, Sigma_4, Sigma_5, Sigma_6, Sigma_7, ...
    Sigma_8, Sigma_9, Sigma_10, Sigma_11, Sigma_12, Sigma_13, Sigma_14, ...
    Sigma_15, Sigma_16, Sigma_17} ;

save SigmaEsts SIGMA ;


% So, we've got an estimated model.  Now we need to use this to generate
% data and make see how different inference procedures work. Have to keep
% all observations lined up by firm, time, state, etc.  
