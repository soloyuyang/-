%% Load Data and Define Time Periods
clear
clc

% Matlab data files
load A
load W1

% *** Number of time periods ***
T=29; % NOTE: This value should be one less 
%             than the total number of time periods in the sample.
% ******************************

N=46; % Number of regions in the sample

%% Row-Normalize W and define Dependent and Explanatory Variables
W=normw(W1); % function of LeSage
y=A(:,3); % column number in the data matrix that corresponds to the dependent variable
x=A(:,[4,6]); % column numbers in the data matrix that correspond to the independent variables

for t=1:T+1
    t1=(t-1)*N+1;
    t2=t*N;
    wx(t1:t2,:)=W*x(t1:t2,:);
    Wy(t1:t2,1)=W*y(t1:t2,1);
end

[nobs,K]=size(x);
xconstant=ones(nobs,1);

%% Set Flags
info.lflag=0;   % no log-determinant approximation
info.tl=1;      % add time lag
info.stl=1;     % add space-time lag?
info.ted=1;     % set ted=0 for model with spatial fixed effects without time dummies in combination with sar_jihai OR
% set ted=1 for model with spatial and time period fixed effects in combination with sar_jihai_time
%info.dyn=0;     % not needed for estimating dynamic spatial panel data model, since lines 478-491 are turned off in sar_panel_FE
info.model=3;   % space and time fixed effects type "help sar_panel_FE" for other options
info.fe=0;      % do not print fixed effects
%info.bc=0;      % 0 = no bias correction OR 1 = bias correction in sar_panel_FE; not needed for estimating dynamic spatial panel data model
% since bias correction is programmed in sar_jihai_time
% file

%% Estimate Model
results=sar_panel_FE(y(N+1:end),[y(1:end-N) Wy(1:end-N) x(N+1:end,:) wx(N+1:end,:)],W,T,info);
vnames=char('logcit','timelag logcit','spacetimelag logcit','logp','logy','W*logp','W*logy');
results1=sar_jihai_time(y(1:nobs),[x(N+1:nobs,:) wx(N+1:nobs,:)],W,info); % use of sar_jihai_time implies time fixed effects next to spatial fixed effects
results.beta=results1.theta1(1:end-2);
results.rho=results1.theta1(end-1);
results.tstat=results1.tstat1(1:end-1);
prt_sp(results,vnames,1);
btemp=results1.theta1;
varcov=results1.varcov;

%% Simulate Parameter Values

NSIM=1000;  % number of simulations for effects estimates
px=size(x,2);
[npar,dummy]=size(btemp);
simresults=zeros(npar-1,NSIM);
simdirst=zeros(px,NSIM);
simindst=zeros(px,NSIM);
simtotst=zeros(px,NSIM);
simdirlt=zeros(px,NSIM);
simindlt=zeros(px,NSIM);
simtotlt=zeros(px,NSIM);
simdirc=zeros(1,NSIM);
simindc=zeros(1,NSIM);
simtotc=zeros(1,NSIM);
stability = zeros(NSIM,1);

for sim=1:NSIM
    parms = chol(varcov)'*randn(size(btemp)) + btemp;
    deltasim = parms(npar-1,1); % coef WY(t)
    betasim = parms(3:npar-2,1);
    tausim = parms(1,1); % Coef Y(t-1)
    etasim = parms(2,1); % Coef WY(t-1)
    simresults(:,sim)=[tausim;etasim;betasim;deltasim];
	stability(sim,1) = simresults(1,sim)+simresults(2,sim)+simresults(npar-1,sim);
    SS=(eye(N)-deltasim*W)\eye(N);
    SC=SS*((tausim-1)*eye(N)+(deltasim+etasim)*W);
    simdirc(1,sim)=sum(diag(SC))/N; % average direct effect
    simindc(1,sim)=sum(sum(SC,2)-diag(SC))/N; % average indirect effect
    simtotc(1,sim)=simdirc(1,sim)+simindc(1,sim);
    for p=1:px
        C=zeros(N,N);
        for i=1:N
            for j=1:N
                if (i==j) C(i,j)=betasim(p);
                else
                    C(i,j)=betasim(p+2)*W(i,j);
                end
            end
        end
        SC=SS*C;
        simdirst(p,sim)=sum(diag(SC))/N; % average direct effect
        simindst(p,sim)=sum(sum(SC,2)-diag(SC))/N; % average indirect effect
        simtotst(p,sim)=simdirst(p,sim)+simindst(p,sim);
        SC=((1-tausim)*eye(N)-(deltasim+etasim)*W)\C;
        simdirlt(p,sim)=sum(diag(SC))/N; % average direct effect
        simindlt(p,sim)=sum(sum(SC,2)-diag(SC))/N; % average indirect effect
        simtotlt(p,sim)=simdirlt(p,sim)+simindlt(p,sim);
    end
end

%% Print Results
fprintf(1,'Convergence Effect of Dependent Variable\n');
convergence_effects = [mean(simdirc,2) mean(simdirc,2)./std(simdirc,0,2) tdis_prb((mean(simdirc,2)./std(simdirc,0,2)),N*T-K-1) ;
    mean(simindc,2) mean(simindc,2)./std(simindc,0,2) tdis_prb((mean(simindc,2)./std(simindc,0,2)),N*T-K-1);
    mean(simtotc,2) mean(simtotc,2)./std(simtotc,0,2) tdis_prb((mean(simtotc,2)./std(simtotc,0,2)),N*T-K-1)];

info.cnames = char('Effect','t-stat','p-value');
info.rnames = char(' ','Direct','Indirect','Total');
mprint(convergence_effects,info)

%% Print Out Stability Statistics
mean_stability = mean(stability);
temp = cr_interval(stability,.95);
interval_stability = [temp(2) temp(1)];
output = [mean_stability interval_stability];
info.rnames = char(' ','delta + eta + tau');
info.cnames = char('mean','Lower 95%','Upper 95%');
fprintf('Stability Statistics \n')
mprint(output,info)

%% Short Term Direct Effects
effects_st_direct = [mean(simdirst,2) mean(simdirst,2)./std(simdirst,0,2) tdis_prb((mean(simdirst,2)./std(simdirst,0,2)),N*T-K-1)];

fprintf(1,'Short term direct effects \n');
info.rnames = char('Variable','log price','log income');
info.cnames = char('direct','t-stat','p-value');
mprint(effects_st_direct,info)


%% Short Term Indirect Effects
effects_st_indirect = [mean(simindst,2) mean(simindst,2)./std(simindst,0,2) tdis_prb((mean(simindst,2)./std(simindst,0,2)),N*T-K-1)];

fprintf(1,'Short term indirect effects \n');
info.rnames = char('Variable','log price','log income');
info.cnames = char('indirect','t-stat','p-value');
mprint(effects_st_indirect,info)

%% Short Term Total Effects
effects_st_total = [mean(simtotst,2) mean(simtotst,2)./std(simtotst,0,2) tdis_prb((mean(simtotst,2)./std(simtotst,0,2)),N*T-K-1)];

fprintf(1,'Short term total effects \n');
info.rnames = char('Variable','log price','log income');
info.cnames = char('total','t-stat','p-value');
mprint(effects_st_total,info)

%% Long Term Direct Effects
effects_lt_direct = [mean(simdirlt,2) mean(simdirlt,2)./std(simdirlt,0,2) tdis_prb((mean(simdirlt,2)./std(simdirlt,0,2)),N*T-K-1)];

fprintf(1,'Long term direct effects \n');
info.rnames = char('Variable','log price','log income');
info.cnames = char('direct','t-stat','p-value');
mprint(effects_lt_direct,info)

%% Long Term Indirect Effects
effects_lt_indirect = [mean(simindlt,2) mean(simindlt,2)./std(simindlt,0,2) tdis_prb((mean(simindlt,2)./std(simindlt,0,2)),N*T-K-1)];

fprintf(1,'Long term indirect effects \n');
info.rnames = char('Variable','log price','log income');
info.cnames = char('indirect','t-stat','p-value');
mprint(effects_lt_indirect,info)


%% Long Term Total Effects
effects_lt_total = [mean(simtotlt,2) mean(simtotlt,2)./std(simtotlt,0,2) tdis_prb((mean(simtotlt,2)./std(simtotlt,0,2)),N*T-K-1)];

fprintf(1,'Long term total effects \n');
info.rnames = char('Variable','log price','log income');
info.cnames = char('total','t-stat','p-value');
mprint(effects_lt_total,info)
