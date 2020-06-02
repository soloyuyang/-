% PURPOSE: An example of using all of the specification tests
%          LM Lag, LM Error, LM Lag Robust, LM Error Robust
%          LM Lag and Error Combined, LM Spatial Error Components
%          and Spatial Hausman Test
%---------------------------------------------------
% USAGE: specification_tests_d
%---------------------------------------------------

clear all;

% W-matrix from Anselin's neigbhorhood crime data set
load anselin.dat; % standardized 1st-order spatial weight matrix
latt = anselin(:,4);
long = anselin(:,5);
[junk W junk] = xy2cont(latt,long);
[n junk] = size(W);
IN = eye(n); 
rho = 0.7;  % true value of rho
sige = 0.1;
k = 3;
x = [ones(n,1) randn(n,k) ];
betav = zeros(k+1,1);
betav(1,1) = 1.0;
betav(2,1) = 2.0;
betav(3,1) = 3.0;
betav(4,1) = 4.0;

%Generate SAR Model

y = inv(IN-rho*W)*x*betav + (IN-rho*W)\randn(n,1)*sqrt(sige); 


% Use the following to generate an OLS model
% y = x*betav + randn(n,1)*sqrt(sige); 

% Use the following to generate a SEM model
% x*betav + (IN-rho*W)\randn(n,1)*sqrt(sige); 

% LM Lag test

res1 = lmlag(y,x,W);
prt_tests(res1);

% LM Error test

res2 = lmerror(y,x,W);
prt_tests(res2);

% LM Lag Robust

res3 = lmlag_robust(y,x,W);
prt_tests(res3);

% LM Error Robust

res4 = lmerror_robust(y,x,W);
prt_tests(res4);

% LM Lag and Error Combined

res5 = lmlagerr(y,x,W);
prt_tests(res5);

% LM Spatial Error Components

res6 = lmsec(y,x,W);
prt_tests(res6);

res7 = spatial_hausman(y,x,W,.05);

% The .05 above will generate a critical value for the 
% chi squared distribution at the 5% level with k degrees of freedom
% where k is the number of regressors including the constant. In this
% example k=4

prt_tests(res7);
