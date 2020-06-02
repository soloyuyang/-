function result = spatial_hausman(y,x,W,sl)
% PURPOSE: Spatial Hausman Test for differences in 
%          OLS and SEM coefficient estimates
% ---------------------------------------------------
%  USAGE: result = spatial_hausman(y,x,W,sl)
%  where: y = dependent variable vector
%         x = independent variables matrix
%         W = contiguity matrix (standardized)
%         sl = significance level for chi-squared probablility 
%               (e.g. .01 for 1%, .05 for 5%, etc.)
% ---------------------------------------------------
%  RETURNS: a structure variable
%         result.meth = 'spatial_hausman'
%         result.sh   = Hausman statistic
%         result.prob = marginal probability
%         result.chi1 = Critical values for chi squared at sig level
%                       and k degrees of freedom
%         result.nobs = # of observations
%         result.nvar = # of variables
% ---------------------------------------------------
% NOTE: sh > chi squared at sig level, => small prob,
%                                       => reject HO: OLS and SEM coefficient
%                                                     estimates differ
% ---------------------------------------------------
% See also:  walds, lratios, moran
% ---------------------------------------------------
% REFERENCES: LeSage and Pace. "A Spatial Hausman Test" Economic Letters
%             vol. 101, no. 3. December 2008, pp. 282-284
% ---------------------------------------------------

% written by:

% Donald J. Lacombe
% Associate Professor
% Department of Personal Financial Planning
% Texas Tech University
% donald.lacombe@ttu.edu

% Based on the lmerror code of James P. LeSage
% James P. LeSage
% jlesage@spatial-econometrics.com

if nargin ~= 4
error('Wrong # of arguments to spatial_hausman');
end;

% Get OLS coefficients
beta_ols = x\y;

% Get SEM Coefficients
options.lflag = 0;
res_sem = sem(y,x,W,options);
beta_sem = res_sem.beta;

[n k] = size(x);
H = inv(x'*x)*x';
In = eye(n);
rho = res_sem.rho;
A = (In-rho*W);
Ainv = inv(In-rho*W);
Binv = inv(In-rho*W');
sig = res_sem.sige;
d = beta_ols - beta_sem;
omega_ols = sig*H*Ainv*Binv*H';
omega_sem = sig*inv(x'*A'*A*x);
diff = omega_ols - omega_sem;
sh = d'*inv(diff)*d;
prob = 1-chis_prb(sh,k);
temp = 1-sl;
chi1 = chis_inv(temp,k);

result.meth = 'spatial_hausman';
result.sh   = sh;
result.prob = prob;
result.chi1 = chi1;
result.nobs = n;
result.nvar = k;
result.sl   = sl;


















