function result = lmerror(y,x,W)
% PURPOSE: LM error statistic for spatial correlation in residuals
%          of a regression model
% ---------------------------------------------------
%  USAGE: result = lmerror(y,x,W)
%  where: y = dependent variable vector
%         x = independent variables matrix
%         W = contiguity matrix (standardized)
% ---------------------------------------------------
%  RETURNS: a structure variable
%         result.meth = 'lmerror'
%         result.lm   = LM statistic
%         result.prob = marginal probability
%         result.chi1 = 6.635 (chi-squared 1 dof at 99% level)
%         result.nobs = # of observations
%         result.nvar = # of variables
% ---------------------------------------------------
% NOTE: lm > 6.635,  => small prob,
%                    => reject HO: of no spatial correlation
% ---------------------------------------------------
% See also:  walds, lratios, moran
% ---------------------------------------------------
% REFERENCES: Anselin (1988), pages 103-104.
% ---------------------------------------------------

% written by:

% Donald J. Lacombe
% Associate Professor
% Department of Personal Financial Planning
% Texas Tech University
% donald.lacombe@ttu.edu

% Based on code written by
% James P. LeSage
% jlesage@spatial-econometrics.com


[n k] = size(x);

xpxi = (x'*x)\eye(k);     % Faster way of computing inv(x'*x)
b = xpxi*(x'*y);          % OLS coefficients
M = eye(n) - x*xpxi*x';
e = M*y;                  % Calculate residuals
sighat = (e'*e)/n;        % Calculate sigma hat


t1 = trace((W+W')*W);
lm1 = (e'*W*e)/sighat;
lmerr = (lm1*lm1)*(1/t1);
prob = 1-chis_prb(lmerr,1);

result.meth = 'lmerror';
result.lm = lmerr;
result.prob = prob;
result.chi1   = 6.635;
result.nobs = n;
result.nvar = k;



