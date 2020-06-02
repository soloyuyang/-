function result = lmlagerr(y,x,W)
% PURPOSE: LM combined test for spatial lag and spatial error correlation
%          in residuals of a regression model
% ---------------------------------------------------
%  USAGE: result = lmlagerr(y,x,W)
%  where: y = dependent variable vector
%         x = independent variables matrix
%         W = contiguity matrix (standardized)
% ---------------------------------------------------
%  RETURNS: a structure variable
%         result.meth = 'lmlagerr'
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
% REFERENCES: Anselin, Bera, Florax, Yoon. "Simple Diagnostic
% Tests for Spatial Dependence" Regional Science and Urban Economics,
% vol. 26, no. 1, February, 1996, pp. 77-104. 
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


if nargin ~= 3
error('Wrong # of arguments to lmlagerr');
end;

[n k] = size(x); 

xpxi = (x'*x)\eye(k);
b = xpxi*(x'*y);
M = eye(n) - x*xpxi*x';
e = M*y;
sighat = (e'*e)/n;

% First calculate LM Lag Robust

T = trace((W'+W)*W);
term1 = [(W*x*b)'*M*(W*x*b)+(T*sighat)];
J = (1/(n*sighat))*term1;

lm1 = (e'*W*e/sighat);
lm2 = (e'*W*y/sighat);
lmr1 = (lm1 - lm2);
lmr2 = lmr1*lmr1;
den = n*J - T;
lmrhor = lmr2/den;

% Next Calculate LM Error

lm1 = (e'*W*e)/sighat;
lmerr = (lm1*lm1)*(1/T);

% Calculate probability

stat = lmrhor + lmerr;
prob = 1-chis_prb(stat,2);

result.meth = 'lmlagerr';
result.lm = stat;
result.prob = prob;
result.chi1   = 9.21;
result.nobs = n;
result.nvar = k;















