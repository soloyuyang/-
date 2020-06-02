function result = lmsec(y,x,W)
% PURPOSE: LM spatial error components test: Tests for the presence
%          of spatial error components
% ---------------------------------------------------
%  USAGE: result = lmsec(y,x,W)
%  where: y = dependent variable vector
%         x = independent variables matrix
%         W = contiguity matrix (standardized)
% ---------------------------------------------------
%  RETURNS: a structure variable
%         result.meth = 'lmsec'
%         result.lm   = lmsec statistic
%         result.prob = marginal probability
%         result.chi1 = 6.635 (chi-squared 1 dof at 99% level)
%         result.nobs = # of observations
%         result.nvar = # of variables
% ---------------------------------------------------
% NOTE: lm > 6.64,  => small prob,
%                    => reject HO: classic regression model with
%                                  uncorrelated and homoskedastic errors
% ---------------------------------------------------
% See also:  walds, lratios, moran
% ---------------------------------------------------
% REFERENCES: Anselin, L. and R. Moreno. Properties of Tests for Spatial
%             Error Components. Regional Science and Urban Economics
%             vol. 33, no. 5, September 2003, pp. 595-618.
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
error('Wrong # of arguments to lmsec');
end;

[n k] = size(x); 

xpxi = (x'*x)\eye(k);
M = eye(n) - x*xpxi*x';
e = M*y;
sighat = (e'*e)/n;

T1 = trace(W*W');
T2 = trace((W*W')*(W*W'));
num = ((e'*(W*W')*e)/sighat - T1)^2;
den = 2*(T2-(T1^2/n));
lmsecstat = num/den;
prob = 1-chis_prb(lmsecstat,1);

result.meth = 'lmsec';
result.lm = lmsecstat;
result.prob = prob;
result.chi1   = 6.64;
result.nobs = n;
result.nvar = k;
