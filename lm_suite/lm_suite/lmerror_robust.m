function result = lmerror_robust(y,x,W)
% PURPOSE: LM robust error test: Test for the presence of a spatial
%          error process when a spatially lagged dependent variable is 
%          present.
% ---------------------------------------------------
%  USAGE: result = lmerror_robust(y,x,W)
%  where: y = dependent variable vector
%         x = independent variables matrix
%         W = contiguity matrix (standardized)
% ---------------------------------------------------
%  RETURNS: a structure variable
%         result.meth = 'lmerror_robust'
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
% REFERENCES: Florax, Folmer, and Rey. "Specification Searches 
% in Spatial Econometrics: The Relevance of Hendry's Methodology"
% Regional Science and Urban Economics, vol. 33, no. 5 September 2003.
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
error('Wrong # of arguments to lmerror_robust');
end;

[n k] = size(x); 

xpxi = (x'*x)\eye(k);
b = xpxi*(x'*y);
M = eye(n) - x*xpxi*x';
e = M*y;
sighat = (e'*e)/n;

T = trace((W'+W)*W);
term1 = [(W*x*b)'*M*(W*x*b)+(T*sighat)];
J = (1/(n*sighat))*term1;

lm1 = (e'*W*e/sighat);
lm2 = T*inv(n*J);
lm3 = (e'*W*y/sighat);
lmr1 = (lm1 - (lm2*lm3));
lmr2 = lmr1*lmr1;
den = T*(1-T*inv(n*J));
lmlambdar = lmr2/den;
prob = 1-chis_prb(lmlambdar,1);

result.meth = 'lmerror_robust';
result.lm = lmlambdar;
result.prob = prob;
result.chi1   = 6.64;
result.nobs = n;
result.nvar = k;
