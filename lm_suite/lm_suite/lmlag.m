function result = lmlag(y,x,W)
% PURPOSE: LM lag statistic for omitted spatial lag
%          of a regression model
% ---------------------------------------------------
%  USAGE: result = lmlag(y,x,W)
%  where: y = dependent variable vector
%         x = independent variables matrix
%         W = contiguity matrix (standardized)
% ---------------------------------------------------
%  RETURNS: a structure variable
%         result.meth = 'lmlag'
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
% jlesage@spatial-econometrics.com

if nargin ~= 3
error('Wrong # of arguments to lmlag');
end;

[n k] = size(x); 

xpxi = (x'*x)\eye(k);     % Faster way of computing inv(x'*x)
b = xpxi*(x'*y);          % OLS coefficients
M = eye(n) - x*xpxi*x';   % M matrix
e = M*y;                  % Calculate residuals
sighat = (e'*e)/n;        % Calculate sigma hat

T = trace((W'+W)*W);                          % T in Florax et al. paper
term1 = [(W*x*b)'*M*(W*x*b)+(T*sighat)];      % Term 2 in equation (6) in Florax et al. paper
J = (1/(n*sighat))*term1;                     % Equation (6) in Florax et al. paper
lm1 = (e'*W*y)/sighat;                        % Numerator of equation (5) in Florax et al. paper
lmlag = (lm1*lm1)*(1/(n*J));                  % Equation (5) in Florax et al. paper
prob = 1-chis_prb(lmlag,1);


result.meth = 'lmlag';
result.lm = lmlag;
result.prob = prob;
result.chi1   = 6.64;
result.nobs = n;
result.nvar = k;

