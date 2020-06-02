function results = log_marginal_panelprob(y,xo,W,N,T,model,out,prior,rmin,rmax,incr)
% PURPOSE: Bayesian log-marginal posterior for all static spatial panel models
%          user should eliminate fixed effects using differencing, de-meaning transformations
%          no priors on beta, sige
%          beta(c,d) prior on rho (see LeSage and Pace 2009)
%-------------------------------------------------------------
% USAGE: results = log_marginal_panelb(y,x,W1,W2,N,T,prior)
% where: y = dependent variable vector (N*T x 1)
%        x = independent variables matrix, WITHOUT INTERCEPT TERM 
%        W1 = N by N spatial weight matrix (for W1*y)
%        W2 = N by N spatial weight matrix (for W2*e)
%             (W1 can be equal to W2)
%        N = # of cross-sectional units
%        T = # of time periods
%    prior = a structure variable with:
%            prior.c    = beta(c,d) prior parameter for rho (default: 1.01)
%            prior.d    = beta(c,d) prior parameter for rho (default: 1.01)
%-------------------------------------------------------------
% RETURNS:  a structure:
%          results.meth   = 'log_margainal_panel'
%          results.nobs   = # of cross-sectional observations
%          results.ntime  = # of time periods
%          results.nvar   = # of variables in x-matrix
%          results.y      = y-vector from input (N*T x 1)
%          results.c      = c prior parameter (from input)
%          results.d      = d prior parameter (from input)
%          results.time   = time for log-marginal posterior calculation
%          results.lmarginal = a 6 x 1 column-vector with [log-marginal]
%          results.probs  = a 6 x 1 column-vector with model probs
%          results.logm_ols
%          results.logm_slx
%          results.logm_sar
%          results.logm_sem
%          results.logm_sdm
%          results.logm_sdem
% --------------------------------------------------------------
% NOTES: - returns only the log-marginal posterior and probabilities for model comparison purposes
%          NO ESTIMATES returned
% - results.lmarginal can be used for model comparison 
% - uses Gauss-Legendre (90-point) integration over -1 to 1 interval for dependence parameters
% --------------------------------------------------------------
% See also: log_marginal_panel.m (that includes the SAC model)

% Koop (2003, p 42): "When comparing models using posterior odds ratios, it is
% acceptable to use noninformative priors over parameters which are common
% to all the models. However, informative, proper priors should be used over all
% other parameters."

% This function uses a beta(c,d) prior on rho, but no priors on beta,sigma
% 
% 
% written by:
% James P. LeSage, last updated 6/2013
% Dept of Finance & Economics
% Texas State University-San Marcos
% 601 University Drive
% San Marcos, TX 78666
% jlesage@spatial-econometrics.com

%
% Adjusted by J.Paul Elhorst
% University of Groningen
% Autumn 2013
% Adjustment 1: program uses less memory space by not blowing up the W
% matrix to an N*T by N*T matrix
% Adjustment 2: correction for degrees of freedom in the fixed effects
% model has been added

timet = clock; % start the timer
tmp = rmin:incr:rmax;
rgrid = tmp';
ngrid = length(rgrid);

% set default prior parameters
c = 1.01;
d = 1.01;

fields = fieldnames(prior);
nf = length(fields);
if nf > 0
    for i=1:nf
        if strcmp(fields{i},'c')
                    c = prior.c;
        elseif strcmp(fields{i},'d')
            d = prior.d;
        end;
    end;
else  % the user has input a blank info structure
      % so we use the defaults
end; 

% Adjustment for degrees of freedom due to included fixed effects 
if (model==0) dofadj=1; % intercept will be included
elseif (model==1) dofadj=N; % correction for N spatial fixed effects
elseif (model==2) dofadj=T; % correction for T time-period fixed effects
elseif (model==3) dofadj=N+T-1; % correction for spatial and time-period fixed effects
else error('model is not correctly specified');
end

[nt,nx] = size(xo);

results.nobs  = N;
results.ntime = T;
results.nvar  = nx;
results.y = y;   
results.c = c;
results.d = d;

IN = speye(N);
IT = speye(T);

%out = lndetfull(W,rmin,rmax);
tt=rmin:.001:rmax; % interpolate a finer grid
outi = interp1(out.rho,out.lndet,tt','spline');
detval = [tt' outi];

% bprior1 = beta_prior(rgrid(1:45,1),c,d);
% bprior2 = beta_prior(rgrid(46:end,1),c,d);
% 
% bprior1 = -bprior1 + 0.5;
% bprior2 = -bprior2 + 0.5;
% 
% bprior1(45,1) = 10000;
% bprior2(1,1) = 10000;
% 
% bprior = [bprior1
%           bprior2];

Wx=xo;
WWx=xo;
Wy=y;
for t=1:T
    t1=1+(t-1)*N;t2=t*N;
    Wx(t1:t2,:)=W*xo(t1:t2,:);
    WWx(t1:t2,:)=W*Wx(t1:t2,:);
    Wy(t1:t2,1)=W*y(t1:t2,1);
end

if (model==0)
    xsar = [ones(nt,1) xo];
    xsdm = [ones(nt,1) xo Wx];
    xsem = [ones(nt,1) xo];
else
    xsar = xo;
    xsdm = [xo Wx];
    xsem = xo;
end
% ====================================================================
% log-marginal for least-squares case based on non-informative priors
x = xsar;
dof = (N*T - nx - dofadj)/2;

lndetx_ols = log(det(x'*x));

logC_ols =  gammaln(dof) - dof*log(2*pi) - 0.5*lndetx_ols;

bhat = (x'*x)\(x'*y);
epe = (y-x*bhat)'*(y-x*bhat);
logm_out = logC_ols - dof*log(epe);
results.logm_ols = logm_out;
% 
% % ====================================================================
% % log-marginal for SLX case based on non-informative priors
% x = xsdm;
% dof = (N*T -2*nx -1)/2;
% 
% lndetx_slx = log(det(x'*x));
% logC_slx = gammaln(dof) - dof*log(2*pi) - 0.5*lndetx_slx;
% 
% bhat = (x'*x)\(x'*y);
% epe = (y-x*bhat)'*(y-x*bhat);
% logm_out = logC_slx - dof*log(epe);
% results.logm_slx = logm_out;

% ====================================================================
% evaluate log-marginal for SAR model over a grid of rho values
% informative beta-prior on rho

x = xsar;
xpx = x'*x;

lndetx_sar = log(det(xpx));
dof = (N*T - nx - dofadj)/2;

D = (1 - 1/rmin); % from uniform prior on rho
logC_sar = -log(D) + gammaln(dof) - dof*log(2*pi) - 0.5*lndetx_sar;

bo = (xpx)\(x'*y);
bd = (xpx)\(x'*Wy);
eo = y - x*bo;
ed = Wy - x*bd;
epeo = eo'*eo;
eped = ed'*ed;
epeod = ed'*eo;

iota = ones(ngrid,1);

Q1 = epeo*iota - 2*rgrid*epeod + (rgrid.*rgrid)*eped;
logm_sar = zeros(ngrid,1);
Q2 = detval(:,2);

logm_sar = -dof*log(Q1) + T*Q2;

adj = max(logm_sar);
results.maxsar = adj;
madj = logm_sar - adj;
xx = exp(madj);

% trapezoid rule integration
yy =rgrid;
isum = sum((yy(2:ngrid,1) + yy(1:ngrid-1,1)).*(xx(2:ngrid,1) - xx(1:ngrid-1,1))/2);

logm_out = isum + adj + logC_sar;

results.logm_sar = logm_out;


% ====================================================================
% evaluate log-marginal for SDM model over a grid of rho values

x = xsdm;
xpx = x'*x;

lndetx_sdm = log(det(xpx));
dof = (N*T - 2*nx - dofadj)/2;

D = (1 - 1/rmin); % from uniform prior on rho

logC_sdm = -log(D) + gammaln(dof) - dof*log(2*pi) - 0.5*lndetx_sdm;

bo = (xpx)\(x'*y);
bd = (xpx)\(x'*Wy);
eo = y - x*bo;
ed = Wy - x*bd;
epeo = eo'*eo;
eped = ed'*ed;
epeod = ed'*eo;

iota = ones(ngrid,1);

Q1 = epeo*iota - 2*rgrid*epeod + rgrid.*rgrid*eped;
Q2 = detval(:,2);

logm_sdm = - dof*log(Q1) + T*Q2;

adj = max(logm_sdm);
results.maxsdm = adj;
madj = logm_sdm - adj;
xx = exp(madj);
% trapezoid rule integration
yy =rgrid;
isum = sum((yy(2:ngrid,1) + yy(1:ngrid-1,1)).*(xx(2:ngrid,1) - xx(1:ngrid-1,1))/2);

logm_out = isum + adj + logC_sdm;

results.logm_sdm = logm_out;

% ====================================================================
% evaluate log-marginal for SEM model over a grid of rho values

dof = (N*T - nx - dofadj)/2;

D = (1 - 1/rmin); % from uniform prior on rho

logC_sem = -log(D) + gammaln(dof) - dof*log(2*pi); 


% do vectorized calculations
x = xsem;

xpx = x'*x;
xpWx = x'*Wx;
xpWpx = Wx'*x;
xpWpWx = Wx'*Wx;
xpy = x'*y;
xpWy = x'*Wy;
xpWpy = Wx'*y;
xpWpWy = Wx'*Wy;
ypy = y'*y;
ypWy = y'*Wy;
ypWpy = Wy'*y;
ypWpWy = Wy'*Wy;

Q1 = zeros(ngrid,1);
Q3 = zeros(ngrid,1);

for iter=1:ngrid;
    rho = rgrid(iter,1);
    
 Axx = xpx - rho*xpWx - rho*xpWpx + rho*rho*xpWpWx;
 
 Q3(iter,1) = log(det(Axx));
 
 Axy = xpy - rho*xpWy - rho*xpWpy + rho*rho*xpWpWy;

 Ayy = ypy - rho*ypWy - rho*ypWpy + rho*rho*ypWpWy;

 b = Axx\Axy;

 Q1(iter,1) = Ayy - b'*Axx*b;

end;

Q2 = detval(:,2);

logm_sem = - dof*log(Q1) + T*Q2 - 0.5*Q3;

adj = max(logm_sem);
results.maxsem = adj;
madj = logm_sem - adj;
xx = exp(madj);

% trapezoid rule integration
yy =rgrid;
isum = sum((yy(2:ngrid,1) + yy(1:ngrid-1,1)).*(xx(2:ngrid,1) - xx(1:ngrid-1,1))/2);

logm_out = isum + adj + logC_sem;

results.logm_sem = logm_out;

% ====================================================================
% evaluate log-marginal for SDEM model over a grid of rho values

dof = (N*T - 2*nx - dofadj)/2;

D = (1 - 1/rmin); % from uniform prior on rho

logC_sdem = -log(D) + gammaln(dof) - dof*log(2*pi); 
% 
% 
% % do vectorized calculations
x = xsdm;
Wxsem=[Wx WWx];

xpx = x'*x;
xpWx = x'*Wxsem;
xpWpx = Wxsem'*x;
xpWpWx = Wxsem'*Wxsem;
xpy = x'*y;
xpWy = x'*Wy;
xpWpy = Wxsem'*y;
xpWpWy = Wxsem'*Wy;
ypy = y'*y;
ypWy = y'*Wy;
ypWpy = Wy'*y;
ypWpWy = Wy'*Wy;

Q1 = zeros(ngrid,1);
Q3 = zeros(ngrid,1);

for iter=1:ngrid;
    rho = rgrid(iter,1);
    
 Axx = xpx - rho*xpWx - rho*xpWpx + rho*rho*xpWpWx;
 
 Q3(iter,1) = log(det(Axx));
 
 Axy = xpy - rho*xpWy - rho*xpWpy + rho*rho*xpWpWy;

 Ayy = ypy - rho*ypWy - rho*ypWpy + rho*rho*ypWpWy;

 b = Axx\Axy;

 Q1(iter,1) = Ayy - b'*Axx*b;

end;

Q2 = detval(:,2);

logm_sdem = - dof*log(Q1) + T*Q2 - 0.5*Q3;

adj = max(logm_sdem);
results.maxsdem = adj;
madj = logm_sdem - adj;
xx = exp(madj);

% trapezoid rule integration
yy =rgrid;
isum = sum((yy(2:ngrid,1) + yy(1:ngrid-1,1)).*(xx(2:ngrid,1) - xx(1:ngrid-1,1))/2);

logm_out = isum + adj + logC_sdem ;

results.logm_sdem = logm_out;

% ====================================================================
% evaluate log-marginal for SAC model over a grid of rho and lambda values

% dof = (N*T -nx -1)/2;
% 
% logC_sac = gammaln(dof) - dof*log(2*pi);
% 
% Q1 = zeros(ngrid,ngrid);
% logm_sac = zeros(ngrid,ngrid);
% 
% x = xsar;
% 
% for i=1:ngrid;
%     rho = rgrid(i,1);
%          
%     detval1 = log(det(speye(N) - rho*W1));
%    
%     for j=1:ngrid;
%         lam = rgrid(j,1);
%                 
%         detval2 = log(det(speye(N) - lam*W2));
% 
%         A = speye(nt) - rho*sparse(W1b); % W1b = kron(eye(T),W1small)
%         B = speye(nt) - lam*sparse(W2b);
%         Bx = (speye(nt) - lam*W2b)*xsar;
%         bhat = (Bx'*Bx)\(Bx'*B*A*y);
%         e = B*(A*y - x*bhat);
%         
%         Q1(i,j) = e'*e;
%         
%         xpx = Bx'*Bx;
%         
%         lndetx_sac = log(det(xpx));
%                                  
%         logm_sac(i,j) = log(bprior(j,1)) + log(bprior(i,1)) -dof*log(Q1(i,j)) + T*detval1 + T*detval2 -0.5*lndetx_sac;
%         
%     end;
% end;
% 
% adj = max(max(logm_sac));
% madj = logm_sac - adj;
% xx = exp(madj);
% 
% isum = zeros(ngrid,1);
% for j=1:ngrid;
% isum = isum + WW*xx(:,j); % integrate over lambda 
% end;
% 
% ijsum = WW*isum;  % integrate over rho
% 
% logm_out = ijsum + adj + logC_sac;
% 
% results.logm_sac = logm_out;

% ===========================================================

time = etime(clock,timet);

results.time = time;

% calculate posterior model probabilities
lmarginal = [results.logm_sar results.logm_sdm results.logm_sem results.logm_sdem ];

adj = max(lmarginal);
madj = lmarginal - adj;

xx = exp(madj);

% compute posterior probabilities
psum = sum(xx);
probs = xx/psum;

results.probs = probs';
results.lmarginal = lmarginal';

function out = beta_prior(rvec,a1,a2)
% PURPOSE: construct beta-prior for rho over -1,1 interval
%-----------------------------------------------------------
% USAGE: out = beta_prior(a1,a2,rvec);
% where:    rvec = grid over rmin,rmax interval, an n x 1 vector
%           a1 = (optional) prior parameter (default = 1.1)
%           a2 = (optional) prior parameter (default = 1.1)
% RETURNS: out = nx1 vector of prior weights for rho
%          in the interval rmin,rmax
%-----------------------------------------------------------
% NOTES: increasing a1,a2 to 1.5,1.5 or 2.0,2.0 increases
%        prior weight placed on zero for rho, and decreases
%        weights on non-zero values for rho
% to see what the prior looks like:
% rvec = -1:0.01:1;
% a1 = 1.1; a2 = 1.1;
% bprior = beta_prior(rvec',a1,a2);
% plot(rvec,bprior);
%-----------------------------------------------------------

% written by:
% James P. LeSage, 4/2003
% University of Toledo
% 2801 W. Bancroft St,
% Toledo, OH 43606
% jlesage@spatial-econometrics.com

if nargin == 1
a1 = 1.01;
a2 = 1.01;
elseif nargin == 2
    a2 = 1.01;
elseif nargin > 3
    error('beta_prior: wrong # of inputs');
end;

B = beta(a1,a2);
num = (1+rvec).^(a1-1);
num = num.*(1-rvec).^(a2-1);
den = 2^(a1+a2-1);
out = (1/B)*num/den;
out(1) = realmin;
out(end) = realmin;

function [a] = logdet(A)

U = chol(A);
a = 2*sum(log(diag(U)));


