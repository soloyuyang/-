load data_cigarette_demand.mat; % stored under the name A
load W_cigarette_demand.mat; % note: stored under the name W1
% Dataset downloaded from www.wiley.co.uk/baltagi/
% Spatial weights matrix constructed by Elhorst
%
% written by: J.Paul Elhorst summer 2010
% University of Groningen
% Department of Economics
% 9700AV Groningen
% the Netherlands
% j.p.elhorst@rug.nl
%
% REFERENCES: 
% Elhorst JP (2010) Matlab Software for Spatial Panels. Under review.
%
% Elhorst JP (2010) Spatial Panel Data Models. In Fischer MM, Getis A (Eds.) 
% Handbook of Applied Spatial Analysis, Ch. C.2. Springer: Berlin Heidelberg New York.
%
% dimensions of the problem
T=30; % number of time periods
N=46; % number of regions
% row-normalize W
W=normw(W1); % function of LeSage
y=A(:,[3]); % column number in the data matrix that corresponds to the dependent variable
x=A(:,[4,6]); % column numbers in the data matrix that correspond to the independent variables
xconstant=ones(N*T,1);
[nobs K]=size(x);
% -------------------------------------------------------------------------
% ols estimation 
results=ols(y,[xconstant x]);
vnames=strvcat('logcit','intercept','logp','logy');
prt_reg(results,vnames,1);
% -------------------------------------------------------------------------
% Use sem_panel_FE(...) to estimate spatial error model
% Use prt_sp(...) to print out coefficient estimates
% -------------------------------------------------------------------------
% Use sar_panel_FE(...) to estimate spatial lag model
% Use prt_sp(...) to print out coefficient estimates
% Set spat_model and use direct_indirect_effects_estimates(...) to print out effects estimates
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
for t=1:T
    t1=(t-1)*N+1;t2=t*N;
    wx(t1:t2,:)=W*x(t1:t2,:);
end
% Use sar_panel_FE(...) to estimate spatial Durbin model
% Change vnames and use prt_sp(...) to print out coefficient estimates
% Set spat_model and use direct_indirect_effects_estimates(...) to print out effects estimates