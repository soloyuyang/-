function llike = f_sacpanel(parm,y,x,W,lambda,T)
rho = parm(1,1);
lam = parm(2,1);
nt = length(y);
n = length(W);
detm=0;
for i=1:n
    detm=detm+log(1-rho*lambda(i))+log(1-lam*lambda(i));
end
Ay=y;
BAy=y;
Bx=x;
for t=1:T
    t1=1+(t-1)*n;t2=t*n;
    Ay(t1:t2,1) = (eye(n) - rho*W)*y(t1:t2,1);
    BAy(t1:t2,1)= (eye(n) - lam*W)*Ay(t1:t2,1);
    Bx(t1:t2,:) = (eye(n) - lam*W)*x(t1:t2,:);
end
b = (Bx)\(BAy);
e = BAy - Bx*b;
epe = e'*e;
llike = (nt/2)*log(epe/nt) - T*detm;