function llike = f2_sacpanel(parm,y,x,W,lambda,T)
k = length(parm); 
n = length(W);
b = parm(1:k-3,1);
rho = parm(k-2,1); 
lam = parm(k-1,1); 
sige = parm(k,1);
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
e = BAy - Bx*b;
epe = e'*e; 
tmp2 = 1/(2*sige);
llike = -(n*T/2)*log(2*pi*sige) + T*detm - tmp2*epe;