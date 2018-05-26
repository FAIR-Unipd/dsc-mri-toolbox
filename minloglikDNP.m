function [l]=minloglikDNP(par,y,G,C_v,alfa,sel)

N=size(G,1);
Nl=size(G,2);

l2=exp(par(1));beta=exp(par(2));o2=1;%exp(par(3));

C=G*exp([zeros(1,alfa) -beta*[1:(Nl-alfa)]])';
[C_z]=l2*creaC_z(Nl,alfa,beta);
E=G*C_z*G';
M=E+o2*C_v;
Mi=inv(M);

A=Mi*(eye(N)-(sel==1)*C*inv(C'*Mi*C)*C'*Mi);

try
    Ec=real(chol(M));
    e=abs(diag(Ec).^2);
    b=sum(log(e))+(sel==1)*log(det(C'*Mi*C));
    l=0.5*b+0.5*y*A*y';
catch
    l=+inf;
end
end