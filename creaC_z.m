function [C_z]=creaC_z(N,alfa,beta)


n=N-alfa;

t=exp(-beta*[1:n]);
t=fliplr(t);

for i=1:n;
    B(i,i:n)=[(t(i).^2)/2].*(t(i:n)-t(i)/3);
    B(i:n,i)=B(i,i:n)';
end

B=rot90(B);B=rot90(B);
A=(1e-12)*eye(alfa);

C_z=[A zeros(alfa,N-alfa);zeros(n,alfa) B];