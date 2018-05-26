function Residual=MRes(CBV,F,LgLambda,t)
% Mouridsen Residual Function (2006, Neuroimage);
% it is a gamma cumulative distribution function from t>Inf
% typical values: CBV=2; CBF=10; Lglambda=0; t=0:99;

% CBV=2; CBF=10;Lglambda=4.65;t=0:99;

lambda=exp(LgLambda);
CBF=F;
alfa=lambda; beta=CBV./(lambda*CBF);

for i1=1:length(t)
    Residual(i1)=1-gamcdf(t(i1),alfa,beta);
%     Residual(i1)=gammainc((t(i1)/beta),alfa,'upper');
end


% plot(t,Residual);hold on

% for i1=1:length(t)
%    h(i1)=(1/beta^alfa/gamma(alfa))*t(i1)^(alfa-1)*exp(-t(i1)/beta);
% end
% 
% for i2=1:length(t)
% %     Sigma=0;i3=1;   
%     Sigma=sum(h(1:i2-1));
%     Residual(i2)=1-Sigma;
% end