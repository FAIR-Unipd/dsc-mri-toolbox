function SIG=mour_Psig_v1(input,CaCBVTime)
% input=[LgFlow,LgLambda,delay]; Ca->AIF+CBV+time

log_F=input(1); LgLambda=input(2); delay=input(3); 
estCBV=CaCBVTime(1,end);
t=CaCBVTime(2,1:length(CaCBVTime)-1);
Ca=CaCBVTime(1,1:length(CaCBVTime)-1);

Sig0=CaCBVTime(3,1); k=CaCBVTime(3,2); TE=CaCBVTime(3,3);
% F(F<0)=0; 
F=exp(log_F); F(F>1e4)=1e4; % max limit of F is 10^6
delay(delay>t(end))=t(end);

t_dCa=t+delay;
dCa=interp1(t_dCa,Ca,t);  dCa(isnan(dCa))=0;

Residual=MRes(estCBV,F,LgLambda,t);
CtconvRes=conv(dCa,Residual);
CTC=F*CtconvRes(1:length(t));
CTC( (k.*CTC.*TE)>100 ) = 100;
CTC( (k.*CTC.*TE)<-100 ) = -100;
SIG=Sig0.*exp(-(k.*CTC.*TE));