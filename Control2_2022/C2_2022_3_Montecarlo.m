%% C2_2022_3_Montecarlo.m
clear ; clc ;

%% DATOS
N=1E7 ; R=2 ; a=2 ; b=2 ; 

%% GENERAR COORDENADAS XYZ aleatorias
% paralelepípedo adecuado por inspección visual/analítica
beta=2 ; alpha=2*sqrt(a) ; %=2.82842712475 
x=alpha*(2*rand(N,1)-1) ; 
y=beta*(2*rand(N,1)-1) ;
z=beta*(1*rand(N,1)) ;

%%
c=(y.^2+z.^2)<R^2 ; % esfera
p=((x.^2)/(a^2)+(y.^2)/(b^2))<z ; % paraboloide

cp=sum(c&p)/N ;
% fprintf('Probabilidad %.3f\n',cp) ; % no se pide
V=(2*alpha)*(2*beta)*(1*beta) ;
fprintf('El volumen común encerrado es %.3f metros cúbicos\n',cp*V) ;
