%% C2_2023_3_Montecarlo.m
clear ; clc ;

%% DATOS (ajustar zona apropiada del paralelepípedo)
N=1E7 ; R=2.3 ; a=3.0 ; b=3.0 ; c=4.0 ;

%% GENERAR COORDENADAS XYZ aleatorias
% paralelepípedo adecuado por inspección visual/analítica
x=a*(2*rand(N,1)-1) ; 
y=R*(2*rand(N,1)-1) ;
z=R*(2*rand(N,1)-1) ;

%% VOLUMEN COMÚN
bc=(y.^2+z.^2)<R^2 ; % cilindro
be=((x.^2)/(a^2)+(y.^2)/(b^2)+(z.^2)/(c^2))<1 ; % elipsoide

pp=sum(bc&be)/N ; % intersección
%fprintf('Probabilidad %.3f\n',pp) ; % no se pide
V=(2*a)*(2*R)*(2*R) ; % volumen paralelepípedo
fprintf('El volumen común encerrado es %.3f metros cúbicos\n',pp*V) ;
