%% C2_2024_3_Montecarlo.m
clear ; clc ;

%% PARALELEPÍPEDO (ajustar tamaño)
N=1E7 ; R=2.0 ;
% generar valores aleatorios
x=(2*rand(N,1)-1) ; 
y=R*(2*rand(N,1)-1) ;
z=R*(2*rand(N,1)-1) ;

%% MONTECARLO
b1=((x-1).^2+y.^2+z.^2)<R^2 ; % dentro de esfera 1
b2=((x+1).^2+y.^2+z.^2)<R^2 ; % dentro de esfera 2

pp=sum(b1&b2)/N ; % encontrar la fracción de los que verifican la intersección
%fprintf('Probabilidad %.3f\n',pp) ; % no se pide
V=(2*R)*(2*R)*(2) ; % volumen del paralelepípedo
fprintf('El volumen común encerrado es %.3f metros cúbicos\n',pp*V) ;
