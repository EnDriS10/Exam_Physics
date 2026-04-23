%% C3_2024_3_vG.m FUNCIONES USUARIO AL FINAL (duración Control = 4 horas)
clear ; clc ; figure(303) ;
load('Control3_datos3.mat') ; % cargar datos

%% C3_A [FUNCIÓN DE USUARIO] ver abajo

%% C3_B [REPRESENTACIÓN GRÁFICA] 
alpha=2.5 ; t=linspace(0,30,1E3) ; beta=linspace(0.001,0.8,100)' ; % rangos t y beta
F=dist_C3(t,alpha,beta) ;
subplot(3,1,1) ; surf(t,beta,F) ; shading interp ; 
title('distribución f_g') ;
xlabel('tiempo (meses)') ; ylabel('\beta (mes^{-1})') ; zlabel('f_g(t;\alpha=2.5,\beta)') ;
subplot(3,2,3) ; plot(t,F(20,:),t,F(70,:),t,F(99,:)) ;
xlabel('tiempo (meses)') ; ylabel('f_g(t;\alpha=2.5,\beta)') ; grid on ;
legend('\beta=0.15 meses^{-1}','\beta=0.56 meses^{-1}','\beta=0.79 mes^{-1}') ;
title('f_g vs t') ;
subplot(3,2,4) ; plot(beta,F(:,500),beta,F(:,850),beta,F(:,999)) ;
xlabel('\beta (meses^{-1})') ; ylabel('f_g(t;\alpha=2.5,\beta)') ; grid on ;
legend('t=15 meses','t=25.5 meses','t=29.97 meses') ; title('f_g vs \beta') ;

%% C3_3C NORMALIZACIÓN 
fg=dist_C3(t,3.2,0.6) ;
area1=sum(diff(t).*(fg(1:end-1)+fg(2:end))/2) ; % integral trapezoidal
maxy=max(fg) ; maxx=max(t) ;
N=1E5 ; x=maxx*rand(1,N) ; y=maxy*rand(1,N) ; % Montecarlo
p=sum(y<dist_C3(x,3.2,0.6))/N ;
area2=maxx*maxy*p ;
fprintf('El área bajo la curva por integración numérica es %.3f\n',area1) ;
fprintf('El área bajo la curva por integración MonteCarlo es %.3f\n',area2) ;

%% C3_3D CARACTERIZACIÓN ESTADÍSTICA 
N=1E5 ; [counts,centers]=hist(muestra,[min(muestra):7/30:max(muestra)]) ; %centers=centers*30/7 ;
subplot(3,1,3) ; bar(centers,counts/N) ; grid on ;
title('distribución de probabilidad con \it{bins} de una semana') ;
xlabel('duración desempleo (meses)') ; ylabel('probabilidad') ;
media=mean(muestra) ; incert=std(muestra) ;
fprintf('i) La duración media de desempleo es %.2f meses y su incertidumbre es %.2f meses\n',media,incert) ;
p18=sum(muestra>18)/N ; 
fprintf('ii) La probabilidad de tiempo de desempleo superior a 18 meses es %.3f\n',p18) ;
p3a15=sum((muestra>3)&(muestra<15))/N ;
fprintf('iii) La probabilidad de estar desempleado entre 3 y 15 meses es %.2f\n',p3a15) ;

cumP=cumsum(counts/N) ; % probabilidad cumulativa
k=find(cumP>=0.9) ; tiempo90=centers(k(1)) ;
fprintf('iv) El tiempo para que (estadísticamente) el 90%% encuentren empleo es %.1f meses\n',tiempo90) ; 

%% FUNCIÓN USUARIO ================================================================
function [fg]=dist_C3(t,alpha,beta)
    t=t(1:end) ; % t son columnas
    beta=beta(1:end) ; % beta son filas
    fg=beta.^alpha.*t.^(alpha-1).*exp(-beta*t)/gamma(alpha) ;
end
