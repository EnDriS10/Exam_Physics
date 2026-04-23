%% C3_2023_3_vG (duración Control = 4 horas)
clear ; clc ; figure(303) ;

%% DATOS
N=1E4 ; tau=1.5E3 ;

%% 3A GENERAR...
t=zeros(1,N) ; for i=1:N ; t(i)=-tau*log(1-rand) ; end
%load('DatosExp_C3_3.mat') ; t=t_desintegracion ;

%% ANÁLISIS...
[counts,tcenter]=hist(t,1E4) ; 
FF=cumsum(counts) ; kk=find(FF>=N/2) ; k=kk(1) ; % basado en histograma cumulativo  
fprintf('tiempo mitad de núcleos desintegrados = %.1f segundos\n',tcenter(k)) ; 
p=sum(t>tau)/N ; fprintf('probabilidad desintegración tiempo>tau = %.3f\n',p) ;
pp=sum((t>=tau/3)&(t<=5*tau/3))/N ;
fprintf('probabilidad desintegración entre tiempo>=tau/3 y tiempo<=5*tau/3 = %.3f\n',pp)

%% Calcular dt, H y gráfica...
tsorted=sort(t) ; dt=diff(tsorted) ;
n=N-2:-1:0 ; H=n.*dt/tau ;
[dtCounts,dtCenter]=hist(H,1E2) ;
Hmean=mean(H) ; Hsigma=std(H,1) ;
bar(dtCenter,dtCounts/N) ;
maximo=max(dtCounts/N)*1.2 ;
hold on ; plot([Hmean,Hmean],[0,maximo],'r','LineWidth',2) ; hold off ;
hold on ; plot([Hmean+Hsigma,Hmean+Hsigma],[0,maximo],'g--','LineWidth',2) ; hold off ;
hold on ; plot([Hmean-Hsigma,Hmean-Hsigma],[0,maximo],'g','LineWidth',2) ; hold off ;
legend('probabilidad H','<H>','<H>+\sigma_H','<H>-\sigma_H') ; grid on ;
xlabel('valor H=n\Deltat/\tau (adimensional)') ; ylabel('probabilidad H') ;
