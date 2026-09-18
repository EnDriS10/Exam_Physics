
%% Generador de num aleatorios

% rand: Números reales con dominio de 0 a 1. DP constante.
% randn: Números realescon dominio de –Inf a +Inf. DP gaussiana o normal.
% randi: Números naturales con dominio de 1 a n (n es un parámetro). DP constante


% VA con distribución de probabilidad constante de números naturales. Ejemplo:
% R=4*randi(10,1,N) ; % números múltiplos de 4 hasta 40
% VA con distribución de probabilidad gaussiana. Ejemplo:
% R=5*randn(1,N)+18 ; % gaussiana con media 18 y desviación estándar 5

%% METODO DE MONTECARLO

% % PARALELEPÍPEDO (ajustar tamaño)
% N=1E7 ; R=2.0 ;
% % generar valores aleatorios
% x=(2*rand(N,1)-1) ; 
% y=R*(2*rand(N,1)-1) ;
% z=R*(2*rand(N,1)-1) ;

% % MONTECARLO
% b1=((x-1).^2+y.^2+z.^2)<R^2 ; % dentro de esfera 1
% b2=((x+1).^2+y.^2+z.^2)<R^2 ; % dentro de esfera 2
% 
% pp=sum(b1&b2)/N ; % encontrar la fracción de los que verifican la intersección
% %fprintf('Probabilidad %.3f\n',pp) ; % no se pide
% V=(2*R)*(2*R)*(2) ; % volumen del paralelepípedo
% fprintf('El volumen común encerrado es %.3f metros cúbicos\n',pp*V) ;


%% MOVIMIENTO BROWNEANO

% d=2 ; N=100 ; % desplazamiento máximo en cada paso y número de pasos
% function [x,y,R]=FunCaminoXY(d,N)
% dx=(2*rand(1,N)-1)*d ; dy=(2*rand(1,N)-1)*d ; % valores de los pasos
% x=[0,cumsum(dx)] ; y=[0,cumsum(dy)]  ; % trayectoria acumulada
% R=sqrt(x(end)^2+y(end)^2) ; % distancia final al origen
% end

%% Genera El dominio con la Distribucion de P dada
% N=1E3 ; % generar N valores
% Dominio=[2,4,7,8] ;
% DistribProb=[0.2,0.1,0.5,0.2] ;
% r=rand(1,N) ; R=zeros(1,N) ; % VA pre-alojar espacio
% ProbCum=cumsum(DistribProb) ; % probabilidad cumulativa
% 
% for i=1:N
%     if r(i)<=ProbCum(1) R(i)=Dominio(1) ; 
%     elseif r(i)<=ProbCum(2) R(i)=Dominio(2) ;
%     elseif r(i)<=ProbCum(3) R(i)=Dominio(3) ; 
%     else
%     R(i)=Dominio(4) ; 
%     end
% end
% subplot(2,1,1) ; hist(R,100);    % etiquetar ejes
% subplot(2,1,2) ; plot(R,'o') ;   % etiquetar ejes

%% Generalizar el algotimo anterior

%Ejemplo

% Dom=1:200 ; % dominio de la variable aleatoria
% 
% Prob=5.2841e-04*Dom.^0.5 ; % distribución de probabilidad del dominio
% R=RandList(1E5,Dom,Prob) ; % función de usuario, ver al final
% hist(R,20) % histograma
% xlabel('Variable aleatoria X (Voltios)') ;
% ylabel('cuentas') ;
% title('histograma (10^5 valores)') ;


%%
% [v,x]=hist(datos,[5:15]) ; % v es la frecuencia y x es la posicion
%bar(x,v/N) ; xlabel('valor suma (número entero)') ; ylabel('probabilidad') ;

% %% ANÁLISIS...
% [counts,tcenter]=hist(t,1E4) ; 
% FF=cumsum(counts) ; kk=find(FF>=N/2) ; k=kk(1) ; % basado en histograma cumulativo  
% fprintf('tiempo mitad de núcleos desintegrados = %.1f segundos\n',tcenter(k)) ;

% p=sum(t>tau)/N ; fprintf('probabilidad desintegración tiempo>tau = %.3f\n',p) ;
% pp=sum((t>=tau/3)&(t<=5*tau/3))/N ;

% fprintf('probabilidad desintegración entre tiempo>=tau/3 y tiempo<=5*tau/3 = %.3f\n',pp)

% %% Calcular dt, H y gráfica...
% tsorted=sort(t) ; dt=diff(tsorted) ;
% n=N-2:-1:0 ; H=n.*dt/tau ;
% [dtCounts,dtCenter]=hist(H,1E2) ;
% Hmean=mean(H) ; Hsigma=std(H,1) ;
% bar(dtCenter,dtCounts/N) ;
% maximo=max(dtCounts/N)*1.2 ;
% hold on ; plot([Hmean,Hmean],[0,maximo],'r','LineWidth',2) ; hold off ;
% hold on ; plot([Hmean+Hsigma,Hmean+Hsigma],[0,maximo],'g--','LineWidth',2) ; hold off ;
% hold on ; plot([Hmean-Hsigma,Hmean-Hsigma],[0,maximo],'g','LineWidth',2) ; hold off ;
% legend('probabilidad H','<H>','<H>+\sigma_H','<H>-\sigma_H') ; grid on ;
% xlabel('valor H=n\Deltat/\tau (adimensional)') ; ylabel('probabilidad H') ;

%% FUNCIONES DE USUARIO ==================================================
function [R]= RandList(N,Dom,Prob)
    ProbCum=cumsum(Prob) ; % probabilidad cumulativa
    Ndom=length(Dom) ; % número de elementos del dominio
    R=zeros(N,1) ; % pre-alojar espacio
    for i=1:N
        r=rand ; 
        for j=1:Ndom           
            if r<ProbCum(j) 
                R(i)=Dom(j) ; 
                break ; % asignado -> salir del bucle for 
            else 
                R(i)=Dom(end) ;
            end
        end
    end 
end

function [Dom , Frec] = Histogr(Arr)
    
    Dom=zeros(1,length(Arr)+1);
    Dom(1)=Arr(1);
    Frec=zeros(1,length(Arr)+1);
    for i = 1:length(Arr)
        if Dom(i) == Arr(i)
        Frec(i) = Frec(i) +1;
        else
            Dom(i)= Arr(i);
            Frec(i)= 1;
        end
    end   
end


%% ================================
 
bar(x) % grafica valores le asigna al primer elemento del array 1 y asi.
bar(x,y,'r')              % x = posiciones, y = altura


bar(y,'FaceColor','g')
bar(y,'EdgeColor','k')
bar(y,'LineWidth',2)
bar(y,'BarWidth',0.5)

bar(y,'stacked')     % barras apiladas
bar(y,'grouped')     % agrupadas


%% ================================
%media

mean(x)
mean(x, 1) %Es el valor por defecto (media de las columnas).
mean(x, 2) %Calcula la media de cada fila. El resultado será un vector columna.
mean(x,'all') % Media de todo los elemntos del array

%VAR (varianza)

var(x)
var(x,1)     % normalizacion por N
var(x,0)     % normalizacion por N-1
var(x,[],2)


%STD (desviacion estandar)

std(x)
std(x,1)
std(x,0)

var(x, [], 1)% (Por defecto) Calcula la varianza de cada columna.
var(x, [], 2)% Calcula la varianza de cada fila.

%MODE (moda)

mode(x)
mode(x, 1) %Es el valor por defecto (media de las columnas).
mode(x, 2) %Calcula la media de cada fila. El resultado será un vector columna.


%% ================================
%% MAX
%% ================================

max(x)
max(x,y)     % max elemento a elemento
[max_val, index] = max(x)
max(x,[],2)


%% ================================
%% MIN
%% ================================

min(x)
min(x,y)
[min_val, index] = min(x)
min(x,[],2)

%% POISSON

M=10*10 ; % celdas
p=zeros(1,M) ; % n partículas en cada celda
N=350 ; % n partículas lanzadas
% LANZAMIENTO
for i=1:N
n_celda=randi(M) ; % ATENCIÓN a randi
p(n_celda)=p(n_celda)+1 ;
end
x=[0:15] ; % bins del histograma: números enteros
[Nz,z_mid]=hist(p,x) ;
subplot(2,1,1) ; plot(p)
subplot(2,1,2) ; bar(z_mid,Nz/sum(Nz)) ; 
xlim([-1 15]) ; % N/M es la media
