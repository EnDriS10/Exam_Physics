%%  C2_2023_2_Cargas.m
clear ; clc ; figure(202) ;

%% DATOS
N=17 ; w=2*pi/N ; V0=-3 ; K=8.99E9 ; % unidades SI

%% POSICIONES de CARGAS y de POTENCIAL
[x,y,z]=deal(zeros(N,1)) ; % pre-alojar espacio
n=1:N ; % etiquetar posiciones por índice
x=5*cos((n-1)*w+pi/2) ; y=2*sin((n-1)*w+pi/2)+3 ; % posiciones cargas
[X,Y,Z]=deal(zeros(N,1)) ;
X=(n-1)-N/2 ; % posiciones potencial
subplot(2,1,1) ;
plot(x,y,'or',X,Y,'ob',x(1),y(1),'ok') ; axis equal ; grid on ;
xlabel('posición x (m)') ; ylabel('posición y (m)') ;
legend('cargas','posiciones V_0') ; title('posiciones de cargas y potenciales') ;

%% CARGAS 
d=zeros(N,N) ;
for i=1:N % posiciones
    for j=1:N % cargas
       d(i,j)=sqrt((x(j)-X(i))^2+(y(j)-Y(i))^2) ; 
    end
end
U=K./d ; % matriz de coeficientes
V=ones(N,1)*V0 ; % columna de términos independientes
q=inv(U)*V ; % resolver
subplot(2,1,2) ;
uC=1E-6 ; % unidades micro-Coulombios
plot(n,q/uC,'o-') ; grid on ; title('valor de las cargas') ;
xlabel('número de carga') ; ylabel('carga (\muC)') ;
