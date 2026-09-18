%%  C2_2024_2_Cargas.m (FUNCIONE DE USUARIO AL FINAL)
clear ; clc ; figure(202) ;

%% DATOS 
N=17 ; V0=-3 ; K=8.99E9 ; % unidades SI

%% ESTABLECER PSOCISIONES CARGAS y V0
[x,y,z]=deal(zeros(N,1)) ; % pre-alojar espacio
n=1:N ; % etiquetar cargas y posiciones por índice
x=3*(n-9.0) ; y=ones(1,N)*10 ; % posiciones de las cargas
[X,Y,Z]=deal(zeros(N,1)) ; % posiciones en las que se conoce el potencial
X=(n-1)-N/2 ; 
subplot(1,2,1) ;
plot(x,y,'or',X,Y,'ob',x(1),y(1),'ok') ; axis equal ; grid on ;
xlabel('posición x (m)') ; ylabel('posición y (m)') ;
legend('posiciones cargas','posiciones V_0') ;

%% CARGAS 
d=zeros(N,N) ; % pre-alojar espacio
% matriz de distancias
for i=1:N % posiciones
    for j=1:N % cargas
       d(i,j)=sqrt((x(j)-X(i))^2+(y(j)-Y(i))^2) ; 
    end
end
U=K./d ; % matriz de coeficientes
V=ones(N,1)*V0 ; % columna de términos independientes
q=inv(U)*V ; % resolver
uC=1E-6 ; % unidades micro-Coulombios
subplot(1,2,2) ; plot(n,q/uC,'o-') ; grid on ; title('valor de las cargas') ;
xlabel('número de carga') ; ylabel('carga (\muC)') ;
