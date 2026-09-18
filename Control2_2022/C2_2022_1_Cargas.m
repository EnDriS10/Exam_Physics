%%  C2_2022_1_Cargas.m
clear ; clc ; figure(201) ;

%% DATOS
N=8 ; w=2*pi/N ; V0=5 ; K=9E9 ; % unidades SI

[x,y,z]=deal(zeros(N,1)) ; % pre-alojar espacio
n=1:N ; % etiquetar cargas y posiciones por índice
x=2*cos((n-1)*w+pi/2) ; y=2*sin((n-1)*w+pi/2)+3 ; % posiciones cargas
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
U=K./d ; % matriz coeficientes
V=ones(N,1)*V0 ; % columna términos independientes
q=inv(U)*V ; % resolver
uC=1E-6 ; % unidades micro-Coulombios
subplot(2,1,2) ; plot(n,q/uC,'o-') ; grid on ; title('valores de cargas') ;
xlabel('número de carga') ; ylabel('carga (\muC)') ;
