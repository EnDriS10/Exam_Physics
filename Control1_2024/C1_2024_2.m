%% C1_2024_2.m [3.5 puntos]
clear ; clc ;

%% 2A [0.5 puntos] generar matrices A y B.
N=100 ; A=zeros(N,N) ; B=zeros(N,N) ;
for n=1:N 
    for m=1:N ; A(n,m)=rem(n,m) ; b=n*m-(n+m) ; B(n,m)=(-1)^b ; end
end

%% 2B [1.50 puntos] matrices C, D y E. Suma diagonal secundaria.
n5=N/5 ; n3=fix(N/3) ; % el fix() es un poco excesivo
D=zeros(n5,n3) ; E=zeros(n3,n5) ; 
for n=1:n5 
    for m=1:n3 ; D(n,m)=A(5*n,3*m) ; end
end
for n=1:n3 
    for m=1:n5 ; E(n,m)=B(3*n,5*m) ; end
end
C=D*E ; % producto matricial
suma=0 ; for n=1:n5 ; suma=suma+C(n,n5-n+1) ; end
%suma=trace(fliplr(C)) ; % código alternativo
fprintf('Suma de elementos de la diagonal secundaria de la matriz C = %.0f\n',suma) ;

%% 2C [1.50 puntos] matriz F y número valores <=10.
F=zeros(N,N) ;
for n=1:N
    for m=1:N
        nm=n+m ;
        if rem(nm,2)~=0 ; F(n,m)=A(n,m) ; % impar            
        elseif (rem(nm,2)==0)&&(~(rem(nm,4)==0)) ; F(n,m)=B(n,m) ;
        elseif rem(nm,4)==0 ; F(n,m)=nm^2 ;
        else ; fprintf('NO debería suceder!\n') ;
        end
    end
end

b=(F<=10) ; nn=sum(b(:)) ; % nn=sum(sum(b)); Acceso mediante matriz lógica y suma de verdaderos
%nn=length(find(F<=10)) ; % alternativamente con find()
fprintf('Número de valores menores o iguales a 10 en la matriz F = %d\n',nn) ;
