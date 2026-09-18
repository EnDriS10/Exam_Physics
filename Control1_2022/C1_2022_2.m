%% C1_2022_2.m (4.25 puntos)
clear all ; clc ; figure(172) ;

%% 2A Función Potencial_C1, ver abajo (1.50 puntos)

%% 2B Potencial (0.50 puntos)
N=1E3 ; d0=0.2 ; d=linspace(0,3*d0,N) ;
E=Potencial_C1(d) ;
subplot(2,1,1) ; plot(d,E) ; grid on ;
xlabel('distancia (nm)') ; ylabel('energía potencial (eV)') ;
title('Energía potencial vs d') ;

%% 2C Energía tetraedro (0.75 puntos)
p0=1/(10*sqrt(2)) ;
p=p0*transpose([1,1,1;-1,-1,1;1,-1,-1;-1,1,-1]) ;
EE=Potencial_C1(Distancia(p)) ;
Etot=0.5*sum(sum(EE)) ;
fprintf('La energía potencial total del tetraedro es %.3f eV\n',Etot) ;

%% 2D energía con partícula (1.50 puntos)
s=linspace(10*p0,0,N) ; % distancia al origen (nm)
EEp=zeros(1,N) ;
for i=1:N        
    D=Distancia([p,[s(i),0,0]']) ;
    P=Potencial_C1(D) ;
    EEp(i)=sum(P(5,:)) ;
end
EEp=EEp+Etot ;
subplot(2,1,2) ; plot(s,EEp) ; grid on ;
xlabel('distancia al origen (nm)') ; ylabel('energía potencial (eV)') ;
title('Energía potencial con partícula vs s') ;


%% FUNCIONES =======================================================
% función intermedia, ver Potencial_C1 más abajo
function [Eint]=Potencial(d)
    d0=0.2 ; V0=1.0 ; d1=d0*(1+sqrt(3)) ; a=V0/d0^2 ; b=-2*a/(3*(d1-d0)) ;
    if d<=d0 Eint=a*(d^2-2*d0*d) ;
    elseif (d>d0)&&(d<=d1) Eint=a*(d^2-2*d0*d)+b*(d-d0)^3 ;
    else Eint=0 ;
    end
end

function [Eint]=Potencial_C1(D)
    [n,m]=size(D) ;
    Eint=zeros(n,m) ;
    for i=1:n for j=1:m Eint(i,j)=Potencial(D(i,j)) ; end ; end
end

function [D]=Distancia(V)
    [~,n]=size(V) ;
    D=zeros(n,n) ;
    for i=1:n 
        for j=1:n 
            D(i,j)=sqrt(sum(((V(:,i)-V(:,j)).^2),1)) ;            
        end
    end
end