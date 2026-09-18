clear; clc;

    V_0= 1;
    d_0= 0.2;


    %el if no funciona cuando hay vectores
% 
% function Eint = Potencial_C1(d)
%     V_0 = 1;
%     d_0 = 0.2;
% 
%     d_1 = d_0*(1+sqrt(3));
%     a = V_0/(d_0^2);
%     b = (-2*a)/(3*(d_1-d_0));
% 
%     Eint = zeros(size(d)); % inicializa el vector de salida
% 
%     % Zona 1
%     mask1 = d <= d_0;
%     Eint(mask1) = a .* (d(mask1).^2 - 2*d_0.*d(mask1));
% 
%     % Zona 2
%     mask2 = (d > d_0) & (d <= d_1);
%     Eint(mask2) = a .* (d(mask2).^2 - 2*d_0.*d(mask2)) + b .* (d(mask2) - d_0).^3;
% end

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

d = linspace(0, 3*d_0, 100);
Eint= Potencial_C1 (d);

subplot(2,1,1) ; plot(d,Eint) ; grid on ;
xlabel('distancia (nm)') ; ylabel('energía potencial (eV)') ;
title('Energía potencial vs d') ;
%%
function [D]=Distancia(V)
    [~,n]=size(V) ;
    D=zeros(n,n) ;
    for i=1:n 
        for j=1:n 
            D(i,j)=sqrt(sum(((V(:,i)-V(:,j)).^2),1)) ;            
        end
    end
end


p0=1/(10*sqrt(2)) ;
p=p0*transpose([1,1,1;-1,-1,1;1,-1,-1;-1,1,-1]) ;


for i= 1:4
    for j= 1:4
        D= Distancia(p);
        E1(j)= Potencial_C1 (D(i,j));
        E2(i)= sum(E1);
        end
end

E_tet=0.5*sum(E2);
fprintf('La energía potencial total del tetraedro es %.3f eV\n',E_tet) ;

s=linspace(10*p0,0,1000) ; % distancia al origen (nm)
EEp=zeros(1,1000) ;
for i=1:1000        
    DD=Distancia([p,[s(i),0,0]']) ;
    P=Potencial_C1(DD); 
    EEp(i)=sum(P(5,:)) ;
end
EEp=EEp+E_tet ;
subplot(2,1,2) ; plot(s,EEp) ; grid on ;
xlabel('distancia al origen (nm)') ; ylabel('energía potencial (eV)') ;
title('Energía potencial con partícula vs s') ;