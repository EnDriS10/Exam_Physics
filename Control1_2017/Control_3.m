clear; clc;


function [E,V] =CampoPotencial (q,rq,r)
    epsilon= 8.8542E-12;
    k=(1/4*pi*epsilon);

    r_=r-rq;

    E = k*q.*r_./((norm(r_))^3);
    V = k*q/(norm(r_));
    
end

dx=0.01;

x= -6.5:dx:1.5;

rq=[0,0,0];
q=1;


for i= 1:length(x);
    r =[x(i);0;0];
    [E,V] = CampoPotencial (q,rq,r);
    for j= 1:3;    
        EE(1,i,j)=E(j);
        VV(1,i)=V;
    end 
end


figure(1)
subplot(1,2,1) ;
plot(x, EE(1,:,1)) ;
xlabel('Desplazamiento en X (m)') ;
ylabel('Campo Electrico (N/C)') ;
title('Campo Electrico a lo largo de una linea') ;
%axis equal;
grid on;

subplot(1,2,2) ;
plot(x, VV) ;
xlabel('Desplazamiento en X (m)') ;
ylabel('Potencial Electrico (N/(C*m))') ;
title('Potencial Electrico a lo largo de una linea') ;
%axis equal;
grid on;