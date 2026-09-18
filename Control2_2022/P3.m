clear;clc;

a=2;
b=2;
R=2;

L=4; 

N=1E7;
x=(2*rand(1,N)-1)*R*a ; y=(2*rand(1,N)-1)*R ; z=(2*rand(1,N)-1)*R; % genera puntos x y en el cubo 
Cond1=((x.^2)/(a^2) + (y.^2)/(b^2))<z ;  % encontrar elementos que están dentro del paraboloide
Cond2=(y.^2+z.^2)<R^2;
Cond3 = Cond1 & Cond2;
p=sum(Cond3)/N ;  

V=p*(L*L*2*R*(a))