clear;clc;

N=1E3;

a=3;
b=3;
c=4;
R=2.3;

xn = a*(2*rand(1,N)-1);
yn = R*(2*rand(1,N)-1);
zn = R*(2*rand(1,N)-1);


p1= (xn./a).^2 + (yn./b).^2 + (zn./c).^2 < 1;
p2= yn.^2 +zn.^2 < R^2;
p = p1 & p2;
Prob= sum(p)/N ;

V= Prob* (8*a*R*R)


