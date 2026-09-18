clear;clc;


N=1E6;
R=2;

xn = R*(2*rand(1,N)-1);
yn = R*(2*rand(1,N)-1);
zn = R*(2*rand(1,N)-1);
sn = R*(2*rand(1,N)-1);

p= sum(xn.^2 + yn.^2 +zn.^2 + sn.^2 < R^2 )/N;

V=p*R^4
