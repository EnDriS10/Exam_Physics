clear;clc;

N=8;
K=1;

n= 1:N;
j=n;

w=2*pi/N;
xn = 2*cos((n-1)*w + pi/2);
yn = 2*sin((n-1)*w+ pi/2) + 3;
zn = zeros(1,N);

Xj =(j-1) - N/2;
Yj= zeros(1,N); Zj= zeros(1,N);

V_0=5* ones(N,1);

A=zeros(N,N);
for i= 1:N
    for j= 1:N

        x1=Xj(i);
        y1=Yj(i);
        z1=Zj(i);

        x2=xn(j);
        y2=yn(j);
        z2=zn(j);


        A(i,j)=K/Dist(x1,y1,z1,x2,y2,z2);
end
end


qn=inv(A')*V_0;



function [Dij] = Dist(x,y,z,x2,y2,z2)
    Dij=sqrt((x-x2)^2 + (y-y2)^2 + (z-z2)^2);
end

 
 subplot(2,1,1) ;
 hold on
 scatter(Xj,Yj,'b') ;
 scatter(xn,yn,'r') ;
 grid on
 axis equal


 subplot(2,1,2) ;
 scatter(Xj,Yj,'b') ;
  axis equal ;    
 % % sólo el subplot 5
  grid on ; % mostrar retícula
  hold off