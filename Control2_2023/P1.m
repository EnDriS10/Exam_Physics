clear; clc;

m=2;
dt=1E-4;
tf=50;

t=linspace(0,tf,1/dt);
x=zeros(1,length(t));
y=zeros(1,length(t));

vx=zeros(1,length(t));
vy=zeros(1,length(t));

ax=zeros(1,length(t));
ay=zeros(1,length(t));

x(1) = 3; % Initial position in x
y(1) = 0; % Initial position in y
vx(1) = 0; % Initial velocity in x
vy(1) = 5; % Initial velocity in y

[ax(1),ay(1)] = acelrt(x(1),y(1),vx(1),vy(1), m);

for i=1:length(t)-1
    vx(i+1)= vx(i) + ax(i)*dt;
    vy(i+1)= vy(i) + ay(i)*dt;

    x(i+1)= x(i) + vx(i)*dt;
    y(i+1)= y(i) + vy(i)*dt;

    [ax(i+1),ay(i+1)] = acelrt(x(i+1),y(i+1),vx(i+1),vy(i+1), m);
end


function [ax,ay] = acelrt(x,y,vx,vy, m)

K=100;
mu=0.01;

    fx = -mu*vx - K* (x/(x^2+y^2)^1.5);
    fy = -mu*vy - K* (y/(x^2+y^2)^1.5);

    ax= fx/m;
    ay=fy/m;
end


plot(x,y)


v_mod = sqrt(vx.^2 + vy.^2); 
[Int_Indef, Int_def ] = Int(t,v_mod, 0);

%plot(t,Int_Indef)

%INTEGRAL
function [Int_Indef, Int_def ] = Int(x,y, y0)

Int_Indef =zeros(size(y)) ;
areas=diff(x).*(y(1:end-1)+y(2:end))/2 ;
Int_def=sum(areas) ;
Int_Indef=[0, cumsum(areas)];

Int_def=y0 + Int_def ;
Int_Indef=y0 + Int_Indef;

end


L= sqrt(x.^2 + y.^2); 