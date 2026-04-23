%% C3_2023_1_vG.m (duración Control = 4 horas) FUNCIONES USUARIO AL FINAL
clear ; clc ; figure(301) ; 

%% DATOS [SI !!!]
dt=1E-3 ; m=140E-3 ; Ka=1.3E-3 ; KM=0.2E-3 ; g=9.8 ; h=1.5 ; v0=40 ; w0=1800*2*pi/60 ;
N=1E4 ;
[x1,y1,z1,x2,y2,z2,x3,y3,z3]=deal(zeros(N,1)) ; % pre-alojar espacio
[vx1,vy1,vz1,vx2,vy2,vz2,vx3,vy3,vz3]=deal(zeros(N,1)) ; % pre-alojar espacio

z1(1)=h ; z2(1)=h ; z3(1)=h ; % condiciones iniciales
vx1(1)=v0 ; vx2(1)=v0 ; vx3(1)=v0 ;
w1=[0,0,0] ; w2=[0,0,w0] ; w3=[0,w0,0] ;
t=linspace(0,dt*(N-1),N) ; % discretización tiempo

%% 1A PARTICULA CASO 1 [TRAYECTORIAS]
for i=1:N-1 % EULER
    F=ForceG(m)+ForceA(vx1(i),vy1(i),vz1(i),Ka)+ForceM(vx1(i),vy1(i),vz1(i),KM,w1) ;   
    ax=F(1)/m ; ay=F(2)/m ; az=F(3)/m ;
    vx1(i+1)=vx1(i)+ax*dt ; vy1(i+1)=vy1(i)+ay*dt ; vz1(i+1)=vz1(i)+az*dt ;
    x1(i+1)=x1(i)+vx1(i+1)*dt ; y1(i+1)=y1(i)+vy1(i+1)*dt ; z1(i+1)=z1(i)+vz1(i+1)*dt ;  
end
k1=find(z1<0) ; k1fin=k1(1) ; % z<0
subplot(2,1,1) ; plot3(x1(1:k1fin),y1(1:k1fin),z1(1:k1fin)) ; grid on ; title('trayectorias') ;
xlabel('posición X (m)') ; ylabel('posición Y (m)') ; zlabel('posición Z (m)') ; %axis equal ;

%% PARTICULA CASO 2
for i=1:N-1 % EULER
    F=ForceG(m)+ForceA(vx2(i),vy2(i),vz2(i),Ka)+ForceM(vx2(i),vy2(i),vz2(i),KM,w2) ;       
    ax=F(1)/m ; ay=F(2)/m ; az=F(3)/m ;
    vx2(i+1)=vx2(i)+ax*dt ; vy2(i+1)=vy2(i)+ay*dt ; vz2(i+1)=vz2(i)+az*dt ;
    x2(i+1)=x2(i)+vx2(i+1)*dt ; y2(i+1)=y2(i)+vy2(i+1)*dt ; z2(i+1)=z2(i)+vz2(i+1)*dt ;  
end
k2=find(z2<0) ; k2fin=k2(1) ; % z<0
hold on ; plot3(x2(1:k2fin),y2(1:k2fin),z2(1:k2fin)) ; hold off ;

%% PARTICULA CASO 3
for i=1:N-1 % EULER
    F=ForceG(m)+ForceA(vx3(i),vy3(i),vz3(i),Ka)+ForceM(vx3(i),vy3(i),vz3(i),KM,w3) ;           
    ax=F(1)/m ; ay=F(2)/m ; az=F(3)/m ;
    vx3(i+1)=vx3(i)+ax*dt ; vy3(i+1)=vy3(i)+ay*dt ; vz3(i+1)=vz3(i)+az*dt ;
    x3(i+1)=x3(i)+vx3(i+1)*dt ; y3(i+1)=y3(i)+vy3(i+1)*dt ; z3(i+1)=z3(i)+vz3(i+1)*dt ;  
end
k3=find(z3<0) ; k3fin=k3(1) ; % z<0
hold on ; plot3(x3(1:k3fin),y3(1:k3fin),z3(1:k3fin)) ; hold off ;

%% 1B pintar y altura máxima 
z3max=max(z3) ; k3max=find(z3-z3max==0) ;
hold on ; plot3(x3(k3max),y3(k3max),z3(k3max),'ok') ; hold off ;
legend('caso 1','caso 2','caso 3','máxima altura','Location','east') ; %axis equal ;

%% FUERZAS, ENERGIAS y TRABAJOS 
Ec2=0.5*m*(vx2.^2+vy2.^2+vy2.^2)' ; % cinética
[FM2,Fa2,Fg2]=deal(zeros(3,N)) ;
for i=1:N
   FM2(:,i)=ForceM(vx2(i),vy2(i),vz2(i),KM,w2) ;
   Fa2(:,i)=ForceA(vx2(i),vy2(i),vz2(i),Ka) ;
   Fg2(:,i)=ForceG(m) ;
end
[FM2dv,Fa2dv,Fg2dv]=deal(zeros(1,N)) ;
for i=1:N
   FM2dv(i)=dot(FM2(:,i),[vx2(i),vy2(i),vz2(i)]) ;
   Fa2dv(i)=dot(Fa2(:,i),[vx2(i),vy2(i),vz2(i)]) ;
   Fg2dv(i)=dot(Fg2(:,i),[vx2(i),vy2(i),vz2(i)]) ;
end

WM2=mi_integral_cum(t,FM2dv,0) ; Wa2=mi_integral_cum(t,Fa2dv,0) ; Wg2=mi_integral_cum(t,Fg2dv,0) ;
Ee2=Ec2-WM2-Wa2-Wg2 ;

ts2=(1:k2fin) ;
subplot(2,1,2) ;  plot(t(ts2),Ec2(ts2),t(ts2),WM2(ts2),t(ts2),Wa2(ts2),t(ts2),Wg2(ts2),t(ts2),Ee2(ts2)) ;
xlabel('tiempo (s)') ; ylabel('energía o trabajo (J)') ; grid on ;
legend('E cinética','W Magnus','W rozamiento','W gravedad','E_T','Location','west') ; title('caso 2') ;

%% =======================================================================
%% FUNCIONES DE USUARIO
%% =======================================================================
function [F]=ForceG(m) ; F=[0,0,-9.8*m] ; end
function [F]=ForceA(vx,vy,vz,Ka) ; F=-Ka*vecnorm([vx,vy,vz])*[vx,vy,vz] ; end
function [F]=ForceM(vx,vy,vz,KM,w) ; F=KM*cross([vx,vy,vz],w) ; end
function [inty]=mi_integral_cum(x,y,y_ini)
    inty=zeros(size(y)) ;
    inty(2:end)=cumsum(diff(x).*(y(1:end-1)+y(2:end))/2) ;
    inty=inty+y_ini ;
end
