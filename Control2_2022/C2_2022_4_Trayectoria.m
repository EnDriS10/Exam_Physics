%% C2_2022_4_Trayectoria.m
clear ; clc ; figure(204) ;

%% DATOS
Dt=1E-4 ; m=3 ; w=0.7 ; % UNIDADES SI

t=0:Dt:20 ; N=length(t) ; % tiempo y número de pasos
[x,y,vx,vy]=deal(zeros(N,1)) ; % pre-alojar espacio
x(1)=1 ; y(1)=0 ; vx(1)=-1 ; vy(1)=-1 ; % condiciones iniciales

%% INTEGRAR POR EL MÉTODO DE EULER
for i=1:N-1
   dt=t(i+1)-t(i) ; % incremento de tiempo
   fx=-x(i)-4*sin(w*t(i))^2 ; % fuerza X
   fy=-y(i)-4*cos(w*t(i))^2 ; % fuerza Y
   ax_now=fx/m ; ay_now=fy/m ; % aceleraciones
   vx(i+1)=vx(i)+ax_now*dt ; vy(i+1)=vy(i)+ay_now*dt ; % componentes velocidad
   x(i+1)=x(i)+vx(i)*dt ; y(i+1)=y(i)+vy(i)*dt ; % posición
end

%% CONTINUACIÓN: INTEGRAR CON MÉTODO DE HEUN... 

%% REPRESENTACIÇON GRÁFICA DE LA TRAYECTORIA
plot(x,y,x(1),y(1),'or') ; axis equal ; title('trayectoria') ;
xlabel('posición x (m)') ; ylabel('posición y (m)') ; grid on ;
legend('trayectoria','posición inicial') ;

