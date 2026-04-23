%% C3_2022_3_vG.M   FUNCIONES USUARIO AL FINAL (duración Control = 4 horas)
clear ; clc ; figure(303) ;

%% DATOS
m=80 ; g=9.8 ; r0=[0;1;0] ; v0=[3;1;4] ; 
rb=[0;-5;0] ; lb=24 ; kb1=50 ; kb2=5 ; 
rc=[0;5;0] ; lc=25 ; kc1=60 ; kc2=10 ; 
beta=40 ;

%% INTEGRAR MOVIMIENTO (EULER)
dt=1E-2 ; % incremento de tiempo
N=25/dt ; t=linspace(0,25,N) ; % discretización de tiempo
[rr,vv,aa,ffb,ffc,ffg,ffr,ff]=deal(zeros(3,N)) ; % pre-alojar espacio
rr(:,1)=r0 ; vv(:,1)=v0 ; % condiciones iniciales

for i=1:N-1
   [ffb(:,i),ffc(:,i)]=fuerza(rr(:,i)) ; % fuerza elástica
   ffg(3,i)=-m*g ; % gravedad
   ffr(:,i)=-beta*vv(:,i) ; % rozamiento
   ff(:,i)=ffb(:,i)+ffc(:,i)+ffg(:,i)+ffr(:,i) ; % suma fuerzas
   % EULER
   aa(:,i)=ff(:,i)/m ; % aceleración
   vv(:,i+1)=vv(:,i)+aa(:,i)*dt ; % velocidad
   rr(:,i+1)=rr(:,i)+vv(:,i+1)*dt ; % posición
end
[ffb(:,end),ffc(:,end)]=fuerza(rr(:,end)) ;
ffg(3,end)=-m*g ;
ffr(:,end)=-beta*vv(:,end) ;
ff(:,end)=ffb(:,end)+ffc(:,end)+ffg(:,end)+ffr(:,end) ;

subplot(2,1,1) ; plot3(rr(1,:),rr(2,:),rr(3,:)) ;
title('trayectoria') ;
xlabel('x(m)') ; ylabel('y(m)') ; zlabel('z(m)') ; grid on ; %axis equal ;

%% TRABAJOS, ENERGÍAS
% trabajos elásticos
Webxyz=zeros(3,N) ;
Webxyz(:,2:end)=0+0.5*(ffb(:,2:end)+ffb(:,1:end-1)).*diff(rr,1,2) ;
Web=cumsum(sum(Webxyz)) ;
Wecxyz=zeros(3,N) ;
Wecxyz(:,2:end)=0+0.5*(ffc(:,2:end)+ffc(:,1:end-1)).*diff(rr,1,2) ;
Wec=cumsum(sum(Wecxyz)) ;
We=Web+Wec ;
% trabajo rozamiento
Wrxyz=zeros(3,N) ;
Wrxyz(:,2:end)=0+0.5*(ffr(:,2:end)+ffr(:,1:end-1)).*diff(rr,1,2) ;
Wr=cumsum(sum(Wrxyz)) ;

% energía elástica
[Ebe,Ece]=deal(zeros(1,N)) ;
for i=1:N [Ebe(i),Ece(i)]=energia(rr(:,i)) ; end
Ee=Ebe+Ece ;
% energía cinética
vmod=sqrt(sum(vv.^2,1)) ;
Ek=0.5*m*sum(vv.^2,1) ;
Etotal=Ee+Ek+m*g*rr(3,:) ; % y potencial gravitatoria

subplot(2,1,2) ; plot(t,Etotal,t,We,t,Wr,t,Etotal-Wr) ; grid on ;
legend('E mecánica','W elástico','W rozamiento','E_{tot}-W_r') ;
xlabel('tiempo (s)') ; ylabel('energía/trabajo (J)') ; title('energías y trabajo') ;


%% FUNCIONES DE USUARIO ==================================================
function [fb,fc]=fuerza(r)
    rb=[0;-5;0] ; lb=24 ; kb1=50 ; kb2=5 ; 
    rc=[0;5;0] ; lc=25 ; kc1=60 ; kc2=10 ;
    
    db=sqrt(sum(((r-rb).^2))) ;
    if db<=lb fb=0 ;
    else fb=-((r-rb)/db)*(kb1*(db-lb)+kb2*(db-lb)^2) ; 
    end 
    dc=sqrt(sum(((r-rc).^2))) ;
    if dc<=lc fc=0 ;
    else fc=-((r-rc)/dc)*(kc1*(dc-lc)+kc2*(dc-lc)^2) ; 
    end 
end

function [Eb,Ec]=energia(r)
    rb=[0;-5;0] ; lb=24 ; kb1=50 ; kb2=5 ; 
    rc=[0;5;0] ; lc=25 ; kc1=60 ; kc2=10 ;
    
    db=sqrt(sum(((r-rb).^2))) ;
    if db<=lb Eb=0 ;    
    else Eb=0.5*kb1*(db-lb)^2+(1/3)*kb2*(db-lb)^3 ;
    end 
    dc=sqrt(sum(((r-rc).^2))) ;
    if dc<=lc Ec=0 ; 
    else Ec=0.5*kc1*(dc-lc)^2+(1/3)*kc2*(dc-lc)^3 ;
    end 
end
