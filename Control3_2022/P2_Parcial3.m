clear;clc;

%%
% ================================
% Datos
% ================================

R1 =2; R2=1; R3=2; R4 =2; R5=1; R6=1;
VA= 3; VB=2; VC=3;

R= [R1,R2,R3,R4,R5,R6]; V= [VA,VB,VC];

%%
% ================================
% 1 A
% ================================
[I , Pv, Pr ] =Circuit_C3(R, V);

%%
% ================================
% 1B
% ================================
dR=1E4;
R6 = (0 : 1/dR :3);

    I_R6 = zeros(6,length(R6));
    Pv_R6 = zeros(3,length(R6));
    Pr_R6 = zeros(6,length(R6));
for i= 1:length(R6)
    R(6) = R6(i);
    [I, Pv, Pr] = Circuit_C3(R, V);
    I_R6(:,i) = I; Pv_R6(:,i) = Pv; Pr_R6(:,i) = Pr;
end

figure (1)
 subplot(2,3,1) ;
 plot(R6,I_R6(1,:),'r') ;
 xlabel('Resistencia 3 (Ohmnios)') ;
 ylabel('Intensidad 1 (A)') ;
 title('Intensidad que pasa por la Resistencia 1') ; 
 grid on ; % mostrar retícula

  subplot(2,3,2) ;
 plot(R6,I_R6(2,:),'r') ;
 xlabel('Resistencia 3 (Ohmnios)') ;
 ylabel('Intensidad 2 (A)') ;
 title('Intensidad que pasa por la Resistencia 2') ;    
 grid on ; % mostrar retícula

  subplot(2,3,3) ;
 plot(R6,I_R6(3,:),'r') ;
 xlabel('Resistencia 3 (Ohmnios)') ;
 ylabel('Intensidad 3 (A)') ;
 title('Intensidad que pasa por la Resistencia 3') ;    
 grid on ; % mostrar retícula

  subplot(2,3,4) ;
 plot(R6,I_R6(4,:),'r') ;
 xlabel('Resistencia 3 (Ohmnios)') ;
 ylabel('Intensidad 4 (A)') ;
 title('Intensidad que pasa por la Resistencia 4') ; 
 grid on ; % mostrar retícula

   subplot(2,3,5) ;
 plot(R6,I_R6(5,:),'r') ;
 xlabel('Resistencia 3 (Ohmnios)') ;
 ylabel('Intensidad 5 (A)') ;
 title('Intensidad que pasa por la Resistencia 5') ;   
 grid on ; % mostrar retícula

   subplot(2,3,6) ;
 plot(R6,I_R6(6,:),'r') ;
 xlabel('Resistencia 3 (Ohmnios)') ;
 ylabel('Intensidad 6 (A)') ;
 title('Intensidad que pasa por la Resistencia 6') ;   
 grid on ; % mostrar retícula

 figure (2)
 subplot(2,3,1) ;
 plot(R6,Pr_R6(1,:),'r') ;
 xlabel('Resistencia 3 (Ohmnios)') ;
 ylabel('Potencia 1 (W)') ;
 title('Potencia disipada por la Resistencia 1') ; 
 grid on ; % mostrar retícula

  subplot(2,3,2) ;
 plot(R6,Pr_R6(2,:),'r') ;
 xlabel('Resistencia 3 (Ohmnios)') ;
 ylabel('Potencia 2 (W)') ;
 title('Potencia disipada por la Resistencia 2') ;    
 grid on ; % mostrar retícula

  subplot(2,3,3) ;
 plot(R6,Pr_R6(3,:),'r') ;
 xlabel('Resistencia 3 (Ohmnios)') ;
 ylabel('Potencia 3 (W)') ;
 title('Potencia disipada por la Resistencia 3') ;    
 grid on ; % mostrar retícula

  subplot(2,3,4) ;
 plot(R6,Pr_R6(4,:),'r') ;
 xlabel('Resistencia 3 (Ohmnios)') ;
 ylabel('Potencia 4 (W)') ;
 title('Potencia disipada por la Resistencia 4') ; 
 grid on ; % mostrar retícula

   subplot(2,3,5) ;
 plot(R6,Pr_R6(5,:),'r') ;
 xlabel('Resistencia 3 (Ohmnios)') ;
 ylabel('Potencia 5 (W)') ;
 title('Potencia disipada por la Resistencia 5') ;   
 grid on ; % mostrar retícula

   subplot(2,3,6) ;
 plot(R6,Pr_R6(6,:),'r') ;
 xlabel('Resistencia 3 (Ohmnios)') ;
 ylabel('Potencia 6 (W)') ;
 title('Potencia disipada por la Resistencia 6') ;   
 grid on ; % mostrar retícula

 %%
% ================================
% 1C
% ================================

for i= 1:length(R6)
 Pot(i)=sum(Pr_R6(:,i));
end
P_VC= Pv_R6(3,:);

y=  (Pot/2) - P_VC;
x= R6;
[zeros]= mi_zeros_autom(x,y)

%%
% ================================
% Funciones
% ================================
function [I , Pv, Pr ] =Circuit_C3(R, V)
    EC1= [1 -1 0 1 -1 0];
    EC2= [0 1 -1 0 1 -1];
    EC3= [-R(1) 0 0 R(4) 0 0];
    EC4= [0 -R(2) 0 0 R(5) 0];
    EC5= [0 0 -R(3) 0 0 R(6)];
    EC6= [0 0 0 R(4) R(5) R(6)];

    A= [EC1;EC2;EC3;EC4;EC5;EC6];
    b1= [0 0 V(1) V(2) V(3) ( V(1) + V(2) + V(3) )];
    b=b1';
    I1= inv(A)*b;
    I=I1';
    Pv = V.*I(4:end);
    Pr = R.*I.^2;
end

 function [zeros]=mi_zeros_autom(x,y)

ff=y(1:end-1).*y(2:end) ; % atención, tiene un elemento menos
 x_mid=(x(1:end-1)+x(2:end))/2 ; % puntos medios de los intervalos
 b_negativos=ff<0 ;      % generar vector lógico
 zeros=x_mid(b_negativos); % extraer b_negativos de x_mid

end