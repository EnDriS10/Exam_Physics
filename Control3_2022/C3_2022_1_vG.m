%% C3_2022_1_vG.m (duración Control = 4 horas)
clear ; clc ; figure(301) ;

%% 1A Extracción 
N=1E4 ; % número de extracciones
r=rand(5,N) ; % extracciones aleatorias (distribución de probabilidad constante)
bolas=zeros(5,N) ; % almacenar extracciones (pre-alojar espacio)
for i=1:N
   caja=[16,13,11] ; % número de bolas con números uno, dos y tres. Total son 40 al empezar
   for j=1:5 % ir sacando cinco bolas y actualizar las que quedan en la caja
      % la probabilidad se obtiene de casos favorables entre casos posibles (bolas que quedan)
      if r(j,i)<(caja(1))/(40-j+1) % bola uno 
         bolas(j,i)=1 ; caja(1)=caja(1)-1 ; % escoger bola y sacar de la caja 
      elseif r(j,i)<((caja(1)+caja(2))/(40-j+1)) % bola dos (atención a la cumulativa) 
         bolas(j,i)=2 ; caja(2)=caja(2)-1 ; % escoger bola y sacar de la caja  
      else
         bolas(j,i)=3 ; caja(3)=caja(3)-1 ;  % bola tres, escoger bola y sacar de la caja 
      end      
   end   
end

%% 1B distribución probabilidad, el dominio son números enteros
suma=sum(bolas,1) ;
[v,x]=hist(suma,[5:15]) ;
bar(x,v/N) ; xlabel('valor suma (número entero)') ; ylabel('probabilidad') ;

%% 1C analizar algunos resultados estadísticos
% que la cuarta sea un valor 3
p43=sum(bolas(4,:)==3)/N ; % probabilidad normalizada
% que no saque ningún uno
p1_=sum((bolas==ones(5,N)),1) ; % cuántas de tipo 1
p11_=(p1_==0) ; % ocasiones en las que no aparece un uno
p1=sum(p1_==0)/N ; % probabilidad normalizada

counter123=0 ;
for i=1:N 
   if bolas(1:3,i)==[1;2;3] counter123=counter123+1 ;  end
   if bolas(2:4,i)==[1;2;3] counter123=counter123+1 ;  end
   if bolas(3:5,i)==[1;2;3] counter123=counter123+1 ;  end
   % no hay else...
end
p123=counter123/N ; % probabilidad normalizada

fprintf('probabilidad cuarta bola sea un tres = %f\n',p43) ;
fprintf('probabilidad ninguna sea un uno = %f\n',p1) ;
fprintf('probabilidad secuencia 1-2-3 = %f\n',p123) ;
