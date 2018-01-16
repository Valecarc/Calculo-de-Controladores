clc ,clear all
syms z

%___________FUNCION DE TRANSFERENCIA___________

%%%%  VANESSA ES BOBA %%%%
num_Planta=[364.4]; 
den_Planta=[1 56.1 522.8];

% num_Planta=[6.202]; 
% den_Planta=[1 8.856];
%___________ Analisis de la planta ___________
[Ts, Zi, Wn]=Fun_AnalisisFun(den_Planta)

%____________ Condiciones ________________ 
Entrada=1; % Escalon
Zd=0.9 ;     Tsd=Ts*0.9;     Wnd=4/(Zd*Tsd);

%______DISCRETIZACION F TRANSFERENCIA___________
T=Tsd/10
F_C=tf(num_Planta,den_Planta) %__Funcion Continua
F_D=c2d(F_C,T,'zoh')          %__Funcion Discreta

% --------------PID -------------------
syms s z q0 q1 q2 q3 So S1 S2

N_PlantaD=vpa(poly2sym(F_D.num{1},z),4);
D_plantaD=vpa(poly2sym(F_D.den{1},z),4);

%------------ TIPO DE CONTROL -------------
%%% ... Falta mejorar seccion ...
% Tipo_Control=1
% N_pidD=(q0*z+q1);                   %____PD
% D_pidD= z;

% N_pidD=(z)*(q0+q1/z);        %____PD con filtro
% D_pidD=(z)*(1-(So/z));

% Tipo_Control=2
% N_pidD=(q0*z+q1);                   %____PI
% D_pidD=z-1;                         % Orden 1 Tipo 0

% Tipo_Control=3
% N_pidD=(q0*(z^2)+q1*z+q2);          %____PID
% D_pidD=(z^2)-z;

Tipo_Control=4;
N_pidD=(q0*(z^2)+q1*z+q2);          %___PID con filtro 
D_pidD=(z-1)*(z+So);                % Orden 2 Tipo 0

% Tipo_Control=5
% N_pidD=(q0*(z^3)+q1*(z^2)+q2*z+q3);   %___PID 2 filtro  
% D_pidD=(z-1)*(z-So)*(z-S1);           % Orden 3 Tipo 0


%------------ POLINOMIO DEL CONTROL -------------
pretty((vpa(N_pidD/D_pidD,4))+(vpa(N_PlantaD/D_plantaD)));
Fun_pidD=vpa(collect((D_pidD*D_plantaD)+(N_pidD*N_PlantaD)),4);
Fun_pidD=vpa(fliplr(coeffs(Fun_pidD,z)),2);
G_pidD=(length(Fun_pidD)-1) % GRADO DEL PID

%_____________ Pol Deseado DISCRETO _________________________
Wd=Wnd*(1-Zd^2)^0.5;
Z_Mag=exp(-T*Zd*Wnd);  Z_ang=T*Wd;
[RE,IM]=pol2cart(Z_ang,Z_Mag);

PD_disc=vpa(collect((z-RE+IM*i)*(z-RE-IM*i)*(z-0.05)^(G_pidD-2)),4);

PD_disc=vpa(fliplr(coeffs(PD_disc)),4)  % Coeficientes POL deseado

%__________________________ CONSTANTES ________________________________

if Tipo_Control==2       %%---------- PI --------------
    [q0 q1 ]=solve(PD_disc==Fun_pidD,[q0 q1])
    Q0=double(q0); Q1=double(q1); 
    Q_Discreto= [Q0 Q1]    
    Q_den= [1 -1]
    
elseif Tipo_Control==3  %%---------- PID ------------
    [q0 q1 q2 ]=solve(PD_disc==Fun_pidD,[q0 q1 q2])
    Q0=double(q0); Q1=double(q1); Q2=double(q2);
    Q_Discreto= [Q0 Q1 Q2]  
    Q_den= [1 -1]
    
elseif Tipo_Control==4   %%------- PID con filtro -----
    [So q0 q1 q2 ]=solve(PD_disc==Fun_pidD,[So q0 q1 q2])
    Q0=double(q0); Q1=double(q1); Q2=double(q2); S0=double(So)
    Q_Discreto= [Q0 Q1 Q2]    
    Q_den= [1 -S0-1 S0]    
elseif Tipo_Control==5   %%------- PID con filtro -----
    fliplr(coeffs(collect((z-1)*(z-So)*(z-S1),z),z))
    [So S1 q0 q1 q2 q3]=solve(PD_disc==Fun_pidD,[So S1 q0 q1 q2 q3])
    Q0=double(q0); Q1=double(q1); Q2=double(q2);  Q3=double(q3);
    S0=double(So); S1=double(S1) ;
    Q0=Q0(1); Q1=Q1(1);Q2=Q2(1);Q3=Q3(1);
    S0=S0(1); S1=S1(1);
    Q_Discreto= [Q0 Q1 Q2 Q3]    
    Q_den=[ 1 -S1-S0-1 S1+S0+S1*S0 -S1*S0]   
end

%___________________ Funcion Trans CONTROL ____________________________
F_PID_D=tf(Q_Discreto,Q_den,T)
Control=F_PID_D
%_______________________ POLOS Y CEROS ________________________________
Planta_Final=(F_PID_D *F_D)/(1+F_PID_D*F_D);
pzmap(Planta_Final)



























% 
% 
% %___________________ FT CONTROL ________________________________
% Control=K*Num_Comp*Den_Comp*(Polos_1^Num_P1)
% 
% roots(Control.num{1})
% roots(Control.den{1})
% %_______________POLOS Y CEROS ________________________________
% Planta_Final=Control*F_D/(1+Control*F_D);
% pzmap(Planta_Final)

