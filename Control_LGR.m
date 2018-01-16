clc ,clear all
syms z

%___________FUNCION DE TRANSFERENCIA___________

num_Planta=[364.4]; 
den_Planta=[1 56.1 522.8];

% num_Planta=[6.202]; 
% den_Planta=[1 8.856];
%___________ Analisis de la planta ___________
[Ts, Zi, Wn]=Fun_AnalisisFun(den_Planta)

%____________ Condiciones ________________ 
Entrada=1;
Zd=0.99 ;     Tsd=Ts*0.9;     Wnd=4/(Zd*Tsd);

%______DISCRETIZACION F TRANSFERENCIA___________
T=Tsd/10
F_C=tf(num_Planta,den_Planta) %__Funcion Continua
F_D=c2d(F_C,T,'zoh')          %__Funcion Discreta
% 
%_______________POLOS Y CEROS FT__________________
Polos=roots(F_D.den{1})  
Ceros=roots(F_D.num{1}) 
n=length(Polos); m=length(Ceros);

%_______________POLOS EN 1 NECESARIOS__________________
Tipo=0; j_i=1; 
while j_i<=n
    if Polos(j_i)==1
        Tipo=Tipo+1   
    end
    j_i=1+j_i;
end
Num_P1=Entrada-Tipo  %Polos en 1 necesarios

%_________ NUEVA FUNCION DE TRANSFERENCIA___________
Polos_1=tf([1],[1 -1],T);
F_D_Nueva=F_D*(Polos_1^Num_P1)
Polos_N=round(roots(F_D_Nueva.den{1}),3)
Ceros_N=roots(F_D_Nueva.num{1}) 
n_N=length(Polos_N);
m_N=length(Ceros_N);

% pzmap(F_D)    % Grafica polos y ceros
% hold on
%_______________P DESEADO_______________
Wd=Wnd*(1-Zd^2)^0.5;
Zd_Mag=exp(-T*Zd*Wnd);  
Zd_ang=T*Wd;
[Zd_RE,Zd_IM]=pol2cart(Zd_ang,Zd_Mag)
% plot(Zd_RE,Zd_IM,'X')
% hold on
ZP_deseado=Zd_RE+Zd_IM*i

%____________Angulo Compensador_______________ 180-AgPolos+AgCeros

j_p=1 ;j_c=1;

if n_N>0         %____Angulos Polos_____ tan=(y2-y1/x2-x1)
    Ang_Polos=0;
    while j_p<=n_N
        if Polos_N(j_p)==1
            Ang_Polos(j_p)=atan2d(Zd_IM-imag(Polos_N(j_p)),Zd_RE-real(Polos_N(j_p)));
        else
            Alfa(j_p)=Polos_N(j_p);
        end
            j_p=1+j_p;
    end  
    Ang_Comp=180-sum(Ang_Polos);
end
Ang_Polos

if m_N>0        %_____Angulos Ceros_____
    Ang_Ceros=0;
    while j_c<=m_N
          Ang_Ceros(j_c)=atan2d(Zd_IM-imag(Ceros(j_c)),Zd_RE-real(Ceros(j_c)))  
          j_c=1+j_c;
    end
    Ang_Comp=Ang_Comp+sum(Ang_Ceros);
end
Ang_Comp
 
%__________________Constantes del compensador_______________

        %______ Se eliminan polos  (ALFA) _________
Alfa         % vector de polos a eliminar
j_alfa=1;
Num_Comp=1;
 while j_alfa<=length(Alfa)
          if Alfa(j_alfa)==0
              j_alfa=j_alfa+1;
          else   
              Num_Comp= Num_Comp*tf([1 -Alfa(j_alfa)],[1],T);
              j_alfa=j_alfa+1;
          end
 end
Num_Comp
        %______ Se Compensa angulo (Beta) _________
Beta=Zd_RE-(Zd_IM/tand(Ang_Comp))
Den_Comp=tf([1],[1 -Beta],T)

Control=Num_Comp*Den_Comp

        %______ Ganancia (K) _________
syms  K_s
F_Comp=vpa(poly2sym(Control.num{1},z)/poly2sym(Control.den{1},z),4)
F_plantaD=vpa(poly2sym(F_D_Nueva.num{1},z)/poly2sym(F_D_Nueva.den{1},z),4)

K=vpa(subs(K_s*F_Comp*F_plantaD,z,ZP_deseado),4);
K=double(solve(abs(K)==1,K_s))

%___________________ FT CONTROL ________________________________
Control=K*Num_Comp*Den_Comp*(Polos_1^Num_P1)

roots(Control.num{1})
roots(Control.den{1})
%_______________POLOS Y CEROS ________________________________
Planta_Final=Control*F_D/(1+Control*F_D);
pzmap(Planta_Final)

