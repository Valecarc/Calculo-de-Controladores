function  [Ts, Zi, Wn]= Fun_AnalisisFun(den_Planta)

%------------------- Analisis de la planta ------------------
syms Z
l_d=length(den_Planta); %Tamaño del vector denominador de la planta
if den_Planta(l_d)==0
    fprintf('La planta es inestable')
    %----- Tipo de sistema ----------
        if l_d==3
            if den_Planta(l_d)==0 
                Tipo=1;
                fprintf(' - Tipo 1')
            end
        elseif l_d==4
            if den_Planta(l_d)==0 && den_Planta(l_d-1)~=0 
                Tipo=1;
                fprintf(' - Tipo 1' )
            end
             if den_Planta(l_d)==0 && den_Planta(l_d-1)==0 
                Tipo=2;
                fprintf(' - Tipo 2' )
            end
        end
        [Ts, Zi, Wn]=[0 0 0];
else
    fprintf('Planta es estable')
    Tipo=0;
    fprintf(' - Tipo 0')
    if l_d==2   
        Tao=(den_Planta(2))^-1;
        Ts=5*Tao;               %-------- Tiempo de establecimiento     
        Wn=0; Zi=0;
    elseif l_d>2       
        Wn=(den_Planta(3))^0.5;
        Zi=vpa(solve(2*Z*Wn==den_Planta(2),'Z'));
        %-------- Tiempo de establecimiento ------------------------
        if  Zi<1    
            %SubAmortiguado
            Ts=vpa(4/(Zi*Wn));           
        elseif Zi==1          
            %Critico 
            Ts=5.83/Wn;           
        elseif Zi>1  && Zi<9         
            %SobreAmortiguado 
%             Ts=(4/Wn)*((1/(Zd+(((Zd^2)+1)^0.5))));   
            Ts=(4/Wn)*((1/(Zi+(((Zi^2)+1)^0.5)))+(1/(((Zi^2)-1)^0.5)));
        end
    end
end
Ts=double(Ts) ;Zi=double(Zi);Wn=double(Wn);
end