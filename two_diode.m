clc
clear
%% datos del catálogo
[File_importpi,PathName_importpi] = uigetfile...
    ('*.dat;*.txt','Seleccione el archivo con los datos');
if isequal(File_importpi,0)
    return
else
    % copiar a la carpeta temp del usuario
    copyfile([PathName_importpi,File_importpi],[tempdir,File_importpi])
    % cargo el pi txt
    filename = fullfile([tempdir,File_importpi]);
    fileID = fopen(filename);
    datosPanelFV = textscan(fileID,'%s %s');
    fclose(fileID);
    fabricante = datosPanelFV{2}{1};
    modelo = datosPanelFV{2}{2};

    Ns = str2double(datosPanelFV{2}{10});
    Np = str2double(datosPanelFV{2}{11});
    largoPV = str2double(datosPanelFV{2}{12});
    anchoPV = str2double(datosPanelFV{2}{13});

    PmaxE = str2double(datosPanelFV{2}{3});
    Vocn = str2double(datosPanelFV{2}{4});
    Vmpp = str2double(datosPanelFV{2}{5});
    Impp = str2double(datosPanelFV{2}{6});
    Iscn = str2double(datosPanelFV{2}{7});
    Ki = str2double(datosPanelFV{2}{8});
    Kv = str2double(datosPanelFV{2}{9});

    Ki = Ki*Iscn/100;
    Kv = Kv*Vocn/100;

end
%% entrada de valores de Temp e Irradiancia
Tt = 25;
T = Tt + 273.15; % temp. de trabajo
Tn = 298.15; % temp. de referencia 25+273=298
dt = Tn - T;
G = 1000; % irradiancia de trabajo del panel
Gn = 1000; % irradiancia de referencia 1000 W/m2
%% constantes físicas
q = 1.60e-19 ; % carga del Electrón
K = 1.38e-23 ; % constante de Boltzman
%% condiciones iniciales
p = 2.2;
n1 = 1; % Se asume
n2 = 1.2; % Se asume
Vt = (Ns*K*T)/q;
Io = (Iscn+Ki*dt)/(exp((Vocn+Kv*dt)/Vt)-1);% Io1 = Io2 % Ecuación 8
Ipv = (Iscn+Ki*dt)*(G/Gn);
TolP = 0.0009;
errorP = 1;
marca = 0;
% global marcaPlot
marcaPlot = 0;
Rs0 = 0.0;
Rs = Rs0;%
Rp0 = (Vmpp/(Iscn-Impp))-((Vocn-Vmpp)/Impp); % Ecuación 11
Rp = Rp0;%
paso = .1;
cont2 = 1;
errorP_ant=100;
%% ecuaciones preliminares
while (errorP >= TolP)
    Rs_final = Rs;
    Rp_final = Rp;
    cont = 1;
    for V = 0: paso: Vocn
        x = Impp; % valor inicial de I
        Tol = 0.0009; % tolerancia o error
        count = 0; % contador de iteraciones
        error = 1; %derivada de la ftn
        %ecuaciones preliminares
        f1 = @(x) Ipv-((Iscn+Ki*dt)/(exp((Vocn+Kv*dt)/Vt)-1))...
            *(exp((V+x*Rs)/(n1*Vt))-1)-((Iscn+Ki*dt)...
            /(exp((Vocn+Kv*dt)/Vt)-1))*(exp((V+x*Rs)/...
            (n2*Vt))-1)-((V+x*Rs)/Rp)-x; % Ecuación 1
        f = f1(x);% evaluo f para el valor inicial de I
        fini_deriv = diff(sym(f1),1);

        fprintf('iter       V           I          f(I)        errorI       P         errorP        Rs          Rp \n')
        while  ((error > Tol || abs(f) > Tol) && x > 0)%
            count = count + 1;
            fprime = eval(fini_deriv);
            xnew = x - (f/fprime);%
            error = abs(x-xnew); 
            x = xnew;
            f = f1(x);% evaluo el nuevo valor de f(I)
            P = V * x;  % calulo la potencia máxima PmaxC 
            fprintf('%3i%12.4f%12.4f%12.4f%12.4f%12.4f%12.4f%12.3f%12.4f \n',...
                count,V,x,f,error,P,errorP,Rs,Rp)
        end
        I_array(cont) = x;
        V_array(cont) = V;
        P_array(cont) = P;
        cont = cont+1;
    end
    PmaxC = max(P_array); % encontrar el valor máximo de PmaxC
    errorP = abs(PmaxC - PmaxE); % Error = |PmaxC-PmaxE|
    %%
    figure(22)
    plot(V_array,P_array,'b');
    axis([0 Vocn+5 0 PmaxE+10]);
    xlabel('Voltage (V)');
    ylabel('Power (W)');
    title('P-V Characteristic');
    hold all
    %%
    figure(23)
    plot(V_array,I_array,'b');
    axis([0 Vocn+5 0 max(I_array)+0.5]);
    xlabel('Voltage (V)');
    ylabel('Current (A)');
    title('I-V Characteristic');
    hold all
    %%
    if marca == 0
        Rs = Rs + 0.001;
    elseif (errorP<0.1)||(errorP_ant~=100)||(errorP_ant-errorP)<0%si el error aumenta se detiene el programa
        Rs = Rs + 0.001;
        msgbox('Se alcanzó el menor error posible','Fin de programa')
        break
        errorP_ant = errorP;% guardo el valor anterior del error
    elseif errorP > 4
        Rs = Rs + 0.1;
    elseif errorP >= 0.1
        Rs = Rs + 0.001;
    end
    marca = 1;
    Rp = (Vmpp+Impp*Rs)/(Iscn-Io*(exp((Vmpp+Impp*Rs)/Vt)+exp((Vmpp+Impp*Rs)/((p-1)*Vt))+2)-(PmaxE/Vmpp)); % Ecuación 10
    if Rp <= 0
        errordlg('Programa detenido Rp <= 0','ERROR')
        break
    end
    PmaxC_array(cont2) = PmaxC;
    Rs_array(cont2) = Rs;
    Rp_array(cont2) = Rp;
    errorP_array(cont2) = errorP;
    cont2 = cont2 + 1;
end
fprintf(' \n')
fprintf('\nPmaxC = %g, errorP = %g, Rs = %g, Rp = %g, Io = %g, Ipv = %g \n',...
    max(P_array),errorP,Rs_final,Rp_final,Io,Ipv)
%% salvo variables de salida V I P
save (['I_',num2str(Tt),'_',num2str(G),'_',modelo,'.mat'], 'I_array')
save (['V_',num2str(Tt),'_',num2str(G),'_',modelo,'.mat'], 'V_array')
save (['P_',num2str(Tt),'_',num2str(G),'_',modelo,'.mat'], 'P_array')
save (['PmaxC_',num2str(Tt),'_',num2str(G),'_',modelo,'.mat'], 'PmaxC_array')
save (['Rs_',num2str(Tt),'_',num2str(G),'_',modelo,'.mat'], 'Rs_array')
save (['Rp_',num2str(Tt),'_',num2str(G),'_',modelo,'.mat'], 'Rp_array')
save (['errorP_',num2str(Tt),'_',num2str(G),'_',modelo,'.mat'], 'errorP_array')
%
load (['I_',num2str(Tt),'_',num2str(G),'_',modelo,'.mat'], 'I_array')
load (['V_',num2str(Tt),'_',num2str(G),'_',modelo,'.mat'], 'V_array')
load (['P_',num2str(Tt),'_',num2str(G),'_',modelo,'.mat'], 'P_array')
load (['PmaxC_',num2str(Tt),'_',num2str(G),'_',modelo,'.mat'], 'PmaxC_array')
load (['Rs_',num2str(Tt),'_',num2str(G),'_',modelo,'.mat'], 'Rs_array')
load (['Rp_',num2str(Tt),'_',num2str(G),'_',modelo,'.mat'], 'Rp_array')
load (['errorP_',num2str(Tt),'_',num2str(G),'_',modelo,'.mat'], 'errorP_array')

%% grafico
figure(1)
plot(V_array,I_array,'LineWidth',1);
grid on
axis([0 Vocn+2 0 Iscn+1]);
xlabel('Voltage (V)');
ylabel('Current (A)');
% title('I-V Characteristic');
title({'I-V Characteristic';...
    ['Temperature: ',num2str(Tt),' °C'];['Irradiance: ',...
    num2str(1000),' W/m^2']});
% text(4,max(I_array)+.3,{[num2str(G),' W/m^2']},'HorizontalAlignment','left')
hold all
%
figure(2)
plot(V_array,P_array,'LineWidth',1);
axis([0 Vocn+2 0 PmaxE+10]);
xlabel('Voltage (V)');
ylabel('Power (W)');
% title('P-V Characteristic');
title({'P-V Characteristic';...
    ['Temperature: ',num2str(Tt),' °C'];['Irradiance:',...
    num2str(1000),' W/m^2']});
% legend('-DynamicLegend','location','Best')
legend('1000 W/m^2','800 W/m^2','600 W/m^2','Location','Best');
% text(Vmpp-5, max(P_array)+3,{[num2str(G),' W/m^2']},'HorizontalAlignment','left')
hold all
%
figure(3)
plot3(V_array,I_array,P_array,'LineWidth',1);
axis([0 Vocn+.5 0 Iscn+.5 0 PmaxE+5]);
xlabel('Voltage (V)');
ylabel('Current (A)');
zlabel('Power (W)');
% title('Curves P-I-V');
title({'P-I-V Characteristic';...
    ['Temperature: ',num2str(Tt),' °C'];['Irradiance:',...
    num2str(1000),' W/m^2']});
legend('-DynamicLegend','location','Best')
grid on
hold on

% % grafico de superficie
figure(4)
[T,G] = meshgrid(0:25:75,200:200:1000);
Ipv = (Iscn+Ki*(T-Tn)).*(G./Gn);
surf(G,T,Ipv);
ylabel('Temperature (°C)');
xlabel('Irradiance (W/m^2)');
zlabel('Voltage (V)');
shading interp
% colorbar
% grafico de P vs Rs
figure(5)
plot(Rs_array,PmaxC_array,'LineWidth',1);
axis([0 max(Rs_array)+.005 0 max(PmaxC_array)]);
xlabel('Rs (\Omega)');
ylabel('Power (W)');
% title('Curve PmaxC vs Rs');
title({'P-V Characteristic';...
    ['Temperature: ',num2str(Tt),' °C'];['Irradiance:',...
    num2str(1000),' W/m^2']});
% legend('-DynamicLegend','location','Best')
legend('1000 W/m^2','800 W/m^2','600 W/m^2','Location','Best');
% text(Vmpp-5, max(P_array)+3,{[num2str(G),' W/m^2']},'HorizontalAlignment','left')
hold all
% % grafico de P vs Rp
figure(6)
plot(Rp_array,PmaxC_array,'LineWidth',1);
xlabel('Rp (\Omega)');
ylabel('Power (W)');
title('Curve PmaxC vs Rp');
title({'P-V Characteristic';...
    ['Tc: ',num2str(Tt),' °C'];['Irradiance:',...
    num2str(1000),' W/m^2']});
legend('-DynamicLegend','location','Best')
legend('1000 W/m^2','800 W/m^2','600 W/m^2','Location','Best');
% text(Vmpp-5, max(P_array)+3,{[num2str(G),' W/m^2']},'HorizontalAlignment','left')
hold all
% % grafico de error vs Rs
figure(7)
plot(Rs_array,errorP_array,'LineWidth',1);
xlabel('Rs (\Omega)');
ylabel('Error');
% title('Curve Error vs Rs');
title({'P-V Characteristic';...
    ['Tc: ',num2str(Tt),' °C'];['Irradiance:',...
    num2str(1000),' W/m^2']});
legend('-DynamicLegend','location','Best')
legend('1000 W/m^2','800 W/m^2','600 W/m^2','Location','Best');
% text(Vmpp-5, max(P_array)+3,{[num2str(G),' W/m^2']},'HorizontalAlignment','left')
hold all