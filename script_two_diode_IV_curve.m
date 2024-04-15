%% Metadata
% Author: Liomnis Osorio Laurencio
% Date: 2024-02-16
% Description: This script for...
% ----------------------------------------------------------------------
%% Cleaning workspace
% ----------------------------------------------------------------------
clc
clear
% ----------------------------------------------------------------------
%% leer datos del módulo FV
% ----------------------------------------------------------------------
datapi = readmatrix('KC200GT.txt');
% ----------------------------------------------------------------------
%% Parámetros del datasheet en condiciones STC
% ----------------------------------------------------------------------
% estos parametros están disponibles en las fichas técnicas de todos los
% módulos fotovoltaicos
% parametros climaticos
G_tr = 1000;            % (W/m2)
T_tr = 25;              % (°C)
% parametros electricos
Pmpp = datapi(1,2);   % (W)
Voc = datapi(2,2);    % (V)
Vmpp = datapi(3,2);   % (V)
Impp = datapi(4,2);   % (A)
Isc = datapi(5,2);    % (A)
Ki = datapi(6,2); 
Kv = datapi(7,2); 
eff = datapi(13,2);   % (%)
% parametros mecanicos
Ns = datapi(8,2);     % (Unidades)
largo = datapi(10,2); % (m)
ancho = datapi(11,2); % (m)
NOCT = datapi(12,2);  % (°C)
% ----------------------------------------------------------------------
%% valores de Temperatura e Irradiancia para obtener la curva I-V
% ----------------------------------------------------------------------
Tt = T_tr;
T = Tt + 273.15; % temp. de trabajo
Tn = 298.15; % temp. de referencia, por defecto se deja igual que
% en el datasheet 25+273=298
dT = T-Tn;
G = G_tr; % irradiancia de trabajo del panel, por defecto se deja igual 
% que en el datasheet 1000 W/m2
Gn = 1000; % irradiancia de referencia 1000 W/m2
% ----------------------------------------------------------------------
%% constantes físicas
% ----------------------------------------------------------------------
q = 1.60e-19 ; % carga del Electrón
K = 1.38e-23 ; % constante de Boltzman
% ----------------------------------------------------------------------
%% Condiciones iniciales
% ----------------------------------------------------------------------
Vt = Ns*((K*T)/q);
% ----------------------------------------------------------------------
%% Parámetros internos que se deben calcular
% estos son los parametros internos desconocidos que se deben calcular, 
% en este caso, ya fueron calculados aplicando otro método, por lo que 
% cual se puede graficar la curva I-V
% ---------- definidos en la página 2 del documento --------------------
Iph = (Isc+Ki*(T-Tn))*(G/Gn);             % (A)
Io1 = datapi(19,2);
Io2 = datapi(20,2);                       % (A)
Rsh = datapi(16,2);                       % (ohm)
Rs = datapi(15,2);                        % (ohm)
n1 = datapi(17,2);                        % (adimensional)
n2 = datapi(18,2);                        % (adimensional)
% ----------------------------------------------------------------------
%% Método Newton-Rapshon para resolver double diode model (DDM)
% ----------------------------------------------------------------------
% en este método se dan valores a V de forma tal que varie cada Δ = 0.1 V
% para ir calculando el valor de I que corresponde a cada punto de V
% de esta forma se hace variar a V en el intervalo de 0 < V < Voc, para
% calcular el valor de I en cada punto de la curva.
cont = 1;
deltaV = .1; % Δ = 0.1 V
pointscurves = 0: deltaV: Voc; % puntos de la curva I-V
% definir arreglos de I donde se guardan los puntos de la curva
I_array = zeros(1, length(pointscurves));
% definir arreglos de V donde se guardan los puntos de la curva
V_array = zeros(1, length(pointscurves));
% definir arreglos de P donde se guardan los puntos de la curva
P_array = zeros(1, length(pointscurves));
for V = pointscurves
    I = Impp;       % valor inicial de I
    Tol = 0.001;    % tolerancia o error
    count = 0;      % contador de iteraciones
    error = 1;      % derivada de la ftn
    % ------------- Inicio - Ecuacion 5 --------------------------------
     f1 = @(I) ...
        Iph - ...
        Io1*(exp((V+I*Rs)/(n1*Vt))-1) - ... % Ecuación 2 --> ID1
        Io2*(exp((V+I*Rs)/(n2*Vt))-1) - ... % Ecuación 3 --> ID2
        ((V+I*Rs)/Rsh)-I;                   % Ecuación 4 --> Ish
    % ------------- Fin - Ecuacion 5 -----------------------------------
    f =  f1(I);
    fini_deriv = diff(sym(f1),1);
    % --------- Método iterativo Newton-Rapshon ------------------------
    fprintf('iter       V           I            f(I)         error \n')
    while  (error > Tol || abs(f) > Tol)
        count = count + 1;
        fprime = eval(fini_deriv);
        xnew = I - (f/fprime);%
        error = abs(I-xnew);%
        I = xnew;
        f = f1(I);  % evaluo el nuevo valor de f(I)
        fprintf('%3i %12.4f %12.4f %12.4f %12.4f \n',count,V,I,f,error)
    end
    % --------- guardo los resultados en arreglos ----------------------
    I_array(1,cont) = I;
    V_array(1,cont) = V;
    P = V * I;
    P_array(1,cont) = P;
    cont = cont + 1;
    if I < 0 % si hay valores negativos de I se ignorarán
        break
    end
end
%% CALCULAR ERROR EN Pmpp
% ------------- determinar el punto de máxima potencia Pmpp ------------
p = P_array;
p = find(p==max(P_array));
Vmpp_calc = V_array(p);
Impp_calc = I_array(p);
Pmpp_calc = Vmpp_calc*Impp_calc;
Isc_calc = max(I_array);
Voc_calc = max(V_array);
ErrorMAE = abs(Pmpp - Pmpp_calc);

% % Crear tablas con los arreglos
results_calc = [Pmpp_calc,Isc_calc,Impp_calc,Voc_calc,Vmpp_calc,ErrorMAE];
results_Exp = [Pmpp, Isc, Impp, Voc, Vmpp, zeros(size(Pmpp))];
results = [results_calc;results_Exp];

% -------------------------------------------------------------------------
% Saving in txt file
% Creating data table
ResultsX = {'Calculated'; 'Experimental'};
varSummary = results;
% Opening file
fileID = fopen('_results_parameters_.txt','w');
% Write the headings
fprintf(fileID, '%s \t %s \t %s \t %s \t %s \t %s \t %s \n', ...
    '  ', 'Pmpp', 'Isc', 'Impp', 'Voc','Vmpp','ErrorP');
% Write the data
for i = 1:numel(ResultsX)
    fprintf(fileID, '%s\t', ResultsX{i});
    fprintf(fileID, '%.4f\t', varSummary(i, :));
    fprintf(fileID, '\n');
end
% Closing file
fclose(fileID);
% ----------------------------------------------------------------------
%% Guardar en txt los puntos de I-V curve
% ----------------------------------------------------------------------
% Obtainig date and time
currentDateTime = datetime('now', 'Format', 'yyyyMMdd_HHmmss');
% Name of txt file
GT = char({[num2str(round(G)),'_',num2str(round(Tt))]});
fileName = ['_points_of_Curve_IV_',GT,char(currentDateTime), '_.txt'];
% Abrir el archivo para escritura
fid = fopen(fileName, 'w');
% Verificar si el archivo se abrió correctamente
if fid == -1
    error('No se pudo abrir el archivo para escritura');
end
% Escribir encabezados
fprintf(fid, 'P (W)\tV (V)\tI (A)\n');
% Escribir los vectores en el archivo
for i = 1:length(P_array)
    fprintf(fid, '%.2f\t%.2f\t%.2f\n', ...
        P_array(i), V_array(i), I_array(i));
end
% Cerrar el archivo
fclose(fid);
% ----------------------------------------------------------------------
%% Plotting I-V curve
% ----------------------------------------------------------------------
figure
plot(V_array, I_array, '-.',...
    'Color',[0.494117647409439 0.184313729405403 0.556862771511078], ...
    'LineWidth',2);
hold on
plot(Voc_calc, 0, 'ko','MarkerSize',8,'MarkerFaceColor',[1 .6 .6])
hold on
plot(0, Isc_calc, 'ko','MarkerSize',8,'MarkerFaceColor',[1 .6 .6])
hold on
plot(Vmpp_calc, Impp_calc, 'ko','MarkerSize',8,'MarkerFaceColor',[1 .6 .6])
hold on
text(2,Isc_calc-mean(I_array),{['G = ',num2str(G_tr),' W/m^2','  ',...
    'T = ',num2str(T_tr),' °C']},...
    'HorizontalAlignment','left','FontSize',11);
hold on
text(Voc+0.2, 0+0.5,{'(Voc, 0)'},...
    'HorizontalAlignment','left','FontSize',11);
hold on
text(0+0.2, Isc_calc+0.5,{'(0, Isc)'},...
    'HorizontalAlignment','left','FontSize',11);
hold on
text(Vmpp+0.2, Impp+0.5,{'(Vmpp, Impp)'},...
    'HorizontalAlignment','left','FontSize',11);
axis([0 Voc+5 0 max(I_array)+1]);
xlabel('Voltage (V)');
ylabel('Current (A)');
title('I-V characteristic curve');
set(gca,'xtick', 0:2:max(V_array));% valores del eje X
set(gca,'fontname','Arial','FontSize', 10);  % Set fontname
hold all
% ----------------------------------------------------------------------