%% traccia Hess Smith (2024)

clc
close all
clear 

addpath mat_functions

%% Input

U_inf = 1;  % Velocità all'infinito [m/s]
alpha = 5;   % Angolo di incidenza [°]
U_inf_x = U_inf * cos(deg2rad(alpha));
U_inf_y = U_inf * sin(deg2rad(alpha));

U_inf = [U_inf_x; U_inf_y];
U_inf_normal = [-U_inf(2); U_inf(1)];
U_inf_normal = U_inf_normal ./ norm(U_inf_normal);

TestCase = 0;

NomeProfilo = 'NACA_0012';
Chord = 1;


LE_X_Position = 0;
LE_Y_Position = 0;

%% Creazione profilo

% numero profilo:
% [x,y]=createProfile(CodiceProfilo,NPannelli,Chord);

Corpo = importXfoilProfile(strcat(NomeProfilo,'.dat'));
% Prima flippa i vettori
x = flipud(Corpo.x);
y = flipud(Corpo.y);
Corpo.x = x.*Chord;
Corpo.y = y.*Chord;
NPannelli = length(x)-1;



figure;
plot(x, y, 'o-')
hold on

axis equal

%% Creazione di una struttura di pannelli

[Centro, Normale, Tangente, Estremo_1, Estremo_2, alpha, lunghezza, L2G_TransfMatrix, G2L_TransfMatrix] = CreaStrutturaPannelli(Corpo);
        
%% Inizializzazione matrici e vettori

% Ora che ho i pannelli, posso inizializzare la matrice ed i vettori

NCols = sum(NPannelli) + 1;
NRows = NCols;
matriceA = zeros(NRows, NCols);
TermineNoto = zeros(NRows, 1);

%% Creazione della matrice quadrata As


for i = 1:NPannelli
    index_i = i; % riga

    Centro_qui = Centro(i, :)';
    Normale_qui = Normale(i, :)';

    indexStart_colonna = 0;

        for j = 1:NPannelli
            index_j = indexStart_colonna + j;  % Colonna

            Estremo_1_qui = Estremo_1(j, :)';
            Estremo_2_qui = Estremo_2(j, :)';

            L2G_TransfMatrix_qui = squeeze(L2G_TransfMatrix(j, :, :));
            G2L_TransfMatrix_qui = squeeze(G2L_TransfMatrix(j, :, :));

            matriceA(index_i, index_j) = dot(ViSorgente(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Normale_qui);

            matriceA(index_i, sum(NPannelli)+1) = matriceA(index_i, sum(NPannelli)+1) + dot(ViVortice(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Normale_qui);


        end

end


%% Creazione delle componenti dei vettori a_v, c_s e c_v


Centro_Start = Centro(1, :)';
Tangente_Start = Tangente(1, :)';

Centro_End = Centro(end, :)';
Tangente_End = Tangente(end, :)';


b = 0;
for j = 1:NPannelli(1)

    index_j = j;

    Estremo_1_qui = Estremo_1(j, :)';
    Estremo_2_qui = Estremo_2(j, :)';
    L2G_TransfMatrix_qui = squeeze(L2G_TransfMatrix(j, :, :));
    G2L_TransfMatrix_qui = squeeze(G2L_TransfMatrix(j, :, :));

    a = dot(ViSorgente(Centro_Start, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_Start);
    b = b + dot(ViVortice(Centro_Start, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_Start);

    a = a + dot(ViSorgente(Centro_End, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_End);
    b = b + dot(ViVortice(Centro_End, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_End);


    matriceA(sum(NPannelli) + 1, index_j) = a;

end

matriceA(sum(NPannelli) + 1, sum(NPannelli) + 1) = b;



%% Creazione del termine noto

for j = 1:NPannelli

    Normale_qui = Normale(j, :)';

    index = j;

    TermineNoto(index) = - dot(U_inf, Normale_qui);
end

Tangente_1 = Tangente(1, :)';
Tangente_end = Tangente(end, :)';
TermineNoto(sum(NPannelli) + 1) = - dot(U_inf, (Tangente_1 + Tangente_end));

%% Risoluzione sistema lineare
Soluzione = linsolve(matriceA,TermineNoto);


%% Calcolo del cp e della velocità sui pannelli
U = repmat(U_inf,1,length(Centro));
Punti = Centro;

for i = 1 : length(Centro(:,1))
    for j = 1 : length(Centro(:,1))
        U(:,i) = U(:,i) + Soluzione(j).*ViSorgente(Punti(i,:)',Estremo_1(j,:)',Estremo_2(j,:)',[L2G_TransfMatrix(j,1,1) , L2G_TransfMatrix(j,1,2); L2G_TransfMatrix(j,2,1) , L2G_TransfMatrix(j,2,2)],[G2L_TransfMatrix(j,1,1) , G2L_TransfMatrix(j,1,2); G2L_TransfMatrix(j,2,1) , G2L_TransfMatrix(j,2,2)])+Soluzione(end).*ViVortice(Punti(i,:)',Estremo_1(j,:)',Estremo_2(j,:)',[L2G_TransfMatrix(j,1,1) , L2G_TransfMatrix(j,1,2); L2G_TransfMatrix(j,2,1) , L2G_TransfMatrix(j,2,2)],[G2L_TransfMatrix(j,1,1) , G2L_TransfMatrix(j,1,2); G2L_TransfMatrix(j,2,1) , G2L_TransfMatrix(j,2,2)]);
    end
end


%%%XFOIL 
% Open the file
fid = fopen('NACA0012.txt', 'r');

% Skip the first few header lines until reaching the data
fgetl(fid); % First line (NACA 0012 ...)
fgetl(fid); % Second line (Alfa = ...)
fgetl(fid); % Third line (#    x        y        Cp)

% Read the data columns as floating-point values
data = textscan(fid, '%f %f %f');

% Close the file
fclose(fid);

% Extract x, y, and Cp values from the cell array
x = data{1};
y = data{2};
Cp = data{3};

% Plot Cp against x or y, depending on the desired plot
figure;
plot(x, Cp, '-o'); % Cp vs x (often useful for airfoil pressure coefficient analysis)
hold on
title('Pressure Coefficient (Cp) along the airfoil surface');
xlabel('x');
ylabel('Cp');
set(gca, 'YDir', 'reverse'); % Flip the y-axis for Cp convention (more negative is up)
grid on;


CP = 1 - ((sqrt(U(1,:).^2+U(2,:).^2))./norm(U_inf)).^2;

plot(Punti(:,1),CP,'LineWidth',2)

legend("XFOIL","HessSmith")

[~, i] = min(Punti(:,1));

figure
plot(Punti(1:i,1),-CP(1:i))
hold on
plot(Punti(i:end,1),-CP(i:end))
legend("Ventre","Dorso")


% CL Dalla circolazione

CL_1 = 2*sum(lunghezza.*Soluzione(end))/(norm(U_inf));

% CL integrando CP

U_inf_normal_CL = repmat(U_inf_normal,1,length(Normale));

CL_2 = -sum(lunghezza.*CP.*dot(Normale',U_inf_normal_CL))/Chord;

% CM 
Centro = [Centro,zeros(length(Centro),1)];
Normale = [Normale,zeros(length(Normale),1)];
z = repmat([0 0 1],length(Normale),1);

CM = sum(lunghezza.*CP.*dot(cross(Centro-[0.25,0,0],Normale)',z'));


disp(["CL CP: " num2str(CL_1)])
disp(["CL GAMMA: ", num2str(CL_2)])
disp(["CM: ", num2str(CM)] )


%% Streamlines 

x_points = linspace(-1,2,1001);
y_points = linspace(-1,1,1001);
u = zeros(length(x_points),length(y_points));
v = zeros(length(x_points),length(y_points));
u_norm = zeros(length(x_points),length(y_points));
Cp = zeros(length(x_points),length(y_points));

for i = 1 : length(x_points)
     for j = 1 : length(y_points)
        if ~(inpolygon(x_points(i),y_points(j),x,y))
            U = U_inf;
            for m = 1 : length(Centro(:,1))
                 U = U + Soluzione(m).*ViSorgente([x_points(i);y_points(j)],Estremo_1(m,:)',Estremo_2(m,:)',[L2G_TransfMatrix(m,1,1) , L2G_TransfMatrix(m,1,2); L2G_TransfMatrix(m,2,1) , L2G_TransfMatrix(m,2,2)],[G2L_TransfMatrix(m,1,1) , G2L_TransfMatrix(m,1,2); G2L_TransfMatrix(m,2,1) , G2L_TransfMatrix(m,2,2)])+Soluzione(end).*ViVortice([x_points(i);y_points(j)],Estremo_1(m,:)',Estremo_2(m,:)',[L2G_TransfMatrix(m,1,1) , L2G_TransfMatrix(m,1,2); L2G_TransfMatrix(m,2,1) , L2G_TransfMatrix(m,2,2)],[G2L_TransfMatrix(m,1,1) , G2L_TransfMatrix(m,1,2); G2L_TransfMatrix(m,2,1) , G2L_TransfMatrix(m,2,2)]);
                 u(i,j)=U(1);
                 v(i,j)=U(2);
                 u_norm(i,j)=norm(U);
                 Cp(i,j) = 1 - norm(U)/norm(U_inf);
            end
        end
     end
end


figure
plot(x, y)
hold on
contourf(x_points,y_points,Cp',100,'LineStyle','None')
streamslice(x_points,y_points,u',v',10,'noarrows');
