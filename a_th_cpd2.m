clear;
clc;
close all

NACA = '6409';
alpha_vett = 3:0.1:4;
n = length(alpha_vett);
numNodes = 160; 

cp_struct = cell(n, 2);
d2_struct = cell(n, 1);

n_points_d2 = 10; % numero  di punti da interpolare per calcolare la derivata seconda

%%

% write input file
fid = fopen('xfoil_input.txt','w');
fprintf(fid,['NACA ' NACA '\n']);
fprintf(fid,'PPAR\n');
fprintf(fid,['N ' num2str(numNodes) '\n']);
fprintf(fid,'\n\n');
fprintf(fid,'OPER\n');

for i = 1:n
    alpha = alpha_vett(i);
    Cp_fileName = ['Cp' num2str(i) '.dat'];
    if (exist(Cp_fileName,'file'))
        delete(Cp_fileName);
    end
    fprintf(fid,['Alfa ' num2str(alpha) '\n']);
    fprintf(fid,['CPWR ' Cp_fileName '\n']);
end

fclose(fid);

% run xfoil
xfoil_path = fullfile('XFOIL\xfoil.exe');
cmd = [xfoil_path, ' < xfoil_input.txt > output.txt'];
system(cmd);

% collect data from files and compute second derivative
for i = 1:n
    alpha = alpha_vett(i);
    Cp_fileName = ['Cp' num2str(i) '.dat'];
    cp_table = import_Cp_plot(Cp_fileName);
    delete(Cp_fileName);

    i_separ = 1;
    for j = 2:numNodes
        if cp_table.x(j) > cp_table.x(j-1)
            i_separ = j;
            break;
        end
    end

    x_top = cp_table.x(1:i_separ-1);
    x_top = flipud(x_top);
    x_bot = cp_table.x(i_separ:end);
    cp_top = cp_table.cp(1:i_separ-1);
    cp_top = flipud(cp_top);
    cp_bot = cp_table.cp(i_separ:end);
    
    cp_struct{i,1} = [x_top, cp_top];
    cp_struct{i,2} = [x_bot, cp_bot];
    
    [~, d2dx] = der2(x_top, -cp_top, n_points_d2);
    d2_struct{i} = [x_top, d2dx'];
end    
    
if (exist('xfoil_input.txt','file'))
    delete('xfoil_input.txt');
end

if (exist(Cp_fileName,'file'))
    delete(Cp_fileName);
end

%% plot

figure(1)
legend_entries = cell(n,1);
for i=1:n
    plot(d2_struct{i}(:,1), d2_struct{i}(:,2))
    hold on
    legend_entries{i} = ['alpha = ' num2str(alpha_vett(i))];
end
plot([0 1], [0 0], 'k')
legend(legend_entries)
ylim([-20 20])
title(['NACA' NACA]);