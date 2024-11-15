function [alpha, top_Xtr, bot_Xtr, n_notConverged] = transizione(NACA, Re, alpha_vett, max_iter, numNodes)

data_fileName = 'data.txt'; 

alfa_min = alpha_vett(1);
alfa_max = alpha_vett(end);
alfa_step = (alpha_vett(end) - alpha_vett(1))/(length(alpha_vett)-1);

if (exist(data_fileName,'file'))
    delete(data_fileName);
end

% run xfoil
fid = fopen('xfoil_input.txt','w');
fprintf(fid,['NACA ' NACA '\n']);
fprintf(fid,'PPAR\n');
fprintf(fid,['N ' num2str(numNodes) '\n']);
fprintf(fid,'\n\n');
fprintf(fid,'OPER\n');
fprintf(fid,['VISC ' num2str(Re) '\n']);
fprintf(fid,['ITER ' num2str(max_iter) '\n']);
fprintf(fid,'PACC\n');
fprintf(fid,[data_fileName '\n\n']);
fprintf(fid,['ASEQ ' num2str(alfa_min) ' ' num2str(alfa_max) ' ' num2str(alfa_step) '\n']);
fprintf(fid,'PACC\n\n');

fclose(fid);

xfoil_path = fullfile('..', 'xfoil.exe');
cmd = [xfoil_path, ' < xfoil_input.txt > output.txt'];
system(cmd);

 % read data file
f_id = fopen(data_fileName);                                                  % Open file for reading
dataBuffer = textscan(f_id,'%f %f %f %f %f %f %f','HeaderLines',12,...                     % Ready data from file
                            'CollectOutput',1,...
                            'Delimiter','');
fclose(f_id);                                                              % Close file
% delete(data_fileName);

% Separate data
alpha  = dataBuffer{1,1}(:,1);
top_Xtr  = dataBuffer{1,1}(:,6);
bot_Xtr  = dataBuffer{1,1}(:,7);

n_notConverged = length(alpha_vett) - length(alpha);


%%

% plot(alpha, top_Xtr, alpha, bot_Xtr);
% legend('top', 'bottom');
% xlabel('alpha');
% ylabel('x tr');
% ylim([0 1])