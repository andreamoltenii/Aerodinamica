function [x_sep_top, x_riatt_top, x_sep_bot, x_riatt_bot] = separazione(NACA, Re, alpha, max_iter, numNodes, plot_flag)

% runna xfoil con naca, re, alpha, max_iter e numNudes specificati come
% parametri
% restituisce le coordinate di separzione e riattacco su dorso e ventre (in
% base al segno del Cf) se avvengono, altrimenti le corrispondenti 
% variabili avranno valore -1

if(nargin == 5)
    plot_flag = 0;
end

data_fileName = 'data.txt'; 
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
fprintf(fid,['ALFA ' num2str(alpha) '\n']);
fprintf(fid, 'VPLO\n');
fprintf(fid, 'CF\n');
fprintf(fid,['DUMP ' data_fileName '\n']);

fclose(fid);

xfoil_path = fullfile('..', 'xfoil.exe');
cmd = [xfoil_path, ' < xfoil_input.txt > output.txt'];
system(cmd);

 % read data file
f_id = fopen(data_fileName);                                                  % Open file for reading
dataBuffer = textscan(f_id,'%f %f','HeaderLines',7,...                     % Ready data from file
                            'CollectOutput',1,...
                            'Delimiter','');
fclose(f_id);                                                              % Close file
delete(data_fileName);

% Separate data
x  = dataBuffer{1,1}(:,1);
cf  = dataBuffer{1,1}(:,2);

%%

n = length(x);

i_top_end = 1;
i_bot_start = 1;
flag = 0;
for i=1:n
    if flag == 0 && x(i) >= 1
        i_top_end  = i;
        flag = 1;
    end
    if flag == 1 && x(i) < 1
        i_bot_start = i;
        flag = 2;
    end
    if flag == 2 && x(i) >= 1
        i_bot_end = i;
        break;
    end
end
if flag == 0
    error('error')
end

%%

x_top = x(1:i_top_end);
cf_top = cf(1:i_top_end);
n_top = length(x_top);

x_bot = x(i_bot_start:i_bot_end);
cf_bot = cf(i_bot_start:i_bot_end);
n_bot = length(x_bot);

%%

flag_sep = 0;
x_sep_top = -1;
x_riatt_top = -1;
for i=1:n_top
    if flag_sep == 0 && cf_top(i) <= 0
        flag_sep = 1;
        x_sep_top = x_top(i);
    end
    if flag_sep == 1 && cf_top(i) > 0
        x_riatt_top = x_top(i);
        break;
    end
end

flag_sep = 0;
x_sep_bot = -1;
x_riatt_bot = -1;
for i=1:n_bot
    if flag_sep == 0 && cf_bot(i) <= 0
        flag_sep = 1;
        x_sep_bot = x_bot(i);
    end
    if flag_sep == 1 && cf_bot(i) > 0
        x_riatt_bot = x_bot(i);
        break;
    end
end

if plot_flag ~= 0
    plot(x_top, cf_top);
    hold on
    plot(x_bot, cf_bot);
    plot([0 1], [0 0], 'k')
    legend('top', 'bottom')
end


