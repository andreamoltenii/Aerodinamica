clear;
clc;
addpath functions 

NACA = '2414';

numNodes = 160;

Re_vett = 8e5:2e5:1.4e6;

alpha_vett = 0:1:5;

max_iter = 100;

l_re = length(Re_vett);
l_alfa  = length(alpha_vett);

x_sep_top = zeros(l_alfa, l_re);
x_sep_bot = zeros(l_alfa, l_re);
x_riatt_top = zeros(l_alfa, l_re);
x_riatt_bot = zeros(l_alfa, l_re);
x_tr_top = zeros(l_alfa, l_re);
x_tr_bot = zeros(l_alfa, l_re);
n_conv = zeros(1, l_re);

for i = 1:l_re
    for j = 1:l_alfa
        [x_sep_top(j,i), x_riatt_top(j,i), x_sep_bot(j,i), x_riatt_bot(j,i)] = separazione(NACA, Re_vett(i), alpha_vett(j), max_iter, numNodes);
    end
    [alpha_conv, x_top, x_bot] = transizione(NACA, Re_vett(i), alpha_vett, max_iter, numNodes);
    n_conv(i) = length(alpha_conv);
    x_tr_top(1:n_conv(i),i) = x_top;
    x_tr_bot(1:n_conv(i),i) = x_bot;
end

%% plot transizione

figure(1)
legend_entries = cell(l_re,1);
for i = 1:l_re
    plot(alpha_vett(1:n_conv(i)), x_tr_top(1:n_conv(i), i));
    hold on
    legend_entries{i} = ['Re = ' num2str(Re_vett(i))];
end
legend(legend_entries)
xlabel('alpha');
ylabel('x transizione sul dorso')
title(['NACA' NACA]);

figure(2)
legend_entries = cell(l_re,1);
for i = 1:l_re
    plot(alpha_vett(1:n_conv(i)), x_tr_bot(1:n_conv(i), i));
    hold on
    legend_entries{i} = ['Re = ' num2str(Re_vett(i))];
end
legend(legend_entries)
xlabel('alpha');
ylabel('x transizione sul ventre')
title(['NACA' NACA]);

%% plot separazione (fa schifo)

figure(3)
legend_entries = cell(l_re,1);
flag = 0;
for i = 1:l_re
    if ~isequal(x_sep_top(:, i), -ones(l_alfa, 1))
        plot(alpha_vett, x_sep_top(:, i), 'o');
        hold on
        legend_entries{i} = ['Re = ' num2str(Re_vett(i))];
        flag = 1;
    end
end
if flag == 1
    legend(legend_entries)
    xlabel('alpha');
    ylabel('x separazione sul dorso')
    ylim([0 1]);
    title(['NACA' NACA]);
else
    close(3)
end
