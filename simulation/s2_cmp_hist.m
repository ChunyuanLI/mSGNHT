

clear all; clc;

% stepsiz
epsl =  0.2;

twowell_sgnht_plot_euler(epsl);
twowell_sgnht_plot_split(epsl);

figure(1)
legend({'SGNHT\_E', 'true distribution'},'FontSize', 14);

figure(3)
legend({'SGNHT\_S', 'true distribution'},'FontSize', 14);

figure(2)
legend({'\xi: SGNHT\_E', '\xi: SGNHT\_S'},'FontSize', 14);
