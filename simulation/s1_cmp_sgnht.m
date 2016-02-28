
function s1_cmp_sgnht()
clear all; clc;
randn('seed', 7);
fgname = 'figures/cmp_sgnht';
xStep = 0.01; T = 1000000; 

Stepsizes = [1e-3,0.01,0.05,0.1,0.15,0.2,0.25,0.3]; 
np = length(Stepsizes); ERRORs = zeros(2,np);

% true plot
xGrid = [-5:xStep:5]; nGrids = length(xGrid); 
for i = 1: nGrids
    y(i) = exp(-twowell(xGrid(i)));
end
y = y / sum(y) / xStep;

for i = 1:np
    epsl = Stepsizes(i);%0.005
    disp(['Step size: ' num2str(epsl)]);

    twowell_sgnht_plot_cmp;
end


figure();  axis manual; 

semilogy(1:np,ERRORs(1,:),'g', 'LineWidth',2); hold on;
semilogy(1:np,ERRORs(2,:),'b', 'LineWidth',2); hold on;

axis([1, np, 0, max(ERRORs(:))]); grid on;

labels = {'1e-3';'0.01';'0.05';'0.1';'0.15';'0.2';'0.25';'0.3'};
set(gca,'XTickLabel',labels)

len = 5;
set(gcf, 'PaperPosition', [0 0 len len/8.0*6.5] )
set(gcf, 'PaperSize', [len len/8.0*6.5] )

legend( 'SGNHT\_E', 'SGNHT\_S');
xlabel('stepsize'); ylabel('KL divergence');
saveas( gcf, [fgname, 'error'] , 'fig');
saveas( gcf, [fgname, 'error'] , 'pdf');


end

function l = twowell(x)
    l = (x+4)*(x+1)*(x-1)*(x-3)/14 + 0.5;
end


function g = grad_twowell(x)
    g = (x+1)*(x-1)*(x-3)/14 + (x+4)*(x-1)*(x-3)/14 + (x+4)*(x+1)*(x-3)/14 + (x+4)*(x+1)*(x-1)/14;
end