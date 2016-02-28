function twowell_sgnht_plot_euler(epsl)
% Euler Integrator

randn('seed', 6);


T = 300000;
B2 = 1;
xi = B2;

q = 0;
p = randn;
E = p*p/2 + twowell(q);
K = p*p/2;
U = twowell(q);

record_p(1) = p;
record_q(1) = q;
record_e(1) = E;
record_k(1) = K;
record_u(1) = U;
record_xi(1) = xi;

for t = 2: T
    
    p = record_p(t-1);
    q = record_q(t-1);
    xi = record_xi(t - 1);
        
    p = p - xi*p*epsl - grad_twowell(q) * epsl - randn* sqrt(2*B2*epsl);
    q = q + epsl * p;
    xi = xi + (p^2 - 1) * epsl;
    
    E = p*p/2 + twowell(q);
    K = p*p/2;
    U = twowell(q);
    
    record_p(t) = p;
    record_q(t) = q;
    record_e(t) = E;
    record_k(t) = K+record_k(t-1);
    record_u(t) = U;
    record_xi(t) = xi;
end

% true plot
xStep = 0.01;
xx = [-5:xStep:5];  L = length(xx);
for i = 1: L
    yy(i) = exp(-twowell(xx(i)));
end

fun = @(x) exp((-x.^4 - x.^3 + 13 * x.^2 + x - 12) / 14 - 0.5);
q = quad(fun, -10, 20);
yy = yy / q;

figure(1)
clf; title(['h = ', num2str(epsl) ]); 
hold on; axis([-5,4,0,1]); grid on;
[f,x] = hist(record_q,L); ye = f/trapz(x,f);
b1 = bar(x,ye,'FaceColor',[.9 .9 .5],'EdgeColor',[.9 .9 .5]); 
set(get(b1,'Children'),'FaceAlpha',.95)
plot (xx, yy,'k','linewidth',2); hold on;


hold off

figure(2)
clf;
hold on
%plot(record_e,'b');
%plot(record_k./[1:T],'g-.','linewidth',2);
%plot(record_u,'r');
record_xi = cumsum(record_xi);
plot(record_xi./[1:T],'g','linewidth',2);
hold off
end

function l = twowell(x)
    l = (x+4)*(x+1)*(x-1)*(x-3)/14 + 0.5;
end


function g = grad_twowell(x)
    g = (x+1)*(x-1)*(x-3)/14 + (x+4)*(x-1)*(x-3)/14 + (x+4)*(x+1)*(x-3)/14 + (x+4)*(x+1)*(x-1)/14;
end