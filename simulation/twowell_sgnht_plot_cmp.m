%% compare different SGNHT-E and SGNHT-S approaches

%% SGNHT-E

B2 = 1; xi = B2;
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


[yhmc,xhmc] = hist(record_q, xGrid);
yhmc = yhmc / sum(yhmc) / xStep;

error = sum(y.*(log(y./(yhmc + 1e-60))))/2 + sum((yhmc + 1e-60).*(log((yhmc + 1e-60)./(y + 1e-60))))/2;
ERRORs(1,i) = error;
disp(['SGNHT-E KL: ' num2str(error)]);


%% SGNHT-S
B2 = 1; xi = B2;
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
        
    q  = q + p * epsl/2;
    xi = xi + (p^2 - 1) * epsl/2;
    p  = exp(-xi*epsl/2) * p -  grad_twowell(q) * epsl/2; 
    p  = p +  randn* sqrt(2*B2*epsl);
    p  = exp(-xi*epsl/2) * p -  grad_twowell(q) * epsl/2; 
%     p  = exp(-xi*epsl/2) * p;
%     p  = p +  randn* sqrt(2*B2*epsl) -  grad_twowell(q) * epsl; 
%     p  = exp(-xi*epsl/2) * p;    
    xi = xi + (p^2 - 1) * epsl/2;
    q  = q + p * epsl/2;
    
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


[yhmc,xhmc] = hist(record_q, xGrid);
yhmc = yhmc / sum(yhmc) / xStep;

error = sum(y.*(log(y./(yhmc + 1e-60))))/2 + sum((yhmc + 1e-60).*(log((yhmc + 1e-60)./(y + 1e-60))))/2;
ERRORs(2,i) = error;
disp(['SGNHT-S KL: ' num2str(error)]);


%% Plot functions and its estimations
% legend( 'True Distribution', 'SGLD',  'SGLD\_rmsprop' );
% 
% len = 5;
% set(gcf, 'PaperPosition', [0 0 len len/8.0*6.5] )
% set(gcf, 'PaperSize', [len len/8.0*6.5] )
% xlabel('\theta');
% saveas( gcf, fgname, 'fig');


%% Plot errors with different #samples
% nc=nsample/10; np = nsample/nc; ERRORs = zeros(4, np);
% for p = 1:np
%     for i = 1:2
%         [yhmc,xhmc] = hist(ys(i, 1: p*nc ), xGrid);
%         yhmc = yhmc / sum(yhmc) / xStep;
%         ERRORs(i,p) = sum(y.*(log(y./(yhmc + 1e-16))));
%         % ERRORs(i,p) = sum(abs(y-yhmc))/length(y);
%     end
% end

