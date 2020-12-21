fsz = 16;
% power law degree distribution p_k = k^(-a)/zeta(a) 
a=2.2;
G0 = @(x)polylog(a,x)/polylog(a,1);
G1 = @(x)polylog(a-1,x)./(x*polylog(a-1,1));

%% (a) Giant Component Size
% u = fzero(@(x)G1(x)-x,0.2); % 0.1963
% S = 1 - G0(u) % 0.8622

%% (b) Epidemic Fraction
% T = 0.4;
% u = fzero(@(x)G1(1-T+T*x)-x,0.3); % 0.2883
% S = 1 - G0(1-T+T*u) % 0.4078

%% (c) Least Vaccinate Rate
T = 0.4; 
nt = 100;
vs = linspace(0,1,nt); % transimissibility
u = zeros(nt,1);
S = zeros(nt,1);
for i = 1 : nt
    v = vs(i);
    u(i) = fzero(@(x) G1((1-T+T*x)*(1-v)+v)-x,0.3);
    S(i) = 1 - (G0((1-T+T*u(i))*(1-v)+v));
end
figure(1);clf;
hold on;
plot(vs,u,'Linewidth',2) 
plot(vs,S,'Linewidth',2) 
legend('u','S'); 
xlabel('v','Fontsize',fsz); 
set(gca,'Fontsize',fsz);
hold off;
