fsz = 16;
%% Question 3 (a),(b) compute the average fraction of giant component and epidemic
% n = 1e4;
% a = 2.2;
% T = 0.4;
% runs = 10;
% 
% gsizes = zeros(1,runs);
% esizes = zeros(1,runs);
% 
% for run = 1 : runs
%     [G,edges,K,p] = MakePowerLawRandomGraph(n,a);
% 
%     gsize = findgiantsize(G, n);
%     gsizes(run) = gsize;
% 
%     E = sparse(n,n);
%     [i,j,s] = find(G);
%     for k = 1 : length(i)
%         if rand() < T
%             E(i(k),j(k)) = s(k);
%         end
%     end
% 
%     esize = findgiantsize(E, n);
%     esizes(run) = esize;
% end

% gfrac = mean(gsizes)/n % 0.8213
% efrac = mean(esizes)/n % 0.3655


%% Question 3 (c) find the critical transmissibility
% n = 1e4;
% a = 2.2;
% runs = 10;
% 
% nt = 100;
% t = linspace(0,1,nt);
% esizes = zeros(1,nt);
% 
% [G,edges,K,p] = MakePowerLawRandomGraph(n,a);
% 
% for iter = 1 : nt
%     E = sparse(n,n);
%     [i,j,s] = find(G);
%     for k = 1 : length(i)
%         if rand() < t(iter)
%             E(i(k),j(k)) = s(k);
%         end
%     end
%     esize = findgiantsize(E, n);
%     esizes(iter) = esize;
% end
% efracs = esizes./n;
% 
% figure(2);
% hold on;
% plot(t, efracs);
% set(gca,'Fontsize',fsz);
% xlabel('T','Fontsize',fsz);
% ylabel('Fraction','Fontsize',fsz);




%% Question 3 (d) find the critical vaccinate fraction
n = 1e4;
a = 2.2;
runs = 10;
T = 0.4;
nv = 100;
v = linspace(0,1,nv);
esizes = zeros(1,nv);

[G,edges,K,p] = MakePowerLawRandomGraph(n,a);

for iter = 1 : nv
    removed = randperm(n,round(n*v(iter)));
    E = sparse(n,n);
    [i,j,s] = find(G);
    for k = 1 : length(i)
        if any(ismember(removed,i(k))) || any(ismember(removed,j(k)))
            continue
        else
            if rand() < T
                E(i(k),j(k)) = s(k);
            end
        end
    end
    esize = findgiantsize(E, n);
    esizes(iter) = esize;
end

efracs = esizes./n;
figure(3);
hold on;
plot(v, efracs);
set(gca,'Fontsize',fsz);
xlabel('T','Fontsize',fsz);
ylabel('Fraction','Fontsize',fsz);

