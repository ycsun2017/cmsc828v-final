%% Question 4. SIR simulation
fsz = 16;
% n = 1e4;
a = 2.2;
T = 0.4;
runs = 100;

% ns = [10,50,1e2,5e2,1e3,2e3,4e3,6e3,8e3,1e4];
% mean_time = zeros(1,length(ns));
ns = [1e3,5e3,1e4,2e4];
for nd = 1 : length(ns)
    n = ns(nd)
    durations = zeros(1,runs);
    figure(nd);clf;
    hold on;
    for run = 1 : runs
        % initialize the graph
        [G,edges,K,p] = MakePowerLawRandomGraph(n,a);

        % infecting nodes
        inf_nodes = [];
        start = randi(n);
        inf_nodes(end+1) = start;

        % recovered nodes
        recovered = sparse(1,n);

        % duration: number of timesteps that infecting nodes exist
        duration = 0;

        % statistics
        num_infecting = [];

        while ~isempty(inf_nodes)
            duration = duration + 1;
            num_infecting(end+1) = length(inf_nodes);
            new_inf = [];
            for idx = 1 : length(inf_nodes)
                % one node can infect neighbors with probability T, then he/she can
                % recover from the disease
                [~,cols,~] = find(G(inf_nodes(idx),:));
                for j = 1 : length(cols)
                    neighbor = cols(j);
                    if ~recovered(neighbor) && rand() < T && ~any(ismember(new_inf, neighbor))
                        % infect new
                        new_inf(end+1) = neighbor;
                    end
                end
                recovered(inf_nodes(idx)) = 1;
            end
            inf_nodes = new_inf;
        end
        plot(1:duration, num_infecting./n);
        durations(run) = duration;
    end   

    set(gca,'Fontsize',fsz);
    xlabel('Timesteps','Fontsize',fsz);
    ylabel('Fraction of Infecting Nodes','Fontsize',fsz);
    hold off;

    fprintf("The mean duration is %d\n", mean(durations))
    mean_time(nd) = mean(durations);
end

% figure(1);clf;
% hold on;
% plot(ns, mean_time);
% set(gca,'Fontsize',fsz);
% xlabel('n','Fontsize',fsz);
% ylabel('Mean duration','Fontsize',fsz);
% hold off;
