%% Question 4. SIR simulation
fsz = 16;
n = 1e4;
a = 2.2;
T = 0.4;
runs = 100;

figure(2);clf;
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
end   
    
set(gca,'Fontsize',fsz);
xlabel('Timesteps','Fontsize',fsz);
ylabel('Number of Infecting Nodes','Fontsize',fsz);
hold off;
