function [giantsize] = findgiantsize(A, n)
% find the size of giant component for a given adjacent matrix A with n nodes

comp = zeros(1,n);
visited = [];
giantsize = 0;
n_comp = 0;

for i = 1 : n
    if comp(i) == 0 % node i's component is not determined yet
        n_comp = n_comp + 1;
        v = [];
        v(1) = i;
        sz = 0;
        % traverse nodes in the list and assign all reachable new nodes to the current component 
        while nnz(v) > 0
            comp(v(1)) = n_comp;
            visited(end+1) = v(1);
            neighbors = find(A(v(1),:));
            new = setdiff(neighbors,visited);
            v = union(v(2:end),new); 
            sz = sz + 1;
        end
        if sz > giantsize
            giantsize = sz;
        end
    end
end

end


