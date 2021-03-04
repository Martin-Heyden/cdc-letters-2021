%% No dist.
nbr_nodes = 4
N = nbr_nodes;
edges = [2:nbr_nodes;1:nbr_nodes-1]'
edge_delay = [2 4 2];


producers = [1 2 3 4]
tau = edge_delay;
r_vec = rand(1,nbr_nodes);
q_vec = rand(1,nbr_nodes);

[ A,B,Q,R ] = generate_graph(edges,producers, q_vec, r_vec,edge_delay );
nbr_states = length(A);
[Xbig,L,G,REPORT] = dare(A,B,Q,R);
G_m = generate_controller(edges, edge_delay, q_vec, r_vec,10);
G_m+G

%% With dist (requires CVX)
edges = [2:nbr_nodes;1:nbr_nodes-1]'
edge_delay = [2 4 2];
nbr_edges = length(edge_delay)
nbr_nodes = nbr_edges+1
N = nbr_nodes;
H = 10
producers = 1:N
tau = edge_delay;
r_vec = rand(1,nbr_nodes);
q_vec = rand(1,nbr_nodes);
[G_m, K_D, gamma_N] = generate_controller(edges, edge_delay, q_vec, r_vec,H);

[ A,B,Q,R ] = generate_graph(edges,producers, q_vec, r_vec,edge_delay );
nbr_states = length(A);
x0 = rand(nbr_states,1);
T = 300;

dist = zeros(nbr_nodes,T);
dist(:,1:H) = 10*randn(N,H);

cvx_begin
    cvx_precision high
    variable x(nbr_states,T+1)
    variable u(nbr_nodes + nbr_edges,T)
    obj = 0;
    for ti = 1:T
        obj = obj + x(:,ti+1)'*Q*x(:,ti+1) + u(:,ti)'*R*u(:,ti);
    end
    minimize obj
    subject to
        x(:,1) == x0;
        for ti = 1:T
            x(:,ti+1) == A*x(:,ti)+B*(u(:,ti)) + [dist(:,ti);zeros(sum(edge_delay),1)];
        end
cvx_end



[v_g,u_g] = calculate_inputs(G_m, K_D, gamma_N/q_vec(N), x(:,1), dist, H,tau);
er = u(:,1) - [v_g; u_g]
[v_g,u_g] = calculate_inputs(G_m, K_D, gamma_N/q_vec(N), x(:,2), dist(:,2:end), H,tau);
er = u(:,2) - [v_g; u_g]