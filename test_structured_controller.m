%No dist
nbr_nodes = 4;
N = nbr_nodes;
edges = [2:nbr_nodes;1:nbr_nodes-1]';
edge_delay = [2 4 2];


producers = [1 2 3 4];
tau = edge_delay;
r_vec = rand(1,nbr_nodes);
q_vec = rand(1,nbr_nodes);

[ A,B,Q,R ] = generate_graph(edges,producers, q_vec, r_vec,edge_delay );
nbr_states = length(A);
[Xbig,L,G,REPORT] = dare(A,B,Q,R);
G_m = generate_controller(edges, edge_delay, q_vec, r_vec,10);


contr = structured_controller(edges,tau,q_vec,r_vec,1);
x = randn(nbr_states,1);
%x = [randn(N,1);[1 1]';[0 0 0 0]';[0 0]'];
[u,v] = contr.calculate_input(x,zeros(N,10));
%-G_m*x
t = G_m*x;
v - t(1:4)
u - t(5:7)
%%
x = randn(nbr_states,1);
H = 10;
dist = zeros(N,H);
dist(:,1) = randn(N,1);
%dist(:,1:H) = 10*randn(N,H);
contr = structured_controller(edges,tau,q_vec,r_vec,H);

[u,v] = contr.calculate_input(x,dist);
[G_m, K_D, gamma_N] = generate_controller(edges, edge_delay, q_vec, r_vec,H);
[v_g,u_g] = calculate_inputs(G_m, K_D, gamma_N/q_vec(N), x(:,1), dist, H,tau);
er = [v;u] - [v_g; u_g]
