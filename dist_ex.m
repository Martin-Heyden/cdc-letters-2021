rng(250212)
%Define graph
nbr_nodes = 5
N = nbr_nodes;
edges = [2:nbr_nodes;1:nbr_nodes-1]'
edge_delay = [3 2 5 4];
nbr_edges = length(edge_delay)
producers = 1:N

%Weightings and H
r_vec = 100*ones(1,nbr_nodes);
q_vec = rand(1,nbr_nodes)*0.4+0.4;
H = 10
%Generate state space
[ A,B,Q,R ] = generate_graph(edges,producers, q_vec, r_vec,edge_delay );
nbr_states = length(A);

%Iniital conditions
x0 = zeros(nbr_states,1);
%HOrizon
T = 100;

%Generate disturbances
d = zeros(nbr_nodes,T);
d(3,10:13) = -0.5;

d2 = zeros(nbr_nodes,T);
d2(2,12:15) = -0.3;

%Generate optimal controller
[Kx, Kd, gamma_N] = generate_controller(edges, edge_delay, q_vec, r_vec,H);
%%% With feedforward %%%
x = zeros(nbr_states,T+1);
for ti = 1:T
    [v,u] = calculate_inputs(Kx, Kd, gamma_N/q_vec(N), x(:,ti), d(:,ti:end)+d2(:,ti:end), H,edge_delay)
    x(:,ti+1) = A*x(:,ti)+B*[v;u]+ [d(:,ti) + d2(:,ti);zeros(sum(edge_delay),1)];
end


%%% Without feedforward %%%
z = zeros(nbr_states,T+1)
for ti = 1:T
    z(:,ti+1) = A*z(:,ti)+B*Kx*z(:,ti)+ [d(:,ti) + d2(:,ti);zeros(sum(edge_delay),1)];
end

%Plotitng
figure(1)
clf
h = 90
h1 = h+1
subplot(2,1,1)
plot(0:h,x(1:5,1:h1),'Linewidth',3)
title('With feedforward','FontSize', 16)
ylim([-1.25,0.5])
xlabel('Samples','FontSize', 14)
ylabel('Node Levels','FontSize', 14)
subplot(2,1,2)
plot(0:h,z(1:5,1:h1),'Linewidth',3)
title('No feedforward','FontSize', 16)
legend('1','2','3','4','5')
ylim([-1.25,0.5])
xlabel('Samples','FontSize', 14)
ylabel('Node Levels','FontSize', 14)

