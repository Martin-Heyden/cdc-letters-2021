rng(2021)
%%% Settings %%%
N_vec =[10,20]
tau_vec = [3 6 12]
nbr_comb = length(N_vec)*length(tau_vec)
max_horizon = 150;
simulation_horizon = 350;
number_of_simulations = 50;

step_limit = [10,30,50,100,200];
step =[1,3,5,10,20];
cur = 0;
horizon_vec = []; 
while cur <max_horizon %Similar to logspace
    horizon_vec = [horizon_vec cur];
    cur = cur + step(find(cur<step_limit,1,'first'))
end
%horizon_vec = min_horizon:horizon_step:max_horizon

%%% Data structures %%%
cost_vec = zeros(length(N_vec),length(tau_vec),length(horizon_vec));
cost_no_F = zeros(length(N_vec),length(tau_vec));
case_id = 0;

for N_ind = 1:length(N_vec)
    N = N_vec(N_ind)
    for i = 1:number_of_simulations
        %%% Generate disturbances %%%
        %Pick 2 nodes at random
        dist_nodes = randsample(1:N,2);
        %randomize total disturbance
        dist_size = -rand(1,2);

        %Randomize start and end time
        dist_start = max_horizon + randi(30,1,2);
        dist_length = randi(5,1,2);

        d = zeros(N,simulation_horizon);
        for j = 1:2
            d(dist_nodes(j),dist_start(j):dist_start(j)+(dist_length(j)-1)) =  ...
                dist_size(j)/dist_length(j);
        end
        %Loop over different delays (same disturbance for all cases)
        for tau_ind = 1:length(tau_vec)       
            tau_val = tau_vec(tau_ind);
            %Generate graph
            edges = [2:N;1:N-1]';
            tau = tau_val*ones(1,N-1);
            r_vec = 10*ones(1,N)*N;
            q_vec = ones(1,N);
            [ A,B,Q,R ] = generate_graph(edges,1:N, q_vec, r_vec,tau );
            nbr_states = length(A);
            %%% Generate controller %%%
            H = horizon_vec(end);
            [Kx, Kd, gamma_N] = generate_controller(edges, tau, q_vec, r_vec,H);
            %%% No horizon %%%
            x = zeros(nbr_states,1);
            cost = 0;
            for ti = 1:simulation_horizon
                u = Kx*x;
                cost = cost + x'*Q*x + u'*R*u;
                x = A*x+B*u+ [d(:,ti);zeros(sum(tau),1)];
            end
            cost_no_F(N_ind,tau_ind) = cost_no_F(N_ind,tau_ind)+cost/number_of_simulations;
            %%% for all horizons
            for h_ind = 1:length(horizon_vec)
                h = horizon_vec(h_ind);
                cost = 0;
                x = zeros(nbr_states,1);
                for ti = 1:simulation_horizon
                    [v,u] = calculate_inputs(Kx, Kd, gamma_N/q_vec(N), x, d(:,ti:min(simulation_horizon,ti+h)), H,tau);
                    cost = cost + x'*Q*x + [v; u]'*R*[v;u];
                    x = A*x+B*[v;u]+ [d(:,ti);zeros(sum(tau),1)];
                end
                cost_vec(N_ind,tau_ind,h_ind) = cost_vec(N_ind,tau_ind,h_ind) + cost/number_of_simulations;
            end
        end
        
    end
end
%% Plotting
clf
subplot(2,1,1)
hold on
%cost_no_F(1,1)
ones_vec = ones(size(horizon_vec))
c1 = [0 0.4470 0.7410];
c2 = [0.8500 0.3250 0.0980];
c3 = [0.9290 0.6940 0.1250];
c4 = [0.4940 0.1840 0.5560];
ms = 5; %marker size

plot([horizon_vec]+1,squeeze(cost_vec(1,1,:))','x','Color',c1,'Linewidth',2,'MarkerSize',ms)
l1 = "N = " + int2str(N_vec(1)) + "  \tau_i = " + int2str(tau_vec(1))
hold on
plot([horizon_vec]+1,squeeze(cost_vec(1,2,:))','x','Color',c2,'Linewidth',2,'MarkerSize',ms)
l2 = "N = " + int2str(N_vec(1)) + "  \tau_i = " + int2str(tau_vec(2))
plot([horizon_vec]+1,squeeze(cost_vec(1,3,:))','x','Color',c3,'Linewidth',2,'MarkerSize',ms)
l3 = "N = " + int2str(N_vec(1)) + "  \tau_i = " + int2str(tau_vec(3))
legend({l1,l2,l3},'FontSize', 12)
% Not using disturbances at all
plot(-1,cost_no_F(1,1),'x','Color',c1,'HandleVisibility','off','Linewidth',2)
plot(-1,cost_no_F(1,2),'x','Color',c2,'HandleVisibility','off','Linewidth',2)
plot(-1,cost_no_F(1,3),'x','Color',c3,'HandleVisibility','off','Linewidth',2)
xlabel('Horizon Length','FontSize', 14)
ylabel('Cost','FontSize', 14)
subplot(2,1,2)
hold on
plot([horizon_vec]+1,squeeze(cost_vec(2,1,:))','x','Linewidth',2,'MarkerSize',ms)
l1 = "N = " + int2str(N_vec(2)) + "  \tau_i = " + int2str(tau_vec(1))
hold on
plot(horizon_vec+1,squeeze(cost_vec(2,2,:))','x','Linewidth',2,'MarkerSize',ms)
l2 = "N = " + int2str(N_vec(2)) + "  \tau_i = " + int2str(tau_vec(2))
cost_no_F(2,3)
plot([horizon_vec]+1,squeeze(cost_vec(2,3,:))','x','Linewidth',2,'MarkerSize',ms)
l3 = "N = " + int2str(N_vec(2)) + "  \tau_i = " + int2str(tau_vec(3))
legend({l1,l2,l3},'FontSize', 12)
% Not using disturbances at all
plot(-1,cost_no_F(2,1),'x','Color',c1,'HandleVisibility','off','Linewidth',2)
plot(-1,cost_no_F(2,2),'x','Color',c2,'HandleVisibility','off','Linewidth',2)
plot(-1,cost_no_F(2,3),'x','Color',c3,'HandleVisibility','off','Linewidth',2)
xlabel('Horizon Length','FontSize', 14)
ylabel('Cost','FontSize', 14)
