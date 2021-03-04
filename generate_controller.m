function [K, K_D, gamma_N] = generate_controller(edges, tau, q_vec, r_vec, H)
    %Returns two matrices and a scalar. K is used on the statespace(levels, and
    %transit states). Can be combined to also handle all D_i for i neq N
    %(see test_generate_controller).
    %K_D is used to multiply D_N.
    %gamma_N is needed for D_N[t+1]
    %Inputs: edges: an array describing the source and destination on each
    %               edge
    %       tau: The delay on each edge
    %       q_vec: a vector of the q values for each node
    %       r_vec: A vector of the r values for each ndoe
    %       H: The horizon for disturbances.
   
    tau_extended = [tau, H];
    max_delay = max(tau_extended);
    N = length(edges)+1;%Nbr of nodes.
    nbr_states = N+sum(tau);
    
    %%%% Gamma, rho,  - requires communication  Sweep 1%%%%%%
     gamma = zeros(N,1);
     rho = zeros(N,1);
     gamma(1) = q_vec(1);
     rho(1) = r_vec(1);
     for i = 2:N
         gamma(i) = q_vec(i)*gamma(i-1)/(q_vec(i) + gamma(i-1));
         rho(i) = r_vec(i)*rho(i-1)/(r_vec(i)+rho(i-1));
     end
     gamma_N = gamma(N);
 
    
     %%%% X - Requires Communication Sweep 2 %%%%
     Xr = gamma(N)/2+sqrt(gamma(N)*rho(N)+gamma(N)^2/4);%same as dare(1,1,gamma(N),rho(N));
     X = zeros(N,max_delay);
     X(N,H+2) = Xr-gamma(N);
     for j = H+1:-1:1
        X(N,j) = (rho(N)*(gamma(N)+X(N,j+1)))/(rho(N)+gamma(N)+X(N,j+1));
     end
     for i = N-1:-1:1
         for j = tau(i):-1:1
            if j == tau(i) %If first for node i
                X(i,j) = (rho(i)*(gamma(i)+X(i+1,1)))/(rho(i)+gamma(i)+X(i+1,1));
            else
                X(i,j) = (rho(i)*(gamma(i)+X(i,j+1)))/(rho(i)+gamma(i)+X(i,j+1));
            end
         end

     end

     %%%%% g and b - very little comm %%%%%
     g = zeros(N,max_delay);
     for i = 1:N
        for j = 2:tau_extended(i)
            g(i,j) = X(i,j)/(X(i,j)+gamma(i));
        end
        if i ~= 1 %If not the first node, calculate g_i(1)
            g(i,1) = X(i,1)/(X(i,1)+gamma(i-1));
        end
     end
     b = zeros(N-1,1); %b_N not needed (and is not defined)
     for i =1:N-1
         if tau(i)>= 2
            b(i) = g(i+1,1)* prod(g(i,2:tau(i)));
         else
             b(i) = g(i+1,1);
         end
     end

     %%%%% P - Local %%%%%%
     P = zeros(N,max_delay,max_delay);
     for i = 1:N
        P(i,1,:) = X(i,1)/rho(i);
        for l = 2:tau_extended(i)
            for m = 1:tau_extended(i)
                if l<= m %then add g term        here
                    P(i,l,m) = (1-X(i,l)/rho(i))*g(i,l)*P(i,l-1,m) + X(i,l)/rho(i);
                else %Else, no g term
                    P(i,l,m) = (1-X(i,l)/rho(i))*P(i,l-1,m) + X(i,l)/rho(i);
                end
            end
        end
     end

     %%% h - Local %%%
     h = zeros(1,N-1); 
     h(1) = P(1,tau(1),tau(1))*g(2,1);
     for i = 2:N-1 %h defined for node N
        h(i) = (1-P(i,tau(i),1))*b(i)*h(i-1)+P(i,tau(i),tau(i))*g(i+1,1);
     end

     
    %%% Calculate phi %%%
    phi = zeros(N,max(tau_extended));
    for i = 1:N
        g_agg = 1; %used to calculate product of g_i
        for delta = 1:tau_extended(i)
            if i == 1 %for first node, as h_0 = 0.
                phi(i,delta) = 1-P(i,tau_extended(i),delta);
            else
                phi(i,delta) = 1-P(i,tau_extended(i),delta) ...
                    -(1-P(i,tau_extended(i),1))*h(i-1)*g_agg;
                if delta~= tau_extended(i) %no need to update on last iteration
                    g_agg = g_agg*g(i,delta+1);
                end
            end
            
        end
    end
    
    %%% calculate a_k and c_k
    a = zeros(1,N);
    c = zeros(1,N);
    for k = 2:N
        a(k) = X(k,1)/r_vec(k) +gamma(k)/q_vec(k) * (1-X(k,1)/rho(k));
        c(k) = -( X(k,1)/r_vec(k)- gamma(k)/q_vec(k)* X(k,1)/rho(k))  * (1-h(k-1)) ...
            +gamma(k)/q_vec(k)*h(k-1);
    end
    
    %%% Calculate delta  excluding D_N %%%
    delta = zeros(N,nbr_states);
    for i = 1:N-1
        if i~= 1 %For all nodes but the first, use the old value of delta
            delta(i,:) = (1-P(i,tau_extended(i),1))*delta(i-1,:);
        end
        delta(i,i) = phi(i,1);
        if i~= N %For all node except N, add u_i[t-tau_i]
            delta(i,N+1+sum(tau(1:i-1)):N+sum(tau(1:i))) = phi(i,1:tau(i));
        end
    end

    %%%Calculate mu (and lambda) excluding D_N %%%
    mu = zeros(N,nbr_states);
    mu(N,N) = 1; %D_N is handled later
    for i = N-1:-1:1
       mu(i,:) = mu(i+1,:)*b(i);
       mu(i,i) = 1;
       g_agg = 1;
       count = 2;
       for j = N+1+sum(tau(1:i-1)):N+sum(tau(1:i))
           mu(i,j) = g_agg;
           if count<=tau(i)
               g_agg = g_agg * g(i,count);
               count = count+1;
           end
       end
    end
    
    %%% Calculate the effect of D_N on epsilon %%%
    epsilon_D = zeros(N,H);
    epsilon_D(N,:) = phi(N,1:H); %only affect node N
    
    %%% Calculate the effect of D_N on lambda %%%
    lambda_D = zeros(N,H);
    g_agg = 1;
    count = 2;
    for j = 1:H %Loop over horion for top node
        lambda_D(N,j) = g_agg;
        if count<=H
            g_agg = g_agg * g(N,count);
            count = count+1;
        end
    end
    for i = N-1:-1:1 %then loop through graph
        lambda_D(i,:) = lambda_D(i+1,:)*b(i);
    end       
        
    %%% Valculate v %%%
    K = zeros(2*length(edges)+1,nbr_states);
    K_D = zeros(2*length(edges)+1,H);
    K(1,:) = -X(1,1)/rho(1)*mu(1,:); %h_0 = 0, delta_0 = 0.
    K_D(1,:) = -X(1,1)/rho(1)*lambda_D(1,:);
    for k = 2:N
        K(k,:) = -X(k,1)/r_vec(k)*(...
             delta(k-1,:)...
             +(1-h(k-1))*mu(k,:));
         K_D(k,:) = -X(k,1)/r_vec(k)*(...
             epsilon_D(k-1,:)...
             +(1-h(k-1))*lambda_D(k,:));
    end
    
    %%% Calculate u_{k-1} %%%
    for k = 2:N
        %zu:The transportation state that will arrive in one timeunite, and the current level
        zu = zeros(1,nbr_states);
        zu(k) = 1;
        if k ~= N %If k == N then no incoming transportation.
            zu(N+sum(tau(1:k-1))+1) = 1;
        end
        K(N+k-1,:) = (1-gamma(k)/q_vec(k))*zu + ...
            -a(k)*delta(k-1,:) ...
            +c(k)*mu(k,:);
        K_D(N+k-1,:) = -a(k)*epsilon_D(k-1,:)+c(k)*lambda_D(k,:);
    end
    

end