classdef structured_controller
	
	properties
		q
		r
		tau
		H
		N
		gamma
		X
		g
		P
		h
		phi
		b
		a
		c
		
	end
	
	methods
		%Corresponds to algorithm 1
		function obj = structured_controller(edges, tau, q_vec, r_vec, H)
			obj.q = q_vec;
			obj.r = r_vec;
			tau_extended = [tau, H];
			obj.tau = tau;
			obj.H = H;
			obj.N = length(tau_extended);
			max_delay = max(tau_extended);
			N = length(edges)+1;%Nbr of nodes.
			
			%%%% Gamma, rho,  - requires communication  Sweep 1%%%%%%
			gamma = zeros(N,1);
			rho = zeros(N,1);
			gamma(1) = q_vec(1);
			rho(1) = r_vec(1);
			for i = 2:N
				gamma(i) = q_vec(i)*gamma(i-1)/(q_vec(i) + gamma(i-1));
				rho(i) = r_vec(i)*rho(i-1)/(r_vec(i)+rho(i-1));
			end
			obj.gamma = gamma;
			
			
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
			obj.X = X;
			
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
			obj.g = g;
			b = zeros(N-1,1); %b_N not needed (and is not defined)
			for i =1:N-1
				if tau(i)>= 2
					b(i) = g(i+1,1)* prod(g(i,2:tau(i)));
				else
					b(i) = g(i+1,1);
				end
			end
			obj.b = b;
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
			obj.P = P;
			
			%%% h - Sweep three %%%
			h = zeros(1,N-1);
			h(1) = P(1,tau(1),tau(1))*g(2,1);
			for i = 2:N-1 %h defined for node N
				h(i) = (1-P(i,tau(i),1))*b(i)*h(i-1)+P(i,tau(i),tau(i))*g(i+1,1);
			end
			 obj.h = h;
			 
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
			 obj.phi = phi;
			 
			 %%% calculate a_k and c_k
			 a = zeros(1,N);
			 c = zeros(1,N);
			 for k = 2:N
				 a(k) = X(k,1)/r_vec(k) +gamma(k)/q_vec(k) * (1-X(k,1)/rho(k));
				 c(k) = -( X(k,1)/r_vec(k)- gamma(k)/q_vec(k)* X(k,1)/rho(k))  * (1-h(k-1)) ...
					 +gamma(k)/q_vec(k)*h(k-1);
			 end
			 obj.a = a;
			 obj.c = c;
		end
		
		%Corresponds to algorithm 2
        function [u, v] = calculate_input(obj,state,d)
            
			%%% calculate D %%%
			n = size(d,2);
			if n< 1+sum(obj.tau)+obj.H %Pad with zero if dist is too small
				d = [d, zeros(obj.N,1+sum(obj.tau)+obj.H-n)];
			end
			%D(i,j) = D_i[t+\sigma_i+j]
			D = zeros(obj.N,max([obj.tau,obj.H]));
			sigma = cumsum([0, obj.tau, obj.H]);
			for i = 1:obj.N
				for j = sigma(i):(sigma(i+1)-1)
					for k = 1:i
						D(i,j-sigma(i)+1) =  D(i,j-sigma(i)+1)+d(k,j+1-sigma(k));
					end
				end
			end
            
            
            %%% Calculate delta   %%%
            delta = zeros(1,obj.N);
            tau_extended = [obj.tau, obj.H];
            for i = 1:obj.N %loop over nodes
                Phi = state(i)*obj.phi(i);
                for Delta = 0:(tau_extended(i)-1)
                    if i ~=obj.N
                        Phi = Phi + obj.phi(i,Delta+1) * (D(i,Delta+1)+state(obj.N+sigma(i)+1+Delta));
					else %i == N, no u_N
                        Phi = Phi + obj.phi(i,Delta+1) * (D(i,Delta+1));
                    end
                end
                if i == 1 %First node
                    delta(i) = Phi;
                else %All but first node
                    delta(i) = Phi + (1-obj.P(i,tau_extended(i),1))*delta(i-1);
                end
            end
            %%% Calculate mu %%%
            mu = zeros(1,obj.N);
            for i = obj.N:-1:1 %loop over nodes
                pi = state(i);
                g_agg = 1; %used to keep track of the product of g_i[j]
                for Delta = 0:(tau_extended(i)-1)
                    if Delta ~= 0
                        g_agg = g_agg*obj.g(i,Delta+1);
                    end
                    if i ~=obj.N
                        pi = pi + g_agg * (D(i,Delta+1)+state(obj.N+sigma(i)+1+Delta));
                    else %For i = N there is no internal flows u_N
                        pi = pi + g_agg * (D(i,Delta+1));
                    end

                end
                if i == obj.N
                    mu(i) = pi;
                else
                    mu(i) = pi + obj.b(i)*mu(i+1);
                end
            end
            
            
            %calculate u
            u = zeros(obj.N-1,1);
            for i = 2:obj.N-1
               u(i-1) = (1-obj.gamma(i)/obj.q(i)) * (state(i)+state(obj.N+sigma(i)+1)+D(i,1))...
                   -obj.a(i)*delta(i-1) + obj.c(i)*mu(i)+d(i,1)-D(i,1);
            end
            u(obj.N-1) = (1-obj.gamma(obj.N)/obj.q(obj.N)) * (state(obj.N)+D(obj.N,1))...
                   -obj.a(obj.N)*delta(obj.N-1) + obj.c(obj.N)*mu(obj.N)+d(obj.N,1)-D(obj.N,1);
            
            %calculate v
            v = zeros(obj.N,1);
            v(1) = -obj.X(1,1)/obj.r(1)*mu(1);
            for i = 2:obj.N
               v(i) = -obj.X(i,1)/obj.r(i)*(delta(i-1)+(1-obj.h(i-1))*mu(i));
            end
        end

    end
end