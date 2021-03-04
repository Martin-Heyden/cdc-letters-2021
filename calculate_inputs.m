function [v,u] = calculate_inputs(Kx, Kd, gamma_Nq, x, dist, H,tau)
%   Calculated shifted disturbances from dist, and uses that to find the
%   optimal v and u.
%   dist is assumed to be of size N x q, where N is the number of nodes
%   


    N = length(tau)+1;
    q = size(dist,2);
    if q< 1+sum(tau)+H %Pad with zero if dist is too small
        dist = [dist, zeros(N,1+sum(tau)+H-q)];
    end
    D = zeros(N,sum(tau)+H); %Stores partial shifted dist vectors
    Dvec = zeros(sum(tau),1); %Shifted disturbances needed for controller
                                %for all nodes except node N
    DN_vec = zeros(H,1);%Shifted disturbances for node N

    sigma = [0 cumsum([tau H])];
    n = N+1; %First N states corresponds to node levels
            %The next corresponds to transit states
            % (The disturbances have the same effect as transit states)
    m = 1;
    %calculate shifted disturbances
    for i = 1:N
        for j = sigma(i):(sigma(i+1)-1)
           for k = 1:i
              D(i,j+1) =  D(i,j+1)+dist(k,j+1-sigma(k));
           end
           if i~=N
                Dvec(n) = D(i,j+1);
                n = n+1;
           else
               DN_vec(m) = D(i,j+1);
               m = m+1;
           end
        end
    end
    v = zeros(N,1);
    u = zeros(N-1,1);
    for i =1:N
         v(i) = Kx(i,:)*(x+Dvec) + Kd(i,:)*DN_vec;
    end
    for i = 2:N-1
         u(i-1) = Kx(N+i-1,:)*(x+Dvec)+dist(i,1)-D(i,sigma(i)+1)+ Kd(N+i-1,:)*DN_vec;
    end
      % For u_{N-1} D_N is not part of Kx, and thus needs special treatment
    u(end) =  Kx(2*N-1,:)*(x+Dvec)+dist(N,1)-gamma_Nq*D(N,sigma(N)+1)+ Kd(2*N-1,:)*DN_vec;


end

