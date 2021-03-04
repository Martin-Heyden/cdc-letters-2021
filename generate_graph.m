 function [ A,B,Q,R ] = generate_graph(edges,prod, Qvec, Rvec,edge_delay )
%Generate graph defined by edges. First node states, then delay
%states. Order correspond to order in prod, then edges.
%All internal transportation is assumed to be delayed!
number_of_prod = length(prod);
number_of_int = length(edges);
number_of_internal_states = sum(edge_delay);


number_of_inputs = number_of_int+number_of_prod;
number_of_nodes = size(edges,1)+1;
number_of_states = number_of_nodes+ number_of_internal_states;
A = zeros(number_of_states,number_of_states);
B = zeros(number_of_states,number_of_inputs);
Q = zeros(number_of_states,number_of_states);
R = zeros(number_of_inputs,number_of_inputs);

A(1:number_of_nodes,1:number_of_nodes) = eye(number_of_nodes,number_of_nodes);
Q(1:number_of_nodes,1:number_of_nodes) = diag(Qvec);
R(1:number_of_prod,1:number_of_prod) = diag(Rvec);

%%% States and inputs related to Production/Consumption %%%
for i = 1:number_of_prod
   dest = prod(i);
   B(dest,i) = 1;
end

%%% States and inputs related to internal flows %%%
current_int_delay = number_of_nodes+1;
for i = number_of_prod+1:number_of_prod+number_of_int
    edge = i - number_of_prod;
    source = edges(edge,1);
    dest = edges(edge,2);
    B(source,i) = -1; %Outflow from source
    if edge_delay(edge)<1
        error("Invalid edge delay")
    else
        A(dest,current_int_delay) = 1; %Delivery to dest
        if edge_delay(edge) == 1
            B(current_int_delay,i) = 1;
        else
            for j = current_int_delay+1:current_int_delay+(edge_delay(edge)-1)
                A(j-1,j) = 1;
            end
            B(j,i) = 1;
        end
           
        current_int_delay = current_int_delay + edge_delay(edge);

    end
end



end