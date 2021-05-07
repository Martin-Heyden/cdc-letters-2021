# cdc-letters-2021

The files calculate_inputs.m and generate_controller.m can be used to generate the optimal 
controller and calculate the optimal input by calculating the optimal K using the results
in the paper.

structured_controller.m calculated the optimal inputs using algorithm 2 in the paper.

dist_ex.m and feedforward_comp.m are used to generate the figures in the paper.

test_generate_controller.m verifies that the results are numerically correct by comparing to DARE and convex optimization (basically MPC).
test_structured_controller.m verifes that the implementation in test_structured_controller is correct by comparing it to generate_controller.

generate_graph.m is used to generate a state space representation for the dynamics.
