function [A,b,C] = c2d_jost(x1_sym,x2_sym,x3_sym,x4_sym,x5_sym,x6_sym,x7_sym,x8_sym,u_sym,dfdx,dfdu,dgdx,x,u_double,t_s)

% This function finds the state matrices in discrete time for the inverted 
% pendulum problem only with n = 4 and m = 2
A = double(subs(dfdx,{x1_sym,x2_sym,x3_sym,x4_sym,x5_sym,x6_sym,x7_sym,x8_sym,u_sym},{x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),u_double}));
b = double(subs(dfdu,{x1_sym,x2_sym,x3_sym,x4_sym,x5_sym,x6_sym,x7_sym,x8_sym,u_sym},{x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),u_double}));
C = double(dgdx);

end