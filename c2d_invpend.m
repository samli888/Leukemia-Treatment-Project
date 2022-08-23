function [A,b,C] = c2d_invpend(x1_sym,x2_sym,x3_sym,x4_sym,u_sym,dfdx,dfdu,dgdx,x,u_double,t_s)

% This function finds the state matrices in discrete time for the inverted 
% pendulum problem only with n = 4 and m = 2
A = double(subs(dfdx,{x1_sym,x2_sym,x3_sym,x4_sym,u_sym},{x(1),x(2),x(3),x(4),u_double}));
b = double(subs(dfdu,{x1_sym,x2_sym,x3_sym,x4_sym,u_sym},{x(1),x(2),x(3),x(4),u_double}));
C = double(dgdx);
% d = [0;0];
% sys_ss = ss(A,b,C,d);
% sys_d = c2d(sys_ss,t_s,'zoh');
% [A,b,C] = dssdata(sys_d);

end