% Simu_symlation Length
t_s = step_size_noisy;
t_end = 21;
N = num_t_noisy;
time = [1:t_s:t_end];

% Dimension
% n: dimension of state vector x
n = 8;
% m: dimension of observation vector y
m = 1;

% Initialize motion function
syms x1_sym x2_sym x3_sym x4_sym x5_sym x6_sym x7_sym x8_sym u_sym

f = [-theta(1)                0 0 0 0 0 0 0; ...
    theta(1) -theta(2)          0 0 0 0 0 0; ...
    0 theta(3)*theta(4) -theta(5) 0 0 0 0 0; ...
    0 0 0 -theta(7)                 0 0 0 0; ...
    0 0 0 theta(7) -theta(7)          0 0 0; ...
    0 0 0 0 theta(7) -theta(7)          0 0; ...
    0 0 0 0 0 theta(7) -theta(7)          0; ...
    0 0 0 0 0 0 theta(7) -theta(10)]*[x1_sym;x2_sym;x3_sym;x4_sym;x5_sym;x6_sym;x7_sym;x8_sym] + ...
    [0.22;0;0;0;0;0;0;0]*u_sym + ...
    [0;...
    0;...
    0;...
    theta(7)*((theta(6)/x8_sym)^theta(9))*x4_sym-theta(7)*theta(8)*((theta(6)/x8_sym)^theta(9))*x3_sym*x4_sym;...
    0;...
    0;...
    0;...
    0];

f = [x1_sym;x2_sym;x3_sym;x4_sym;x5_sym;x6_sym;x7_sym;x8_sym] + t_s * f;

dfdx = jacobian(f,[x1_sym;x2_sym;x3_sym;x4_sym;x5_sym;x6_sym;x7_sym;x8_sym]);
dfdu = jacobian(f,u_sym);

% Initialize observation function
g = [x8_sym];
dgdx = jacobian(g,[x1_sym;x2_sym;x3_sym;x4_sym;x5_sym;x6_sym;x7_sym;x8_sym]);

% Initialize disturbance xi and covariance Q
xi = zeros([n 1 N]);
Q = zeros([n n N]);
sig = 0.01;
Q_const = sig * eye(n);
for iter = 1:N
    Q(:,:,iter) = Q_const;
    xi(:,:,iter) = Q(:,:,iter) * randn([n,1]);
end
clear sig Q_const

% Initialize measurement noise n and covariance R
nn = zeros([m 1 N]);
R = zeros([m m N]);
sig = 0.1;
R_const = sig * eye(m);
for iter = 1:N
    R(:,:,iter) = R_const;
    nn(:,:,iter) = R(:,:,iter) * randn([m,1]);
end
clear sig R_const

% Initialize nominal trajectory and control to track
% nominal control is set as linear for test purposes
A_nominal = zeros([n n N]);
b_nominal = zeros([n N]);

% convert and reshape for Sam's code
u_nominal = transpose(u_ref_flattened(1:N,1));

x_nominal = zeros([n 1 N]);
% Simulate System for the nominals
x_o = x0;
for k = 1:N
    % Compute linearized c2d-ed state matrices, we use this to simulate the
    % system because we do not know the discrete time nonlinear equation
    % for the inverse pendulum problem
    if k == 1
        x_tmp = x_o;
    else
        x_tmp = x_nominal(:,:,k-1);
    end
    [A_actual,b_actual,~] = c2d_jost(x1_sym,x2_sym,x3_sym,x4_sym,x5_sym,x6_sym,x7_sym,x8_sym,u_sym,dfdx,dfdu,dgdx,x_tmp,u_nominal(:,k),t_s);
    A_nominal(:,:,k) = A_actual;
    b_nominal(:,k) = b_actual;
    % Time step
    x_nominal(:,:,k) = double(subs(f,{x1_sym,x2_sym,x3_sym,x4_sym,x5_sym,x6_sym,x7_sym,x8_sym,u_sym},{x_tmp(1),x_tmp(2),x_tmp(3),x_tmp(4),x_tmp(5),x_tmp(6),x_tmp(7),x_tmp(8),u_nominal(:,k)}));
end
clear k x_tmp A_actual b_actual


% Initialize Cost Parameters
W = zeros([n n N+1]);
rho = zeros([n 1 N+1]);
lambda = zeros([1 N]);
W_const = 100 * eye(n);
lambda_const = 0.1;
for iter = 1:N+1
    W(:,:,iter) = W_const;
    if iter == 1
        rho(:,:,iter) = x_o;
    else
        rho(:,:,iter) = x_nominal(:,:,iter-1);
    end
end
for iter = 1:N
    lambda(:,iter) = lambda_const;
end
clear W_const lambda_const

% Note: Iteration variables does not correspond to the notes due to matlab
% indexing start with 1 instead of 0

x_o = x0;
x_hat_o = x0+[-0.3;0.2;0.05;0.01;-0.5;0.01;0;-0.32]; % CHANGE this?
sigma_o = 0.2 * eye(n); % CHANGE THIS?

% Simulate System with control
is_controlled = 1;
[x_c,y_c,x_hat_p_c,u_c] = simulate_jost(N,n,m,t_s,x1_sym,x2_sym,x3_sym,x4_sym,x5_sym,x6_sym,x7_sym,x8_sym,u_sym,f,dfdx,dfdu,dgdx,xi,Q,nn,R,W,rho,lambda,x_o,x_hat_o,sigma_o,A_nominal,b_nominal,x_nominal,u_nominal,is_controlled);


% Process Result
x1_c_plot = [x_o(1,1);reshape(x_c(1,1,:),[N 1])];
x2_c_plot = [x_o(2,1);reshape(x_c(2,1,:),[N 1])];
x3_c_plot = [x_o(3,1);reshape(x_c(3,1,:),[N 1])];
x4_c_plot = [x_o(4,1);reshape(x_c(4,1,:),[N 1])];
x5_c_plot = [x_o(5,1);reshape(x_c(5,1,:),[N 1])];
x6_c_plot = [x_o(6,1);reshape(x_c(6,1,:),[N 1])];
x7_c_plot = [x_o(7,1);reshape(x_c(7,1,:),[N 1])];
x8_c_plot = [x_o(8,1);reshape(x_c(8,1,:),[N 1])];
y1_c_plot = [x_o(1,1);reshape(y_c(1,1,:),[N 1])];
x1_hat_c_plot = [x_hat_o(1,1);reshape(x_hat_p_c(1,1,:),[N 1])];
x2_hat_c_plot = [x_hat_o(2,1);reshape(x_hat_p_c(2,1,:),[N 1])];
x3_hat_c_plot = [x_hat_o(3,1);reshape(x_hat_p_c(3,1,:),[N 1])];
x4_hat_c_plot = [x_hat_o(4,1);reshape(x_hat_p_c(4,1,:),[N 1])];
x5_hat_c_plot = [x_hat_o(5,1);reshape(x_hat_p_c(5,1,:),[N 1])];
x6_hat_c_plot = [x_hat_o(6,1);reshape(x_hat_p_c(6,1,:),[N 1])];
x7_hat_c_plot = [x_hat_o(7,1);reshape(x_hat_p_c(7,1,:),[N 1])];
x8_hat_c_plot = [x_hat_o(8,1);reshape(x_hat_p_c(8,1,:),[N 1])];
