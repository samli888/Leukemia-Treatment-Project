clear all; clc; close all; 
format longE

% Simulation Length
N = 10;
time = [0:1:N];

% Dimension
% n: dimension of state vector x
n = 3;
% m: dimension of observation vector y
m = 3;

% Initialize state matrix A
A = zeros([n n N]);
dt = 1;
A_const = [1 dt 0.5*dt^2; 0 1 dt; 0 0 1];
for iter = 1:N
    A(:,:,iter) = A_const;
end
clear A_cosnt iter

% Initialize control matrix B
% Note: in our case control is scalar so B becomes a vector
b = zeros([n N]);
b_const = [0;0;1];
for iter = 1:N
    b(:,iter) = b_const;
end
clear b_cosnt iter

% Initialize observe matrix C
C = zeros([m n N]);
C_const = eye(m);
for iter = 1:N
    C(:,:,iter) = C_const;
end
clear C_cosnt iter

% Initialize disturbance xi and covariance Q
xi = zeros([n 1 N]);
Q = zeros([n n N]);
sig = 1;
Q_const = sig * eye(n);
for iter = 1:N
    Q(:,:,iter) = Q_const;
    xi(:,:,iter) = Q(:,:,iter) * randn([n,1]);
end
clear sig Q_const

% Initialize measurement noise n and covariance R
nn = zeros([m 1 N]);
R = zeros([m m N]);
sig = 5;
R_const = sig * eye(m);
for iter = 1:N
    R(:,:,iter) = R_const;
    nn(:,:,iter) = R(:,:,iter) * randn([m,1]);
end
clear sig R_const


% Initialize Cost Parameters
W = zeros([n n N+1]);
rho = zeros([n 1 N+1]);
lambda = zeros([1 N]);
W_const = [1 0 0; 0 1 0; 0 0 1];
lambda_const = 0.001;
target_acceleration = -15;
start_point = [100;0;-10];
for iter = 1:N+1
    W(:,:,iter) = W_const;
    t = (iter - 1) * dt;
    rho_1 = start_point(1,1) + start_point(2,1)*t + 0.5*(t^2)*target_acceleration;
    rho_2 = start_point(2,1) + t*target_acceleration;
    rho_3 = target_acceleration;
    rho(:,:,iter) = [rho_1;rho_2;rho_3];
end
for iter = 1:N
    lambda(:,iter) = lambda_const;
end
clear W_const lambda_const dt


% Note: Iteration variables does not correspond to the notes due to matlab
% indexing start with 1 instead of 0

x_o = start_point;

% Simulate System with control
is_controlled = 1;
[x_c,y_c,x_hat_p_c,u_c] = simulate_system(N,n,m,A,b,C,xi,Q,nn,R,W,rho,lambda,x_o,is_controlled);

% Compute cost
J_with_control = compute_cost(N,x_c,x_o,u_c,W,rho,lambda)


% Simulate System without control
is_controlled = 0;
[x_uc,y_uc,x_hat_p_uc,u_uc] = simulate_system(N,n,m,A,b,C,xi,Q,nn,R,W,rho,lambda,x_o,is_controlled);


% Compute cost
J_without_control = compute_cost(N,x_uc,x_o,u_uc,W,rho,lambda)



% Process Result
x_c_plot = [x_o(1,1);reshape(x_c(1,1,:),[N 1])];
dx_c_plot = [x_o(2,1);reshape(x_c(2,1,:),[N 1])];
ddx_c_plot = [x_o(3,1);reshape(x_c(3,1,:),[N 1])];
y1_c_plot = [x_o(1,1);reshape(y_c(1,1,:),[N 1])];
y2_c_plot = [x_o(2,1);reshape(y_c(2,1,:),[N 1])];
y3_c_plot = [x_o(3,1);reshape(y_c(3,1,:),[N 1])];
x_hat_c_plot = [x_o(1,1);reshape(x_hat_p_c(1,1,:),[N 1])];
dx_hat_c_plot = [x_o(2,1);reshape(x_hat_p_c(2,1,:),[N 1])];
ddx_hat_c_plot = [x_o(3,1);reshape(x_hat_p_c(3,1,:),[N 1])];

x_uc_plot = [x_o(1,1);reshape(x_uc(1,1,:),[N 1])];
dx_uc_plot = [x_o(2,1);reshape(x_uc(2,1,:),[N 1])];
ddx_uc_plot = [x_o(3,1);reshape(x_uc(3,1,:),[N 1])];
y1_uc_plot = [x_o(1,1);reshape(y_uc(1,1,:),[N 1])];
y2_uc_plot = [x_o(2,1);reshape(y_uc(2,1,:),[N 1])];
y3_uc_plot = [x_o(3,1);reshape(y_uc(3,1,:),[N 1])];
x_hat_uc_plot = [x_o(1,1);reshape(x_hat_p_uc(1,1,:),[N 1])];
dx_hat_uc_plot = [x_o(2,1);reshape(x_hat_p_uc(2,1,:),[N 1])];
ddx_hat_uc_plot = [x_o(3,1);reshape(x_hat_p_uc(3,1,:),[N 1])];

figure(1)
plot(time,x_c_plot)
hold on
plot(time,y1_c_plot)
hold on
plot(time,x_hat_c_plot)
legend('actual','observed','estimated')

figure(2)
rho1_plot = reshape(rho(1,1,:),[N+1 1]);
rho2_plot = reshape(rho(2,1,:),[N+1 1]);
rho3_plot = reshape(rho(3,1,:),[N+1 1]);
plot(time,x_c_plot)
hold on
plot(time,x_uc_plot)
hold on
plot(time,rho1_plot)
legend('controlled','without control','target')



