clear all; clc; close all; 
format longE

% Simu_symlation Length
t_s = 0.02;
t_end = 5;
N = t_end / t_s;
time = [0:t_s:t_end];

% Parameter Definition
m_cart = 0.6;
m_pend = 0.3;
l_pend = 0.35;
g_grav = 9.81;

% Dimension
% n: dimension of state vector x
n = 4;
% m: dimension of observation vector y
m = 2;

% Initialize motion function
syms x1_sym x2_sym x3_sym x4_sym u_sym
f = [x2_sym;
    (l_pend*(x4_sym^2)*sin(x3_sym) - g_grav*cos(x3_sym)*sin(x3_sym) + u_sym/m_pend) / (m_cart/m_pend + (sin(x3_sym))^2);
    x4_sym;
    ((m_cart+m_pend)/(m_pend*l_pend)*g_grav*sin(x3_sym) - (x4_sym^2)*cos(x3_sym)*sin(x3_sym) - cos(x3_sym)*u_sym/(m_pend*l_pend)) / (m_cart/m_pend + (sin(x3_sym))^2)];
f = [x1_sym;x2_sym;x3_sym;x4_sym] + t_s * f;
dfdx = jacobian(f,[x1_sym;x2_sym;x3_sym;x4_sym]);
dfdu = jacobian(f,u_sym);

% Initialize observation function
g = [x1_sym;x3_sym];
dgdx = jacobian(g,[x1_sym;x2_sym;x3_sym;x4_sym]);

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
u_start = 0;
u_end = 5;
u_nominal = zeros([1 N]);
x_nominal = zeros([n 1 N]);
for iter = 1:N
    u_nominal(:,iter) = u_start + (u_end-u_start)/(N-1)*(iter-1);
end
% Simulate System for the nominals
x_o = [-1;0;0.1;0];
for k = 1:N
    % Compute linearized c2d-ed state matrices, we use this to simulate the
    % system because we do not know the discrete time nonlinear equation
    % for the inverse pendulum problem
    if k == 1
        x_tmp = x_o;
    else
        x_tmp = x_nominal(:,:,k-1);
    end
    [A_actual,b_actual,~] = c2d_invpend(x1_sym,x2_sym,x3_sym,x4_sym,u_sym,dfdx,dfdu,dgdx,x_tmp,u_nominal(:,k),t_s);
    A_nominal(:,:,k) = A_actual;
    b_nominal(:,k) = b_actual;
    % Time step
    x_nominal(:,:,k) = double(subs(f,{x1_sym,x2_sym,x3_sym,x4_sym,u_sym},{x_tmp(1),x_tmp(2),x_tmp(3),x_tmp(4),u_nominal(:,k)}));
end
clear k x_tmp A_actual b_actual
% plot(time(2:N+1),reshape(x_nominal(1,1,:),[N 1]))


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

x_o = [-1;0;0.1;0];
x_hat_o = [-0.9;0.5;0.21;-0.1];
sigma_o = 100 * eye(n);

% Simulate System without control
is_controlled = 0;
[x_uc,y_uc,x_hat_p_uc,u_uc] = simulate_invpend(N,n,m,t_s,x1_sym,x2_sym,x3_sym,x4_sym,u_sym,f,dfdx,dfdu,dgdx,xi,Q,nn,R,W,rho,lambda,x_o,x_hat_o,sigma_o,A_nominal,b_nominal,x_nominal,u_nominal,is_controlled);

% Simulate System with control
is_controlled = 1;
[x_c,y_c,x_hat_p_c,u_c] = simulate_invpend(N,n,m,t_s,x1_sym,x2_sym,x3_sym,x4_sym,u_sym,f,dfdx,dfdu,dgdx,xi,Q,nn,R,W,rho,lambda,x_o,x_hat_o,sigma_o,A_nominal,b_nominal,x_nominal,u_nominal,is_controlled);



% Process Result
x1_uc_plot = [x_o(1,1);reshape(x_uc(1,1,:),[N 1])];
x2_uc_plot = [x_o(2,1);reshape(x_uc(2,1,:),[N 1])];
x3_uc_plot = [x_o(3,1);reshape(x_uc(3,1,:),[N 1])];
x4_uc_plot = [x_o(4,1);reshape(x_uc(4,1,:),[N 1])];
y1_uc_plot = [x_o(1,1);reshape(y_uc(1,1,:),[N 1])];
y2_uc_plot = [x_o(2,1);reshape(y_uc(2,1,:),[N 1])];
x1_hat_uc_plot = [x_hat_o(1,1);reshape(x_hat_p_uc(1,1,:),[N 1])];
x2_hat_uc_plot = [x_hat_o(2,1);reshape(x_hat_p_uc(2,1,:),[N 1])];
x3_hat_uc_plot = [x_hat_o(3,1);reshape(x_hat_p_uc(3,1,:),[N 1])];
x4_hat_uc_plot = [x_hat_o(4,1);reshape(x_hat_p_uc(4,1,:),[N 1])];

x1_c_plot = [x_o(1,1);reshape(x_c(1,1,:),[N 1])];
x2_c_plot = [x_o(2,1);reshape(x_c(2,1,:),[N 1])];
x3_c_plot = [x_o(3,1);reshape(x_c(3,1,:),[N 1])];
x4_c_plot = [x_o(4,1);reshape(x_c(4,1,:),[N 1])];
y1_c_plot = [x_o(1,1);reshape(y_c(1,1,:),[N 1])];
y2_c_plot = [x_o(2,1);reshape(y_c(2,1,:),[N 1])];
x1_hat_c_plot = [x_hat_o(1,1);reshape(x_hat_p_c(1,1,:),[N 1])];
x2_hat_c_plot = [x_hat_o(2,1);reshape(x_hat_p_c(2,1,:),[N 1])];
x3_hat_c_plot = [x_hat_o(3,1);reshape(x_hat_p_c(3,1,:),[N 1])];
x4_hat_c_plot = [x_hat_o(4,1);reshape(x_hat_p_c(4,1,:),[N 1])];


figure(1)
plot(time,x1_uc_plot)
hold on
plot(time,y1_uc_plot)
hold on
plot(time,x1_hat_uc_plot)
legend('actual','observed','estimated')

figure(2)
plot(time,x3_uc_plot)
hold on
plot(time,y2_uc_plot)
hold on
plot(time,x3_hat_uc_plot)
legend('actual','observed','estimated')

figure(3)
rho1_plot = reshape(rho(1,1,:),[N+1 1]);
rho3_plot = reshape(rho(3,1,:),[N+1 1]);
plot(time,x1_c_plot)
hold on
plot(time,x1_uc_plot)
hold on
plot(time,rho1_plot)
legend('controlled','without control','target')

figure(4)
plot(time,x3_c_plot)
hold on
plot(time,x3_uc_plot)
hold on
plot(time,rho3_plot)
legend('controlled','without control','target')

figure(5)
plot(time(1:N),u_c)
legend('control sequence')

