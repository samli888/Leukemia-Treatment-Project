function [x,y,x_hat_p,u] = simulate_system(N,n,m,A,b,C,xi,Q,nn,R,W,rho,lambda,x_o,is_controlled)

% Initialize state vector x
x = zeros([n 1 N]);
%x_o = [100;0;-10];

% Initialize control u
% Note: in our ase control is scalar
u = zeros([1 N]);

% Initialize observation y
y = zeros([m 1 N]);

% Initialize predicted/updated state x_hat
x_hat_a = zeros([n 1 N]);
x_hat_p = zeros([n 1 N]);
x_hat_o = x_o;

% Initialize predicted/updated covariance sigma
sigma_a = zeros([n n N]);
sigma_p = zeros([n n N]);
sigma_o = [1 0 0; 0 1 0; 0 0 1]; % small initial covariance because we assuming initial state known

% Initilaize Kalman Gain kg
kg = zeros([n m N]);

% Compute the solution to Riccati Equation
[K_tilda,p_tilda,mu] = solve_riccati(N,n,W,rho,lambda,A,b);

% Simulate System
for k = 1:N
    % Control Search Algorithm
    if k == 1
        x_hat_km1 = x_hat_o;  % Choose the predicted state to feed into algorithm
    else
        x_hat_km1 = reshape(x_hat_p(:,:,k-1),[n 1]);
    end

    if is_controlled == 0
        u(:,k) = 0;
    else
        u(:,k) = calc_control(A,b,K_tilda,p_tilda,mu,x_hat_km1,k);
    end

    % Kalman Filter Prediction
    if k == 1
        x_tmp = x_hat_o;
        sigma_tmp = sigma_o;
    else
        x_tmp = x_hat_p(:,:,k-1);
        sigma_tmp = sigma_p(:,:,k-1);
    end
    x_hat_a(:,:,k) = A(:,:,k) * x_tmp + b(:,k) * u(:,k);
    sigma_a(:,:,k) = A(:,:,k) * sigma_tmp * A(:,:,k)' + Q(:,:,k);
    kg(:,:,k) = (C(:,:,k)*sigma_a(:,:,k)*C(:,:,k)'+R(:,:,k)) \ (sigma_a(:,:,k)*C(:,:,k)');
    % Time step
    if k == 1
        x_tmp = x_o;
    else
        x_tmp = x(:,:,k-1);
    end
    x(:,:,k) = A(:,:,k) * x_tmp + b(:,k) * u(:,k) + xi(:,:,k);
    % Observer
    y(:,:,k) = C(:,:,k) * x(:,:,k) + nn(:,:,k);
    % Kalman Filter Update
    x_hat_p(:,:,k) = x_hat_a(:,:,k) + kg(:,:,k) * (y(:,:,k) - C(:,:,k) * x_hat_a(:,:,k));
    sigma_p(:,:,k) = (eye(m) - kg(:,:,k) * C(:,:,k)) * sigma_a(:,:,k);
end
clear x_tmp sigma_tmp x_hat_k_km1 sigma_kp1_k


end