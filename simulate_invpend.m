function [x,y,x_hat_p,u] = simulate_invpend(N,n,m,t_s,x1_sym,x2_sym,x3_sym,x4_sym,u_sym,f,dfdx,dfdu,dgdx,xi,Q,nn,R,W,rho,lambda,x_o,x_hat_o,sigma_o,A_nominal,b_nominal,x_nominal,u_nominal,is_controlled)

% Initialize state vector x
x = zeros([n 1 N]);

% Initialize control u
% Note: in our ase control is scalar
u = zeros([1 N]);

% Initialize observation y
y = zeros([m 1 N]);

% Initialize predicted/updated state x_hat
x_hat_a = zeros([n 1 N]);
x_hat_p = zeros([n 1 N]);

% Initialize predicted/updated covariance sigma
sigma_a = zeros([n n N]);
sigma_p = zeros([n n N]);

% Initilaize Kalman Gain kg
kg = zeros([n m N]);

% Compute the solution to Riccati Equation
% Solve for linearized state matrices for the target trajectory
% A = zeros([n n N]);
% b = zeros([n 1 N]);
% [A_const,b_const] = c2d_invpend(x1_sym,x2_sym,x3_sym,x4_sym,u_sym,dfdx,dfdu,dgdx,x_hat_o,0,t_s);
% for j = 1:N
%     A(:,:,j) = A_const;
%     b(:,:,j) = b_const;
% end
% [K_tilda,p_tilda,mu] = solve_riccati(N,n,W,rho,lambda,A,b);

[K_tilda,p_tilda,mu] = solve_riccati(N,n,W,rho,lambda,A_nominal,b_nominal);
p_tilda = zeros([n 1 N+1]);


% Simulate System
% x_nominal = zeros([n 1 N]);
for k = 1:N
    % Calculate for control
    if is_controlled == 0
        u(:,k) = 0;
    else
        % Update the nominal trajectory based on current estiamte
        if k == 1
            x_hat_km1 = x_hat_o;  % Choose the predicted state to feed into algorithm
        else
            x_hat_km1 = reshape(x_hat_p(:,:,k-1),[n 1]);
        end
%         for j = k:N
%             if j == k
%                 x_tmp = x_hat_km1;
%             else
%                 x_tmp = x_nominal(:,:,j-1);
%             end
%             [A_new,b_new,~] = c2d_invpend(x1_sym,x2_sym,x3_sym,x4_sym,u_sym,dfdx,dfdu,dgdx,x_tmp,u_nominal(:,j),t_s);
%             A_nominal_new(:,:,j) = A_new;
%             b_nominal_new(:,j) = b_new;
%             x_nominal(:,:,j) = A_new * x_tmp + b_new * u_nominal(:,j);
%         end
%         [K_tilda,p_tilda,mu] = solve_riccati(N,n,W,rho,lambda,A_nominal,b_nominal);

% 2
        delta_x_hat_km1 = x_hat_km1 - rho(:,:,k);
        delta_u = calc_control(A_nominal,b_nominal,K_tilda,p_tilda,mu,delta_x_hat_km1,k);
        u(:,k) = u_nominal(:,k) + delta_u;
    end

    % Compute linearized c2d-ed state matrices, we use this to simulate the
    % system because we do not know the discrete time nonlinear equation
    % for the inverse pendulum problem
%     if k == 1
%         x_tmp = x_o;
%     else
%         x_tmp = x(:,:,k-1);
%     end
%     [A_actual,b_actual,C_actual] = c2d_invpend(x1_sym,x2_sym,x3_sym,x4_sym,u_sym,dfdx,dfdu,dgdx,x_tmp,u(:,k),t_s);

    % Kalman Filter Prediction
    if k == 1
        x_tmp = x_hat_o;
        sigma_tmp = sigma_o;
    else
        x_tmp = x_hat_p(:,:,k-1);
        sigma_tmp = sigma_p(:,:,k-1);
    end
    [A_ekf,b_ekf,C_ekf] = c2d_invpend(x1_sym,x2_sym,x3_sym,x4_sym,u_sym,dfdx,dfdu,dgdx,x_tmp,u(:,k),t_s);
    % x_hat_a(:,:,k) = A_actual * x_tmp + b_actual * u(:,k);
    x_hat_a(:,:,k) = double(subs(f,{x1_sym,x2_sym,x3_sym,x4_sym,u_sym},{x_tmp(1),x_tmp(2),x_tmp(3),x_tmp(4),u(:,k)}));
    sigma_a(:,:,k) = A_ekf * sigma_tmp * A_ekf' + Q(:,:,k);
    kg(:,:,k) = (sigma_a(:,:,k)*C_ekf') * inv(C_ekf*sigma_a(:,:,k)*C_ekf'+R(:,:,k));
    % Time step
    if k == 1
        x_tmp = x_o;
        % x_nom_tmp = x_o;
    else
        x_tmp = x(:,:,k-1);
        % x_nom_tmp = x_nominal(:,:,k-1);
    end
    % x(:,:,k) = A_actual * x_tmp + b_actual * u(:,k) + xi(:,:,k);
    x(:,:,k) = double(subs(f,{x1_sym,x2_sym,x3_sym,x4_sym,u_sym},{x_tmp(1),x_tmp(2),x_tmp(3),x_tmp(4),u(:,k)})) + xi(:,:,k);
    % Observer
    % y(:,:,k) = C_actual * x(:,:,k) + nn(:,:,k);
    y(:,:,k) = double(dgdx) * x(:,:,k) + nn(:,:,k);
    % Kalman Filter Update
    x_hat_p(:,:,k) = x_hat_a(:,:,k) + kg(:,:,k) * (y(:,:,k) - double(dgdx) * x_hat_a(:,:,k));
    sigma_p(:,:,k) = (eye(n) - kg(:,:,k) * C_ekf) * sigma_a(:,:,k);
end
clear x_tmp x_nom_tmp sigma_tmp x_hat_k_km1 sigma_kp1_k


end