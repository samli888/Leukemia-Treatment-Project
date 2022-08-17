function [J] = compute_cost(N,x,x_o,u,W,rho,lambda)

% Compute cost
J = (x(:,:,N) - rho(:,:,N+1))' * W(:,:,N+1) * (x(:,:,N) - rho(:,:,N+1));
for k = 1:N
    J = J + lambda(:,k)*(u(:,k)^2);
    if k == 1
        x_tmp = x_o;
    else
        x_tmp = x(:,:,k-1);
    end
    J = J + (x_tmp - rho(:,:,k))' * W(:,:,k) * (x_tmp - rho(:,:,k));
end
J = J*0.5;

end