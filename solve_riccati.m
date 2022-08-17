function [K_tilda,p_tilda,mu] = solve_riccati(N,n,W,rho,lambda,A,b)

% Solving Riccati Equation
K_tilda = zeros([n n N+1]);
p_tilda = zeros([n 1 N+1]);
mu = zeros([1 N]);
for j = N+1:-1:1
    if j == N+1
        K_tilda(:,:,j) = W(:,:,N+1);
        p_tilda(:,:,j) = -1*W(:,:,N+1)*rho(:,:,N+1);
    else
        mu(:,j) = 1/(lambda(:,j) + b(:,j)'*K_tilda(:,:,j+1)*b(:,j));
        K_tilda(:,:,j) = A(:,:,j)' * (eye(n)-mu(:,j)*K_tilda(:,:,j+1)*b(:,j)*b(:,j)') * K_tilda(:,:,j+1) * A(:,:,j) + W(:,:,j);
        p_tilda(:,:,j) = A(:,:,j)' * (eye(n)-mu(:,j)*K_tilda(:,:,j+1)*b(:,j)*b(:,j)') * p_tilda(:,:,j+1) - W(:,:,j)*rho(:,:,j);
    end
end

end