function [u] = calc_control(A,b,K_tilda,p_tilda,mu,x_hat_km1,k)

% Calculate for control for time k
u = -1*mu(:,k)*b(:,k)'*(K_tilda(:,:,k+1)*A(:,:,k)*x_hat_km1 + p_tilda(:,:,k+1));

end