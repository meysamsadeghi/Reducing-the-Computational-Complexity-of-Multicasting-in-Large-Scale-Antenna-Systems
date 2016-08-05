function [c_Huristic] = HuristicPrecoder(K,G,eta,tau,channel_Eqe)
c_Huristic = zeros(tau,G);
%
Increasing_order=zeros(G,K);
Chan_eq_rearraneged=zeros(tau,K,G);
eta_rearranged = zeros(K,G);
% Finding the weakest user
for g=1:G
    Chnl_Qulty = zeros(K,1); 
    for k=1:K
        Chnl_Qulty(k,1)= eta(k,1) / (norm(channel_Eqe(:,k,g))^2) ; % the simple ordering measure introduced in line 2 of Alg. 3
    end
    [~,Increasing_order(g,:)]=sort(Chnl_Qulty,'descend'); % line 2 of Alg. 3
% Channel Rearrangement
    for k=1:K
        Chan_eq_rearraneged(:,k,g) = channel_Eqe(:,Increasing_order(g,k),g);
        eta_rearranged(k,g)=eta(Increasing_order(g,k),g);
    end
    [Q,~]=gramSmithGen(Chan_eq_rearraneged(:,:,g)); % Perform Gram-Smith procedure as line 3 of Alg. 3
    % add the weakest user
    weakest_user = Increasing_order(g,1); % line 4 of Alg. 3
    c_Huristic(:,g) = ( (sqrt(eta(weakest_user,g))) / (norm(Chan_eq_rearraneged(:,1,g)'*Q(:,1))) ) * Q(:,1); % line 4 of Alg. 3
    for k = 2:K % line 5 of Alg. 3
        if abs(Chan_eq_rearraneged(:,k,g)'*c_Huristic(:,g))^2 < eta_rearranged(k,g) % line 6 of Alg. 3
            % update
            InsideReal = c_Huristic(:,g)'*Chan_eq_rearraneged(:,k,g)*Chan_eq_rearraneged(:,k,g)'*Q(:,k);
            Ang_alpha = - angle(InsideReal);
            A_quadratic = abs(Q(:,k)'*Chan_eq_rearraneged(:,k,g))^2;
            B_quadratic = 2*real(exp(1i*Ang_alpha)*InsideReal);
            C_quadratic = abs(c_Huristic(:,g)'*Chan_eq_rearraneged(:,k,g))^2 - eta_rearranged(k,g); %eta(k,g);
            Abs_alpha = (-B_quadratic + sqrt(B_quadratic^2 - 4*A_quadratic*C_quadratic)) / (2*A_quadratic);
            alpha = Abs_alpha * exp(1i*Ang_alpha);
            c_Huristic(:,g) = c_Huristic(:,g) + (alpha * Q(:,k)); % line 8 of Alg. 3
        end
    end
end









