function [F,channel_Eqe] = SVD_preliminaries(channel,M,K,G,sigma_sqrd,tau)

%% 1. F generator
F = zeros(M,tau,G); % the last G represent each group
for g = 1:G
    countr = 1;
    SAI = zeros(K*(G-1),M);
    for g_in=1:G
        if g_in ~= g
            SAI( ((countr-1)*K)+1:countr*K ,:) = channel(:,:,g_in)';
            countr = countr +1;
        end
    end
    %[~,~,V] = svd(SAI);
    [V_test, ~] = qr(SAI');
    F(:,:,g) = V_test(:,(K*G)-K+1:end);%V(:,(K*G)-K+1:end);
end

%% 2. g equivalent <=> gbar
channel_Eqe = zeros(tau,K,G);
for g=1:G
    for k=1:K
        channel_Eqe(:,k,g) = (1/sqrt(sigma_sqrd(k,g))) * (F(:,:,g)' * channel(:,k,g));
    end
end
