function [SINR_Huri] = SINR_Huri_Gen(W_Huri,G,K,channel,sigma_sqrd)
%%   SDR
SINR_Huri=zeros(G,K);
for g=1:G
    for k=1:K
        %% SDR precoding Test
        % signal power
        sig_pow = (abs(channel(:,k,g)' * W_Huri(:,g)))^2;
        % intra cell interference
        anbar = 0;
        for g_intra = 1:G
            if g_intra ~= g
                anbar = anbar + (abs(channel(:,k,g)' * W_Huri(:,g_intra)))^2;
            end
        end
        % SINR 
        SINR_Huri(g,k) = sig_pow / (anbar + sigma_sqrd(k,g));     
    end
end