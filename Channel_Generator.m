function [channel,BigChannel] = Channel_Generator(M,K,G,PathLoss)
channel = zeros(M,K,G); 
BigChannel = zeros(M,K*G);

%%
cntr=1; % to make big G for unicast
for g=1:G
    for k=1:K
        channel(:,k,g) = sqrt(PathLoss(k,g))  * sqrt(1/2) * (randn(M,1) + 1i*randn(M,1)) ; %* randn(M,1);%
    end
    BigChannel(:,((cntr-1)*K+1):(cntr*K)) = channel(:,:,g);
    cntr=cntr+1;
end
