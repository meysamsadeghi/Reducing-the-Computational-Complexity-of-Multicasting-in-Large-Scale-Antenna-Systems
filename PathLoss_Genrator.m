function [PathLoss] = PathLoss_Genrator(K,G,Terminal_pos)

PathLoss = zeros(K,G); % path-loss array

for g = 1:G %
    for k = 1:K % user
        Distance = sqrt( Terminal_pos(k,g,1)^2 + Terminal_pos(k,g,2)^2 );
        pathloss_in_db = 128.1 + 37.6 * log10(Distance/1000);
        PathLoss(k,g) = 10^(-pathloss_in_db /10);
    end
end


