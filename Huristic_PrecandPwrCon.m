function [W_Huri,PowHuristic] = Huristic_PrecandPwrCon(F,c_Huristic,M,G)
%% Precoder
W_Huri =zeros(M,G);
for g=1:G
    W_Huri(:,g)  = F(:,:,g) * c_Huristic(:,g);
end
%% Power Consumption for comparison puprposes
PowHuristic = 0;
for g=1:G
    PowHuristic = PowHuristic + (W_Huri(:,g)' * W_Huri(:,g));
end