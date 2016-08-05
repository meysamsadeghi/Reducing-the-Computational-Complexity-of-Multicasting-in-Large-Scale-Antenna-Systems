%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Multicasting in the context of Massive MIMO
% Version: 1
% The Aim of Code: to find an approximate solution of the QoS problem with
% the proposed Heuristic Algorith, i.e., BDZF + Algorithm 3
% System Model: as considered in the manuscript
% This code enables you to check the time and power consumption of QoS
% problem while using BDZF + Alg. 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear all;
close all;
%% General Initialization
MnoteCarlo_LSF = 100;                                         % Number of MonteCarlos when the large scale fading is changing  
MnoteCarlo_SSF = 1;                                           % Number of MonteCarlos when just the small scale fading is changing
r = 900;                                                      % Cell radius
Marray = [40:10:90] ;                                         % Number o antenna at BS, changing from 40 to 90
G = 3;                                                        % number of groups
K = 10;                                                       % number of users per group
L = 1;                                                        % number of cells
eta =  255 *  ones(K,G);                                      % The Prescribed SINR, can be anything, e.g., randi(512,K,G)
sigma_sqrd = 20 * 10^(-14.4) * ones(K,G);                     % Noise of a 20MHz BW channel !
Pwr_Huri=zeros(MnoteCarlo_LSF,MnoteCarlo_SSF,length(Marray)); % Power Consumption
timeHuri=zeros(MnoteCarlo_LSF,MnoteCarlo_SSF,length(Marray)); % Time required
%%
for Mindex=1:length(Marray)
    M = Marray(Mindex);
    tau = M - (K*G) + K ;
    % 
    for MC_LSF=1:MnoteCarlo_LSF
        disp(['M is ',num2str(M),' and LSF is ',num2str(MC_LSF)])
        [Terminal_pos] = Terminal_Position (K,G,r); % Terminal_pos is an array of size (K,G,2)
        [PathLoss] = PathLoss_Genrator(K,G,Terminal_pos); % PathLoss is is an array of size (K,G)
        %
        for MC_SSF=1:MnoteCarlo_SSF
            [channel,BigChannel] = Channel_Generator(M,K,G,PathLoss) ; % channel is (M,K,G) and BigChanne is (M,K*G)
            %% Outer Layer - BDZF Part of the algorithm. Here it is implemented by SVD, it also can be implemented by QR decomposition.
            tic
            [F,channel_Eqe] = SVD_preliminaries(channel,M,K,G,sigma_sqrd,tau);% generate F, equivalent channels, a feasible precoding to start, and the power consumption with this feasible point
            %% Inner Layer - Algorithm 3
            [c_Huristic] = HuristicPrecoder(K,G,eta,tau,channel_Eqe);
            [W_Huri,PowHuristic] = Huristic_PrecandPwrCon(F,c_Huristic,M,G);
            timeHuri(MC_LSF,MC_SSF,Mindex)=toc;
            Pwr_Huri(MC_LSF,MC_SSF,Mindex)=PowHuristic;
            %% Calculate SINR - Verification of methods
            [SINR_Huri] = SINR_Huri_Gen(W_Huri,G,K,channel,sigma_sqrd);
        end
    end
end
%%
figure
MeanPowQoS = zeros(1,length(Marray)); %Lwr bound, sdr, sca, fpp, Huristic
MeanTimeQoS = zeros(1,length(Marray)); %sdr, sca, fpp, Huristic
%%
for Mindex = 1:length(Marray) 
    % POWER
    MeanPowQoS(1,Mindex) = mean(Pwr_Huri(:,:,Mindex));% QR
    % Time
    MeanTimeQoS(1,Mindex) = mean(timeHuri(:,:,Mindex));
end
%%
plot(Marray,MeanPowQoS(1,:),'-*r') % Lower Bound
ylabel('Power Consumption [Watts]')
xlabel('number of antennas')
legend('BDZF + Alg. 3')
title('Power required by BDZF + Alg.3 to meet all the requested SINRs for QoS problem')
grid on
figure
semilogy(Marray,MeanTimeQoS(1,Mindex),'-*')
grid on
title('Time required to find an approximate solution for QoS problem')
ylabel('Time [seconds]')
xlabel('number of antennas')



