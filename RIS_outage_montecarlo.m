clc;
clear all;
close all;

%% Define the variables
trial = 1;
gamma_th = 10;
for kappa = [50] % parameter for von-mises phase estimation error distribution
for Kmin = [0.5]
steepness = 1.5; % steepness of the curve between amplititude & phase shift
for N = [10] % number of RIS antenna elements
%kappa = 20; 
%
epcilon = (((1-Kmin)/pi)* beta(steepness+0.5,0.5))+Kmin;
delta = (((1- Kmin)^2)/pi)* beta(2*steepness+0.5,0.5) - (((1- Kmin)^2)/pi^2)* (beta(steepness+0.5,0.5))^2;
mean_Ar = ((N*pi*epcilon)*besseli(1,kappa))/(4*besseli(0,kappa));
mean_Ai = 0;
util = (delta + (epcilon^2)*(1-pi*pi/16));
func = @(x) (exp(kappa*cos(x))).*cos(x).^2;
w = integral(func, -pi,pi);
variance_Ar = N*((util)*((w/(2*pi*besseli(0,kappa)))-(besseli(1,kappa)/besseli(0,kappa))^2) ...
             + (util)*((besseli(1,kappa)/besseli(0,kappa))^2) ...
             + ((w/(2*pi*besseli(0,kappa)))-(besseli(1,kappa)/besseli(0,kappa))^2)*(pi*pi*epcilon*epcilon/16));
func1 = @(x) (exp(kappa*cos(x))).*sin(x).^2;
v = integral(func1, -pi,pi);        
variance_Ai = N*((util)*(v/(2*pi*besseli(0,kappa))) + (v/(2*pi*besseli(0,kappa)))*(pi*pi*epcilon*epcilon/16));
%% generate gaussian Ar and Ai with above mean and variance
outage_trials = zeros(65,trial);
vector = ones(trial,1);
num_sample = 1000000;
for MC_trial=1:trial
    count = 1;
    for transmit_dbm=1.6:.1:8.0
        gamma_t = (10^(transmit_dbm))/1000;
        data = generate_gamma_samples(kappa,Kmin,steepness,N,gamma_t,num_sample);
        less_data = data(data < gamma_th);
        outage_montecarlo(count) = size(less_data,2)/num_sample;
        count = count+1;
    end
    %histogram(data,'normalization','pdf');
    outage_trials(:,MC_trial) = outage_montecarlo;
end
averaged_outage_mc = outage_trials*vector./trial;
% %% Do it using derivation
% count=1;
% for transmit_dbm=1.6:.1:8.0
%     gamma_t = (10^(transmit_dbm))/1000;
%     T = gamma_th/gamma_t;
%     pdf_bsq = @(b) (1/(sqrt(2)*(variance_Ai)*gamma(0.5))).*((b/variance_Ai).^-0.5).*exp(-b/(2*variance_Ai));
%     intpdfbsqa = @(a) integral(pdf_bsq,0,T-a);
%     pdf = @(a) intpdfbsqa(a)*(1/(2*variance_Ar)).*((a/mean_Ar.^2).^-0.25).*exp(-(a+mean_Ar.^2)/(2*variance_Ar)).*(besseli(-0.5,sqrt(a*mean_Ar.^2)/variance_Ar));
%     outage_derivation(count) = integral(pdf,0,T,'ArrayValued', true);
%     count=count+1;
% end
% %%
 semilogy(16:80,averaged_outage_mc,'--o','linewidth', 2);
% %semilogy(16:80,outage_montecarlo,'--o','linewidth', 2);
 hold on
% semilogy(16:80,outage_derivation,'--o','linewidth', 2);
% 
 xlabel('Transmit SNR(dBm)');
 ylabel('Outage Probability');
 xlim([16 80])
% grid on;
% hold on;
 end
 end
 end
% 
% legend('MC Simulation:\kappa=25','Derivation:\kappa=25','MC Simulation:\kappa=50','Derivation:\kappa=50');
%%
function data =generate_gamma_samples(kappa,kmin,steepness,N,gamma_t,out_data_size)

        for k=1:out_data_size
        % i = [1:N]
        % gamma = |Sum i: alpha(i)*beta(i)*rho(i)*exp(j*ph_error(i))|^2*gamma_th
        % generate data elements
        phase_error_samples = vmrand(0,kappa,[N 1]);
        phi = unifrnd(-pi,pi,[N 1]);
        u = phi+pi/2;
        rho = ((1-kmin)*((sin(phi-u)+1)/2).^steepness)+kmin;
        alpha = raylrnd(sqrt(0.5),N,1);
        beta  = raylrnd(sqrt(0.5),N,1);
        % Form gamma by summing all the terms
        A_sq=0;
        for index =1:N
            % generate alpha,beta,phase_error,a uniform rv between -pi to pi
            A_sq = A_sq + alpha(index)*beta(index)*rho(index)*exp(1i*phase_error_samples(index));    
        end
        data(k) = (abs(A_sq)^2)*gamma_t;
        end
end
%%



