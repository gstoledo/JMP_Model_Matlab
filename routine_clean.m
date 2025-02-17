clear all
close all
location="local";
% location="hpc";
if location == "local"
    cd('/Users/gabrieltoledo/Library/CloudStorage/GoogleDrive-gabrielstoledo.gt@gmail.com/My Drive/PHD NYU/Labor Firm Structure/JMP_Model_Matlab')
    addpath('/Users/gabrieltoledo/Library/CloudStorage/GoogleDrive-gabrielstoledo.gt@gmail.com/My Drive/PHD NYU/Labor Firm Structure/JMP_Model_Matlab')
    addpath('/Users/gabrieltoledo/Library/CloudStorage/GoogleDrive-gabrielstoledo.gt@gmail.com/My Drive/PHD NYU/Labor Firm Structure/Python_Firm Structure/replication_HLMP/matlab/run_A4_revisit')
end
load('xmin_results_combined.mat')
x=xmin; %Use some values of HLPM to start 




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Toggles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tg.use_guess='y'; 
tg.zero_tol=1e-10;                  %Tolerance for zero
tg.nm_penal=0;                            %Manager penalty toggle
tg.speed=1;                         %Convergenge speed
tg.update_speed_v=1;                %Update speed for value functions
tg.update_speed=1;                  %Update speed for wages

%% Simulation parameters
sp.seed=1;                          %Seed for random number generator
sp.n_firms=10;                      %Number of firms
sp.n_workers=sp.n_firms;             %Number of workers
sp.n_months=24;                     %Number of months
sp.burn_years=1;                  %Number of years to burn
sp.t_burn=12*sp.burn_years;        %Number of months to burn
sp.n_years=sp.n_months/12;        %Number of years total
sp.stable_years = sp.n_years - sp.burn_years; %Number of years to keep the simulation


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pre-set parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ps.tpts = 3;                        %Type points space
ps.ats = 2;                         %Productivity space
ps.wpts=100;                        %Wage type
ps.bt=1/(1.1^(1/12));               %beta, discount factor
ps.death=1/(35*12);                 %Probability agent dies
ps.bpw=1-xmin(11);                  %bpw is bargaining power of worker
ps.bpf=1-ps.bpw;                    %Firm bargaining power (1-gamma)
ps.n=1;                             %Measure of firms

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Main parameters (to be calibrated)
p.aup=0.1;                          %Probability of moving up
p.adown=0.3;                        %Probability of moving down
p.qup=0.3*ones(ps.ats,1);           %Probability of moving up in q
p.qdown=0.25*ones(ps.ats,1);        %Probability of moving down in q
p.ulose=x(6);                       %Probability of moving donw in u
p.A=x(14);                          %Aggegate productivity parameter, still unsure hor to add this here
p.fcomp=x(3);                       %Complementarity in F
p.alpha_m=0.7;                      %Manager share
p.homeprod=x(12);                   %Home production
p.mnew=x(9);                        %Paramter of exponetial dist for newborn workers
p.lamu=x(1);                        %Prob of U find a firm
p.lam=x(2);                         %Prob of firm find another firm
p.del=x(5);                         %Separation rate
p.cost_d=1;                         %cost of demoting a manager to non manager  
p.cost_p=1;                         %cost of promoting a non manager to manager



%Joint loop 
tic;
[v,e,w]=joint_loop(p,ps,tg,"spec3");
toc;

% Simulations
n_months=30;
n_firms=10;
n_workers=n_firms;
fs=SimulateFirm_cl(p,ps,tg,sp,v,e,w);
ws=SimulateWorker_cl(p,ps,sp,tg,v,e,w,fs);

save('simulations.mat','fs','ws','sp')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




