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


%Attempt 0 at calibration for very few parameters

% Learning parameters q_up
%Cost of demotion and promotion

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Toggles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tg.use_guess='y'; 
tg.zero_tol=1e-10;                  %Tolerance for zero
tg.tr=0;                            %Manager penalty toggle
tg.speed=1;                         %Convergenge speed
tg.update_speed_v=1;                %Update speed for value functions
tg.update_speed=1;                  %Update speed for wages
tg.display='y';                      %Display iterations

%% Simulation parameters
sp.seed=1;                          %Seed for random number generator
sp.n_firms=20000;                      %Number of firms
sp.n_workers=sp.n_firms;             %Number of workers
sp.n_months=240;                     %Number of months
sp.burn_years=1;                  %Number of years to burn
sp.t_burn=12*sp.burn_years;        %Number of months to burn
sp.n_years=sp.n_months/12;        %Number of years total
sp.stable_years = sp.n_years - sp.burn_years; %Number of years to keep the simulation


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pre-set parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ps.tpts = 7;                        %Type points space
ps.ats = 5;                         %Productivity space
if location=="hpc"
    ps.ats=5;
    ps.tpts=7;
end
ps.wpts=100;                        %Wage type
ps.bt=1/(1.1^(1/12));               %beta, discount factor
ps.death=1/(35*12);                 %Probability agent dies
ps.bpw=1-xmin(11);                  %bpw is bargaining power of worker
ps.bpf=1-ps.bpw;                    %Firm bargaining power (1-gamma)
ps.n=1;                             %Measure of firms

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ps.aup=0.1;                          %Probability of moving up
ps.adown=0.3;                        %Probability of moving down

ps.qdown=0*ones(ps.ats,1);        %Probability of moving down in q
ps.ulose=x(6);                       %Probability of moving donw in u
ps.A=x(14);                          %Aggegate productivity parameter, still unsure hor to add this here
ps.fcomp=x(3);                       %Complementarity in F
ps.alpha_m=0.7;                      %Manager share
ps.homeprod=x(12);                   %Home production
ps.mnew=x(9);                        %Paramter of exponetial dist for newborn workers
ps.lamu=x(1);                        %Prob of U find a firm
ps.lam=x(2);                         %Prob of firm find another firm
ps.del=x(5);                         %Separation rate

ps.qup=0.3*ones(ps.ats,1);           %Probability of moving up in q. The qup for the last level of productivity is irrelevant

%% Main parameters (to be calibrated)
p.cost_d=1;                         %cost of demoting a manager to non manager  
p.cost_p=1;                         %cost of promoting a non manager to manager




% %Function to wrap the SMM function
% [v,e,w]=joint_loop(p,ps,tg,"spec3");
% save('solved_model.mat','v','e','w')

% % Simulations 
% fs=SimulateFirm_cl(p,ps,tg,sp,v,e,w);
% ws=SimulateWorker_cl(p,ps,sp,tg,v,e,w,fs);
% % save('simulations.mat','fs','ws','sp')
% %Model moments
% mm=model_moments(ps,sp,fs,ws);
% numMoments = length(mm);
% distance=wrapperSMM(p_vec, ps, tg, sp, dm);


% lb = [0.5; 0.5 ; 0.1; 0.1 ];
% ub = [2; 2; 0.5; 0.5];
% % %Number of parameters
% % nvars=length(p_vec);
% % options = optimoptions('ga', 'Display', 'iter', ...
% %                        'UseParallel', true, 'PopulationSize', 100);

% %Run the optimization
% tic;
% [theta_opt, fval] = ga(@(p_vec) wrapperSMM(p_vec, ps, tg, sp, dm), nvars, ...
%                 [], [], [], [], lb, ub, [], [1 2 3 4], options);
% toc;               


%Lets try ionly with the 2 costs
%Import data moments
run data_moments.m
grid_cost_d=linspace(0.5,1.5,5);
grid_cost_p=linspace(0.5,1.5,5);

% parpool('local');
[CostD, CostP] = ndgrid(grid_cost_d, grid_cost_p);
combinations = [CostD(:), CostP(:)];
nCombinations = size(combinations, 1);
distances=zeros(nCombinations,1);
moments=zeros(nCombinations,2);

parfor idx = 1:nCombinations
    % Get the current combination
    cost_d = combinations(idx, 1);
    cost_p = combinations(idx, 2);
    p_vec=[cost_d; cost_p];
    [distances(idx), moments(idx,:)]=wrapperSMM(p_vec, ps, tg, sp, dm);
end
save('calibrationv0.mat','distances','moments','combinations','dm')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Lets get a sense of how these parameters affect the model
%Create a grid (10 points for each parameter) arounf the inital "guesses"
% grid_cost_d=linspace(0.5,1.5,5);
% grid_cost_p=linspace(0.5,1.5,5);
% grid_qup=linspace(0.1,0.5,5);


% parpool('local');
% if ps.ats==2
%     [CostD, CostP, Qup1, Qup2 ] = ndgrid(grid_cost_d, grid_cost_p, grid_qup, grid_qup);
%     combinations = [CostD(:), CostP(:), Qup1(:), Qup2(:)];
%     nCombinations = size(combinations, 1);
%     distances=zeros(nCombinations,1);
%     moments=zeros(nCombinations,4);
%     parfor idx = 1:nCombinations
%         % Get the current combination
%         cost_d = combinations(idx, 1);
%         cost_p = combinations(idx, 2);
%         qup1 = combinations(idx, 3);
%         qup2 = combinations(idx, 4);
    
%         p_vec=[cost_d; cost_p; qup1; qup2];
%         [distances(idx), moments(idx,:)]=wrapperSMM(p_vec, ps, tg, sp, dm);
%     end

% elseif ps.ats==5
%     [CostD, CostP, Qup1, Qup2, Qup3, Qup4, Qup5 ] = ndgrid(grid_cost_d, grid_cost_p, grid_qup, grid_qup, grid_qup, grid_qup, grid_qup);
%     combinations = [CostD(:), CostP(:), Qup1(:), Qup2(:), Qup3(:), Qup4(:), Qup5(:)];
%     nCombinations = size(combinations, 1);
%     distances=zeros(nCombinations,1);
%     moments=zeros(nCombinations,6); %When it is right has to be the same or more than the number of parameters
%     parfor idx = 1:nCombinations
%         % Get the current combination
%         cost_d =combinations(idx, 1);
%         cost_p =combinations(idx, 2);
%         qup1 = combinations(idx, 3);
%         qup2 = combinations(idx, 4);
%         qup3 = combinations(idx, 5);
%         qup4 = combinations(idx, 6);
%         qup5 = combinations(idx, 7);
%         p_vec=[cost_d; cost_p; qup1; qup2; qup3; qup4; qup5];
%         [distances(idx), moments(idx,:)]=wrapperSMM(p_vec, ps, tg, sp, dm);
%     end
% end
% save('calibrationv0.mat','distances','moments','combinations','dm')


