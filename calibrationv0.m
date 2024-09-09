location = "local";
location="hpc";
if location == "local"
    %Clear all but location variable
    clearvars -except location
    close all
    cd('/Users/gabrieltoledo/Library/CloudStorage/GoogleDrive-gabrielstoledo.gt@gmail.com/My Drive/PHD NYU/Labor Firm Structure/JMP_Model_Matlab')
    addpath('/Users/gabrieltoledo/Library/CloudStorage/GoogleDrive-gabrielstoledo.gt@gmail.com/My Drive/PHD NYU/Labor Firm Structure/JMP_Model_Matlab')
    addpath('/Users/gabrieltoledo/Library/CloudStorage/GoogleDrive-gabrielstoledo.gt@gmail.com/My Drive/PHD NYU/Labor Firm Structure/Python_Firm Structure/replication_HLMP/matlab/run_A4_revisit')
else
    % Add your working directory to MATLAB's path
    addpath('/home/gst247/HPC_Model_Matlab/JMP_Model_Matlab')
end
load('xmin_results_combined.mat')
x=xmin; %Use some values of HLPM to start 


%Attempt 0 at calibration for very few parameters

% Learning parameters q_up
%Cost of demotion and promotion

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Toggles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tg.use_guess='no'; 
tg.zero_tol=1e-10;                  %Tolerance for zero
tg.tr=0;                            %Manager penalty toggle
tg.speed=1;                         %Convergenge speed
tg.update_speed_v=1;                %Update speed for value functions
tg.update_speed=1;                  %Update speed for wages
tg.display='y';                      %Display iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulation parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sp.seed=1;                          %Seed for random number generator
sp.n_firms=20000;                      %Number of firms
sp.n_workers=sp.n_firms;             %Number of workers
sp.n_months=240;                     %Number of months
sp.burn_years=1;                  %Number of years to burn
sp.t_burn=12*sp.burn_years;        %Number of months to burn
sp.n_years=sp.n_months/12;        %Number of years total
sp.stable_years = sp.n_years - sp.burn_years; %Number of years to keep the simulation


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data Moments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dm.NiMi                             =0.047+0.009;  %Rate internal promotions to managers (as a fraction of total managers)
dm.MiNi                             =0.021+0.001;  %Rate internal demotions to non-managers (as a fraction of total non-managers)
dm.Q5Q1nm_wage                      =6.84;          %Non Manager Wage Ratio quintile 5 to quintile 1 of total avg wage
dm.Q5Q2nm_wage                      =3.35;          %Non Manager Wage Ratio quintile 5 to quintile 2 of total avg wage
dm.Q5Q3nm_wage                      =2.30;          %Non Manager Wage Ratio quintile 5 to quintile 3 of total avg wage
dm.Q5Q4nm_wage                      =1.68;          %Non Manager Wage Ratio quintile 5 to quintile 4 of total avg wage
dm.Q5Q5nm_wage                      =1;           %Non Manager Wage Ratio quintile 5 to quintile 4 of total avg wage

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pre-set parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta.tpts = 3;                        %Type points space
theta.ats = 2;                         %Productivity space
if location == "hpc"
    theta.ats=5;
    theta.tpts=7;
end
theta.wpts=100;                        %Wage type
theta.bt=1/(1.1^(1/12));               %beta, discount factor
theta.death=1/(35*12);                 %Probability agent dies
theta.bpw=1-xmin(11);                  %bpw is bargaining power of worker
theta.bpf=1-theta.bpw;                    %Firm bargaining power (1-gamma)
theta.n=1;                             %Measure of firms

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theta.aup=0.1;                          %Probability of moving up
theta.adown=0.3;                        %Probability of moving down
theta.qdown=0*ones(theta.ats,1);        %Probability of moving down in q
theta.ulose=x(6);                       %Probability of moving donw in u
theta.A=x(14);                          %Aggegate productivity parameter, still unsure hor to add this here
theta.fcomp=x(3);                       %Complementarity in F
theta.alpha_m=0.7;                      %Manager share
theta.homeprod=x(12);                   %Home production
theta.mnew=x(9);                        %Paramter of exponetial dist for newborn workers
theta.lamu=x(1);                        %Prob of U find a firm
theta.lam=x(2);                         %Prob of firm find another firm
theta.del=x(5);                         %Separation rate
theta.qup=0.3*ones(theta.ats,1);           %Probability of moving up in q. The qup for the last level of productivity is irrelevant
theta.cost_d=0.8;                         %cost of demoting a manager to non manager  
theta.cost_p=0.3;                         %cost of promoting a non manager to manager

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Grids
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cgrid.cost_d=linspace(0.5,1.4,15);
cgrid.cost_p=linspace(0.05,0.4,15);
cgrid.qup = repmat({linspace(0.1, 0.5, 5)},length(theta.qup),1);
cgrid.lamu=linspace(0.1,0.6,10);
cgrid.lam=linspace(0.1,0.6,10);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
%Import data moments

% run data_moments.m

% tic
% distance=wrapperSMM(p_vec, ps, tg, sp, dm);
% toc;


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




%% Manual calibration for some combinations
%Only costs
p_selection={'cost_p','cost_d',};
moments_selection = {'NiMi', 'MiNi'};
cgrid.cost_d=linspace(0.5,1.4,15);
cgrid.cost_p=linspace(0.01,0.2,15);
manual_calib(moments_selection, p_selection, cgrid, theta, tg, sp, dm,location,'costs_lowcp');

% %Only qup
% p_selection={'qup'};
% moments_selection = {'Q5Q1nm_wage', 'Q5Q2nm_wage', 'Q5Q3nm_wage', 'Q5Q4nm_wage', 'Q5Q5nm_wage'};
% manual_calib(moments_selection, p_selection, cgrid, theta, tg, sp, dm,location,'qup');


%Costs amd Lambadas to get only NiMi and MiNi
p_selection={'cost_p','cost_d','lamu','lam'};
moments_selection = {'NiMi', 'MiNi'};
manual_calib(moments_selection, p_selection, cgrid, theta, tg, sp, dm,location,'costs_lambdas');





% % With minimization
% p_selection={'cost_p','cost_d',};
% [p_vec, ps] = param_selection(theta, p_selection);
% moments_selection = {'NiMi', 'MiNi'};
% %Number of parameters
% nvars=length(p_selection);
% lb = [cgrid.cost_p(1); cgrid.cost_d(1)];
% ub = [cgrid.cost_p(end); cgrid.cost_d(end)];

% %Initial guess
% pvec0=[0.1,1];
% options = optimoptions('ga', 'Display', 'iter', ...
%                        'UseParallel', true, 'PopulationSize', 10, 'InitialPopulationMatrix', pvec0);

% %Run the optimization
% tic;
% [theta_opt, fval] = ga(@(p_vec) wrapperSMM(p_vec, p_selection, theta, ps, tg, sp, dm, moments_selection), nvars, ...
%                 [], [], [], [], lb, ub, [], [], options);
% toc;               
% save('calibrationv0_mincosts.mat','theta_opt','fval')
% if location=="hpc"
%     save('/home/gst247/HPC_Model_Matlab/Calibration_outcomes/min_calib.mat','theta_opt','fval')
% end