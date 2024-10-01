location = "local";
% location="hpc";
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load baseline parameters
run baseline_param.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tg.display_iter=1; %Display iterations
tg.display_iter_v=1;
tg.display_iter_dist=1;
% theta.ats=5; %Number of productivity levels
% theta.tpts=7; %Number of type points
theta.qup=0.05*ones(theta.ats,1);       %Probability of moving up in q. The qup for the last level of productivity is irrelevant


% %Function to wrap the SMM function
p_selection={'cost_p','cost_d','lamu','lam','alpha_m'};
[p_vec, ps] = param_selection(theta, p_selection);
p=p_vec_to_struct(p_vec, p_selection, theta);
[a_trans,q_trans,u_trans,fteam,fman,fnman,fe,b,mnew_high,typebirth,wmin,wmax,wgrid] = derivatated_p(p,ps,split_top_bot);
tg.speed=1;
tg.update_speed_v=1;
tg.speed_dist=1;
[v,e,w]=joint_loop(p,ps,tg,"spec3");
save('solved_model.mat','v','e','w')




% % Simulations 
% fs=SimulateFirm_cl(p,ps,tg,sp,v,e,w);
% ws=SimulateWorker_cl(p,ps,sp,tg,v,e,w,fs);
% save('simulations.mat','fs','ws','sp')
% % %Model moments
% mm=model_moments(ps,sp,fs,ws);
% numMoments = length(mm);
%Import data moments

% run data_moments.m


%% Manual calibration for some combinations
%% Only costs
p_selection={'cost_p','cost_d'};
[p_vec, ps] = param_selection(theta, p_selection);
moments_selection = {'NiMi', 'MiNi'}; 
cgrid.cost_d=linspace(0.5,1.4,12);
cgrid.cost_p=linspace(0.01,0.2,12);
cgrid.alpha_m=linspace(0.6,0.9,5);
manual_calib(moments_selection, p_selection, cgrid, theta, tg, sp, dm,cw,location,'costs');

%% Only costs with qup=0.06
p_selection={'cost_p','cost_d'};
theta.qup=0.06*ones(theta.ats,1);       %Probability of moving up in q. The qup for the last level of productivity is irrelevant
[p_vec, ps] = param_selection(theta, p_selection);
moments_selection = {'NiMi', 'MiNi', 'ManWorkerRatio'};
cgrid.cost_d=linspace(0.7,2,12);
cgrid.cost_p=linspace(0.1,0.5,12);
manual_calib(moments_selection, p_selection, cgrid, theta, tg, sp, dm,cw,location,'costsq06');


%% Only costs with qup=0.055
p_selection={'cost_p','cost_d'};
theta.qup=0.055*ones(theta.ats,1);       %Probability of moving up in q. The qup for the last level of productivity is irrelevant
[p_vec, ps] = param_selection(theta, p_selection);
moments_selection = {'NiMi', 'MiNi', 'ManWorkerRatio'};
cgrid.cost_d=linspace(1,3,12);
cgrid.cost_p=linspace(0.2,1,12);
manual_calib(moments_selection, p_selection, cgrid, theta, tg, sp, dm,cw,location,'costsq055');


%% With alpha
p_selection={'cost_p','cost_d','alpha_m'};
[p_vec, ps] = param_selection(theta, p_selection);
moments_selection = {'NiMi', 'MiNi', 'ManWorkerRatio'};
cgrid.cost_d=linspace(0.5,1.4,12);
cgrid.cost_p=linspace(0.01,0.2,12);
cgrid.alpha_m=linspace(0.6,0.9,5);
manual_calib(moments_selection, p_selection, cgrid, theta, tg, sp, dm,cw,location,'costs_alpha');


%% Costs amd Lambadas to get only NiMi and MiNi
p_selection={'cost_p','cost_d','lamu','lam'};
moments_selection = {'NiMi', 'MiNi'};
cgrid.cost_d=linspace(0.5,1.4,10);
cgrid.cost_p=linspace(0.01,0.2,10);
cgrid.lamu=linspace(0.1,0.6,5);
cgrid.lam=linspace(0.1,0.6,5);
manual_calib(moments_selection, p_selection, cgrid, theta, tg, sp, dm,cw,location,'costs_lambdas');


%% Test of Wage moments
% Set qup to zero and alpha to 1 
theta.qup=0.05*ones(theta.ats,1);       %Probability of moving up in q. The qup for the last level of productivity is irrelevant
theta.alpha_m=0.8;                      %Manager share
p_selection={'lamu','lam'};
moments_selection = {'ManWorkerRatio','bcross'};
manual_calib(moments_selection, p_selection, cgrid, theta, tg, sp, dm,cw,location,'No_NM');


%% Focus on the rates
p_selection={'cost_p','cost_d', 'lamu', 'lam', 'alpha_m'};
moments_selection = {'NiMi', 'MiNi', 'ManWorkerRatio'};
%Make grids smaller
cgrid.cost_d=linspace(0.5,1.4,4);
cgrid.cost_p=linspace(0.02,0.2,4);
cgrid.lamu=linspace(0.2,0.6,5);
cgrid.lam=linspace(0.2,0.6,5);
cgrid.alpha_m=linspace(0.5,0.9,2);
manual_calib(moments_selection, p_selection, cgrid, theta, tg, sp, dm,cw,location,'rates');


%% Only qup
p_selection={'qup'};
moments_selection = {'NiMi', 'MiNi'};
theta.cost_p=0.25;
theta.cost_d=1.227;
cgrid.qup = repmat({linspace(0.05, 0.07, 3)},length(theta.qup),1);
manual_calib(moments_selection, p_selection, cgrid, theta, tg, sp, dm,cw,location,'qup');


%% Only lambdas
p_selection={'lamu','lam'};
moments_selection = {'NiMi', 'MiNi', 'ManWorkerRatio'};
theta.qup=0.05*ones(theta.ats,1);       
theta.cost_p=0.25;
theta.cost_d=1.5;
cgrid.lamu=linspace(0.4,0.65,10);
cgrid.lam=linspace(0.2,0.35,10);
manual_calib(moments_selection, p_selection, cgrid, theta, tg, sp, dm,cw,location,'lambdas');



%% Only lambdas with cost_d=2
p_selection={'lamu','lam'};
moments_selection = {'NiMi', 'MiNi', 'ManWorkerRatio'};
theta.qup=0.05*ones(theta.ats,1);       
theta.cost_p=0.25;
theta.cost_d=2.5;
theta.alpha_m=0.75;
cgrid.lamu=linspace(0.4,0.65,10);
cgrid.lam=linspace(0.2,0.35,10);
manual_calib(moments_selection, p_selection, cgrid, theta, tg, sp, dm,cw,location,'lmdcd2');


% %% Lets put alpha to the mix
% p_selection={'alpha_m'};
% moments_selection = {'ManWorkerRatio','b_cross'};
% manual_calib(moments_selection, p_selection, cgrid, theta, tg, sp, dm,cw,location,'alpha_qup_Bcross');





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% With minimization
% p_selection={'cost_p','cost_d', 'lamu', 'lam', 'alpha_m'};
% [p_vec, ps] = param_selection(theta, p_selection);
% moments_selection = {'NiMi', 'MiNi', 'ManWorkerRatio'};
% %Number of parameters
% nvars=length(p_selection);
% lb = [cgrid.cost_p(1); cgrid.cost_d(1); cgrid.lamu(1); cgrid.lam(1); cgrid.alpha_m(1)];
% ub = [cgrid.cost_p(end); cgrid.cost_d(end); cgrid.lamu(end); cgrid.lam(end); cgrid.alpha_m(end)];

% %Initial guess
% pvec0=[0.02,0.82,x(1),x(2), 0.7];
% options = optimoptions('ga', 'Display', 'iter', ...
%                        'UseParallel', true, 'PopulationSize', 10, 'InitialPopulationMatrix', pvec0);

% %Run the optimization
% tic;
% [theta_opt, fval] = ga(@(p_vec) wrapperSMM(p_vec, p_selection, theta, ps, tg, sp, dm, moments_selection), nvars, ...
%                 [], [], [], [], lb, ub, [], [], options);
% toc;
% [v,e,w]=joint_loop(theta_opt,ps,tg,"spec3");               
% save('calibrationv0_mincosts.mat','theta_opt','fval','p_selection','moments_selection','lb','ub','v','e','w')
% if location=="hpc"
%     save('/home/gst247/HPC_Model_Matlab/Calibration_outcomes/min_calib_all.mat','theta_opt','fval','p_selection','moments_selection','lb','ub','v','e','w')
% end

