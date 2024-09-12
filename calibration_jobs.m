function calibration_jobs(job_name)
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
    tg.use_guess='no';                %Use guess for value functions
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
    dm.ManWorkerRatio                   =1.53;           % Ratio between managers and workers wage/employment ratio
    dm.bcross                           =0.137;         %Cross Section Wage Premium

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
    theta.qup=0.05*ones(theta.ats,1);       %Probability of moving up in q. The qup for the last level of productivity is irrelevant
    theta.cost_d=0.8;                       %cost of demoting a manager to non manager  
    theta.cost_p=0.02;                       %cost of promoting a non manager to manager

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Grids
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cgrid.cost_d=linspace(0.5,1.4,15);
    cgrid.cost_p=linspace(0.02,0.2,15);
    cgrid.qup = repmat({linspace(0.08, 0.4, 5)},length(theta.qup),1);
    cgrid.lamu=linspace(0.2,0.4,10);
    cgrid.lam=linspace(0.2,0.4,10);
    cgrid.alpha_m=linspace(0.5,0.9,10);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calibration jobs 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(job_name, 'costs')
        %% Only costs
        p_selection={'cost_p','cost_d'};
        [p_vec, ps] = param_selection(theta, p_selection);
        moments_selection = {'NiMi', 'MiNi'}; 
        cgrid.cost_d=linspace(0.5,1.4,12);
        cgrid.cost_p=linspace(0.01,0.2,12);
        cgrid.alpha_m=linspace(0.6,0.9,5);
        manual_calib(moments_selection, p_selection, cgrid, theta, tg, sp, dm,location,'costs');
    elseif strcmp(job_name, 'costsq06')
        %% Only costs with qup=0.06
        p_selection={'cost_p','cost_d'};
        theta.qup=0.06*ones(theta.ats,1);       %Probability of moving up in q. The qup for the last level of productivity is irrelevant
        [p_vec, ps] = param_selection(theta, p_selection);
        moments_selection = {'NiMi', 'MiNi', 'ManWorkerRatio'};
        cgrid.cost_d=linspace(0.7,2,12);
        cgrid.cost_p=linspace(0.1,0.5,12);
        manual_calib(moments_selection, p_selection, cgrid, theta, tg, sp, dm,location,'costsq06');
    elseif strcmp(job_name, 'costsq055')
        %% Only costs with qup=0.055
        p_selection={'cost_p','cost_d'};
        theta.qup=0.055*ones(theta.ats,1);       %Probability of moving up in q. The qup for the last level of productivity is irrelevant
        [p_vec, ps] = param_selection(theta, p_selection);
        moments_selection = {'NiMi', 'MiNi', 'ManWorkerRatio'};
        cgrid.cost_d=linspace(1,3,12);
        cgrid.cost_p=linspace(0.2,1,12);
        manual_calib(moments_selection, p_selection, cgrid, theta, tg, sp, dm,location,'costsq055');
    elseif strcmp(job_name, 'costs_alpha')
        %% With alpha
        p_selection={'cost_p','cost_d','alpha_m'};
        [p_vec, ps] = param_selection(theta, p_selection);
        moments_selection = {'NiMi', 'MiNi', 'ManWorkerRatio'};
        cgrid.cost_d=linspace(0.5,1.4,12);
        cgrid.cost_p=linspace(0.01,0.2,12);
        cgrid.alpha_m=linspace(0.6,0.9,5);
        manual_calib(moments_selection, p_selection, cgrid, theta, tg, sp, dm,location,'costs_alpha');
    elseif strcmp(job_name, 'costs_lambdas')
        %% Costs amd Lambadas to get only NiMi and MiNi
        p_selection={'cost_p','cost_d','lamu','lam'};
        moments_selection = {'NiMi', 'MiNi'};
        cgrid.cost_d=linspace(0.5,1.4,10);
        cgrid.cost_p=linspace(0.01,0.2,10);
        cgrid.lamu=linspace(0.1,0.6,5);
        cgrid.lam=linspace(0.1,0.6,5);
        manual_calib(moments_selection, p_selection, cgrid, theta, tg, sp, dm,location,'costs_lambdas');
    elseif strcmp(job_name, 'rates')
        %% Focus on the rates (didnt run)
        p_selection={'cost_p','cost_d', 'lamu', 'lam', 'alpha_m'};
        moments_selection = {'NiMi', 'MiNi', 'ManWorkerRatio'};
        %Make grids smaller
        cgrid.cost_d=linspace(0.5,1.4,4);
        cgrid.cost_p=linspace(0.02,0.2,4);
        cgrid.lamu=linspace(0.2,0.6,5);
        cgrid.lam=linspace(0.2,0.6,5);
        cgrid.alpha_m=linspace(0.5,0.9,2);
        manual_calib(moments_selection, p_selection, cgrid, theta, tg, sp, dm,location,'rates');
    elseif strcmp(job_name, 'qup')
        %% Only qup
        p_selection={'qup'};
        moments_selection = {'NiMi', 'MiNi'};
        theta.cost_p=0.25;
        theta.cost_d=1.227;
        cgrid.qup = repmat({linspace(0.05, 0.07, 3)},length(theta.qup),1);
        manual_calib(moments_selection, p_selection, cgrid, theta, tg, sp, dm,location,'qup');

    elseif strcmp(job_name, 'lambdas')
        %% Only lambdas
        p_selection={'lamu','lam'};
        moments_selection = {'NiMi', 'MiNi', 'ManWorkerRatio'};
        theta.qup=0.05*ones(theta.ats,1);       
        theta.cost_p=0.25;
        theta.cost_d=1.5;
        cgrid.lamu=linspace(0.4,0.65,10);
        cgrid.lam=linspace(0.2,0.35,10);
        manual_calib(moments_selection, p_selection, cgrid, theta, tg, sp, dm,location,'lambdas');
    elseif strcmp(job_name, 'lmdcd2')
        %% Only lambdas with cost_d=2
        p_selection={'lamu','lam'};
        moments_selection = {'NiMi', 'MiNi', 'ManWorkerRatio'};
        theta.qup=0.05*ones(theta.ats,1);       
        theta.cost_p=0.25;
        theta.cost_d=2.5;
        theta.alpha_m=0.75;
        cgrid.lamu=linspace(0.4,0.65,10);
        cgrid.lam=linspace(0.2,0.35,10);
        manual_calib(moments_selection, p_selection, cgrid, theta, tg, sp, dm,location,'lmdcd2');
    else
        error('Unknown job name: %s', job_name);
    end
    


