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
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Load baseline parameters
    run baseline_param.m
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        manual_calib(moments_selection, p_selection, cgrid, theta, tg, sp, dm,cw,location,'costs');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif strcmp(job_name, 'costsq06')
        %% Only costs with qup=0.06
        p_selection={'cost_p','cost_d'};
        theta.qup=0.06*ones(theta.ats,1);       %Probability of moving up in q. The qup for the last level of productivity is irrelevant
        theta.lamu=0.45; %Tweaking a bit the rates
        theta.lam=0.3;
        [p_vec, ps] = param_selection(theta, p_selection);
        moments_selection = {'NiMi', 'MiNi', 'ManWorkerRatio'};
        cgrid.cost_d=linspace(1.5,4,12);
        cgrid.cost_p=linspace(0.01,0.5,12);
        manual_calib(moments_selection, p_selection, cgrid, theta, tg, sp, dm,cw,location,'costsq06');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif strcmp(job_name, 'costsq055')
        %% Only costs with qup=0.055
        p_selection={'cost_p','cost_d'};
        theta.qup=0.055*ones(theta.ats,1);       %Probability of moving up in q. The qup for the last level of productivity is irrelevant
        theta.lamu=0.45; %Tweaking a bit the rates
        theta.lam=0.3;
        [p_vec, ps] = param_selection(theta, p_selection);
        moments_selection = {'NiMi', 'MiNi', 'ManWorkerRatio'};
        cgrid.cost_d=linspace(1.5,4,12);
        cgrid.cost_p=linspace(0.01,1,12);
        manual_calib(moments_selection, p_selection, cgrid, theta, tg, sp, dm,cw,location,'costsq055');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif strcmp(job_name, 'costsq05')
        %% Only costs with qup=0.05
        p_selection={'cost_p','cost_d'};
        theta.qup=0.05*ones(theta.ats,1);   
        theta.lamu=0.45; %Tweaking a bit the rates
        theta.lam=0.3;
        [p_vec, ps] = param_selection(theta, p_selection);
        moments_selection = {'NiMi', 'MiNi', 'ManWorkerRatio'};
        cgrid.cost_d=linspace(1.5,4,12);
        cgrid.cost_p=linspace(0.01,1,12);
        manual_calib(moments_selection, p_selection, cgrid, theta, tg, sp, dm,cw,location,'costsq05');
    elseif strcmp(job_name, 'costsq04')
        %% Only costs with qup=0.04
        p_selection={'cost_p','cost_d'};
        theta.qup=0.04*ones(theta.ats,1);   
        theta.lamu=0.45; %Tweaking a bit the rates
        theta.lam=0.3;
        [p_vec, ps] = param_selection(theta, p_selection);
        moments_selection = {'NiMi', 'MiNi', 'ManWorkerRatio'};
        cgrid.cost_d=linspace(1.5,4,12);
        cgrid.cost_p=linspace(0.01,1,12);
        manual_calib(moments_selection, p_selection, cgrid, theta, tg, sp, dm,cw,location,'costsq04');
    elseif strcmp(job_name, 'costsq02')
        %% Only costs with qup=0.02
        p_selection={'cost_p','cost_d'};
        theta.qup=0.02*ones(theta.ats,1);
        theta.lamu=0.475; %Tweaking a bit the rates
        theta.lam=0.255;
        cgrid.cost_d=linspace(0.5,3.3,12);
        cgrid.cost_p=linspace(0.3,1.2,12);
        [p_vec, ps] = param_selection(theta, p_selection);
        moments_selection = {'NiMi', 'MiNi', 'ManWorkerRatio'};
        manual_calib(moments_selection, p_selection, cgrid, theta, tg, sp, dm,cw,location,'costsq02');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % elseif strcmp(job_name, 'costs_alpha')
    %     %% With alpha
    %     p_selection={'cost_p','cost_d','alpha_m'};
    %     [p_vec, ps] = param_selection(theta, p_selection);
    %     moments_selection = {'NiMi', 'MiNi', 'ManWorkerRatio'};
    %     cgrid.cost_d=linspace(0.5,1.4,12);
    %     cgrid.cost_p=linspace(0.01,0.2,12);
    %     cgrid.alpha_m=linspace(0.6,0.9,5);
    %     manual_calib(moments_selection, p_selection, cgrid, theta, tg, sp, dm,cw,location,'costs_alpha');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif strcmp(job_name, 'costs_lambdas')
        %% Costs amd Lambadas to get only NiMi and MiNi
        p_selection={'cost_p','cost_d','lamu','lam'};
        moments_selection = {'NiMi', 'MiNi'};
        cgrid.cost_d=linspace(0.5,1.4,10);
        cgrid.cost_p=linspace(0.01,0.2,10);
        cgrid.lamu=linspace(0.1,0.6,5);
        cgrid.lam=linspace(0.1,0.6,5);
        manual_calib(moments_selection, p_selection, cgrid, theta, tg, sp, dm,cw,location,'costs_lambdas');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif strcmp(job_name, 'cp_alpha')
        %% Promotion costs and alpha
        p_selection={'cost_p','alpha_m'};
        moments_selection = {'NiMi', 'MiNi', 'ManWorkerRatio'};
        theta.qup=0.045*ones(theta.ats,1);
        theta.cost_d=3.25;
        cgrid.cost_p=linspace(0.01,0.4,12);
        cgrid.alpha_m=linspace(0.63,0.9,12);
        manual_calib(moments_selection, p_selection, cgrid, theta, tg, sp, dm,cw,location,'cp_alpha');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif strcmp(job_name, 'qup')
        %% Only qup
        p_selection={'qup'};
        moments_selection = {'NiMi', 'MiNi'};
        theta.cost_p=0.3;
        theta.cost_d=1.227;
        cgrid.qup = repmat({linspace(0.05, 0.07, 3)},length(theta.qup),1);
        manual_calib(moments_selection, p_selection, cgrid, theta, tg, sp, dm,cw,location,'qup');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif strcmp(job_name, 'lambdas')
        %% Only lambdas
        p_selection={'lamu','lam'};
        moments_selection = {'NiMi', 'MiNi', 'ManWorkerRatio'};
        theta.qup=0.05*ones(theta.ats,1);       
        theta.cost_p=0.25;
        theta.cost_d=1.5;
        cgrid.lamu=linspace(0.4,0.65,10);
        cgrid.lam=linspace(0.2,0.35,10);
        manual_calib(moments_selection, p_selection, cgrid, theta, tg, sp, dm,cw,location,'lambdas');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif strcmp(job_name, 'lmdcd2')
        %% Only lambdas with cost_d=2
        p_selection={'lamu','lam'};
        moments_selection = {'NiMi', 'MiNi', 'ManWorkerRatio'};
        theta.qup=0.045*ones(theta.ats,1);       
        theta.cost_p=0.4;
        theta.cost_d=2;
        theta.alpha_m=0.75;
        cgrid.lamu=linspace(0.3,0.45,10);
        cgrid.lam=linspace(0.2,0.35,10);
        manual_calib(moments_selection, p_selection, cgrid, theta, tg, sp, dm,cw,location,'lmdcd2');
    elseif strcmp(job_name, 'lmdcd1')
        %% Only lambdas with cost_d=1
        p_selection={'lamu','lam'};
        moments_selection ={ 'NiMi', 'MiNi', 'ManWorkerRatio'};
        theta.qup=0.04*ones(theta.ats,1);
        theta.cost_p=0.4;
        theta.cost_d=1;
        theta.alpha_m=0.65;
        cgrid.lamu=linspace(0.3,0.45,10);
        cgrid.lam=linspace(0.2,0.35,10);
        manual_calib(moments_selection, p_selection, cgrid, theta, tg, sp, dm,cw,location,'lmdcd1');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif strcmp(job_name, 'cost_a05')
        %% Costs tweatking a_up
        p_selection={'cost_p','cost_d'};
        moments_selection = {'NiMi', 'MiNi', 'ManWorkerRatio'};
        theta.qup=0.05*ones(theta.ats,1);
        theta.aup=0.05; %Tweaking a bit the rates
        theta.adown=0.06;
        theta.lamu=0.475; %Tweaking a bit the rates
        theta.lam=0.255;
        cgrid.cost_d=linspace(0.5,3.3,10);
        cgrid.cost_p=linspace(0.5,1.5,10);
        manual_calib(moments_selection, p_selection, cgrid, theta, tg, sp, dm,cw,location,'cost_a05');
    elseif strcmp(job_name, 'alpha_A')
        %% Costs tweatking a_up
        p_selection={'alpha_m', 'A'};
        moments_selection = {'NiMi', 'MiNi', 'ManWorkerRatio'};
        theta.qup=0.045*ones(theta.ats,1);
        theta.cost_p=0.3;
        theta.cost_d=3.227;
        cgrid.alpha_m=linspace(0.6,0.8,10);
        cgrid.A=linspace(1,3,10);
        manual_calib(moments_selection, p_selection, cgrid, theta, tg, sp, dm,cw,location,'alpha_A');
    elseif strcmp(job_name, 'alpha_rho')
        %% Costs tweatking a_up
        p_selection={'alpha_m', 'fcomp'};
        moments_selection = {'NiMi', 'MiNi', 'ManWorkerRatio'};
        theta.qup=0.045*ones(theta.ats,1);
        theta.cost_p=0.3;
        theta.cost_d=3.227;
        cgrid.alpha_m=linspace(0.6,0.8,10);
        cgrid.fcomp=linspace(0.7,1.2,10);
        manual_calib(moments_selection, p_selection, cgrid, theta, tg, sp, dm,cw,location,'alpha_rho');
    else
        error('Unknown job name: %s', job_name);
    end



