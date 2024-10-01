function min_jobs(job_name)
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Min Calibration Jobs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    if strcmp(job_name, 'min_cp_al')
        p_selection={'cost_p','alpha_m'};
        moments_selection = {'NiMi', 'ManWorkerRatio'};
        lb.cost_p=0.01;
        ub.cost_p=0.5;
        lb.alpha_m=0.6;
        ub.alpha_m=0.9;
        %Initial guess
        pvec0=[0.1, theta.alpha_m];
        min_routine(pvec0,p_selection, theta, lb,ub ,tg, sp, dm, moments_selection,cw,job_name)
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
    elseif strcmp(job_name, 'min_alpha')
        tg.use_guess='no';
        p_selection={'alpha_m'};
        moments_selection = {'ManWorkerRatio'};
        cw.ManWorkerRatio=1; 
        lb.alpha_m=0.6;
        ub.alpha_m=0.95;
        %Inital guess
        pvec0=[0.6654];
        min_routine(pvec0,p_selection, theta, lb,ub ,tg, sp, dm, moments_selection,cw,job_name)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif strcmp(job_name, 'min_al_Bcross')
        tg.use_guess='no';
        %Selections
        p_selection={'alpha_m'};
        moments_selection = {'bcross'};
        cw.bcross=1;
        %Bounds
        lb.alpha_m=0.6;
        ub.alpha_m=0.95;
        %Inital guess
        pvec0=[theta.alpha_m];
        min_routine(pvec0,p_selection, theta, lb,ub ,tg, sp, dm, moments_selection,cw,job_name)
    else 
        disp('Job name not recognized')
    end





