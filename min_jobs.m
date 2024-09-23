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
        [~, ps] = param_selection(theta, p_selection);
        moments_selection = {'NiMi', 'ManWorkerRatio'};
        %Number of parameters
        nvars=length(p_selection);
        lb.cost_p=0.01;
        ub.cost_p=0.5;
        lb.alpha_m=0.6;
        ub.alpha_m=0.9;

        lb=[lb.cost_p, lb.alpha_m];
        ub=[ub.cost_p, ub.alpha_m];
        %Initial guess
        pvec0=[0.1, theta.alpha_m];
        options = optimoptions('ga', 'Display', 'iter', ...
                            'UseParallel', true, 'PopulationSize', 10, 'InitialPopulationMatrix', pvec0);
        %Run the optimization
        [pvec_opt, fval] = ga(@(p_vec) wrapperSMM(p_vec, p_selection, theta, ps, tg, sp, dm, moments_selection), nvars, ...
                        [], [], [], [], lb, ub, [], [], options);
        p_opt=p_vec_to_struct(pvec_opt, p_selection, theta);
        [v,e,w]=joint_loop(p_opt,ps,tg,"spec3");
        hpc_path='/home/gst247/HPC_Model_Matlab/Calibration_outcomes/'+string(job_name);    
        save(fullfile(hpc_path, 'min_calib_all.mat'),'p_opt','fval','p_selection','moments_selection','lb','ub','v','e','w')
        graphs_selected_model(p_opt,ps,v,e,w,job_name)
    else 
        disp('Job name not recognized')
    end





