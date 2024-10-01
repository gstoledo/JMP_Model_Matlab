function min_routine(pvec0,p_selection, theta, lb,ub ,tg, sp, dm, moments_selection,cw,job_name)
    nvars=length(p_selection);
    [~, ps] = param_selection(theta, p_selection);

    % Convert structure to a cell array
    cellArray = struct2cell(lb);
    lb_vec = cell2mat(cellArray)';
    cellArray = struct2cell(ub);
    ub_vec = cell2mat(cellArray)';  

    options = optimoptions('ga', ...
            'Display', 'iter', ...               % Display each generation's output
            'PopulationSize', 32, ...            % Set population size to 50
            'InitialPopulationMatrix', pvec0, ... % Use initial guess
            'MaxGenerations', 100, ...           % Maximum of 100 generations
            'MaxStallGenerations', 20, ...       % Stop if no improvement for 20 generations
            'CrossoverFraction', 0.8, ...        % 80% of the population are crossover children
            'SelectionFcn', @selectiontournament, ... % Use tournament selection
            'CrossoverFcn', @crossoverintermediate, ... % Use intermediate crossover
            'MutationFcn', @mutationadaptfeasible, ... % Use adaptive feasible mutation
            'UseParallel', true); % Use two plot functions
    %Run the optimization
    [pvec_opt, fval] = ga(@(p_vec) wrapperSMM(p_vec, p_selection, theta, ps, tg, sp, dm, moments_selection,cw), nvars, ...
        [], [], [], [], lb_vec, ub_vec, [], [], options);

    p_opt=p_vec_to_struct(pvec_opt, p_selection, theta);
    [v,e,w]=joint_loop(p_opt,ps,tg,"spec3");
    fs=SimulateFirm_cl(p_opt,ps,tg,sp,v,e,w);
    ws=SimulateWorker_cl(p_opt,ps,tg,sp,v,e,w,fs);
    mm_vec= model_moments(ps,sp,fs,ws,moments_selection);
    mm=p_vec_to_struct(mm_vec,moments_selection,theta);



    hpc_path='/home/gst247/HPC_Model_Matlab/Calibration_outcomes/'+string(job_name);
    save(fullfile(hpc_path, 'min_calib.mat'),'p_opt','fval','mm','p_selection','moments_selection','lb','ub','v','e','w')
    graphs_selected_model(p_opt,ps,v,e,w,job_name)
    
    fprintf('\n');
    fprintf('The optimal parameters are:\n');
    fields = fieldnames(p_opt);  % Get all field names in the struct
    for i = 1:length(fields)
        fprintf('%s: %.4f\n', fields{i}, p_opt.(fields{i}));
    end
    fprintf('\n');

    fprintf('The model moments are:\n');
    fields = fieldnames(mm);  % Get all field names in the struct
    for i = 1:length(fields)
        fprintf('%s: %.4f\n', fields{i}, mm.(fields{i}));
    end
    
    fprintf('\n');
    fprintf('The data moments are:\n');
    for i = 1:length(moments_selection)
        fprintf('%s: %.4f\n', moments_selection{i}, dm.(moments_selection{i}));
    end
    
    fprintf('The objective function value is: %.4f\n', fval);
    
