function manual_calib(moments_selection, p_selection, cgrid, theta, tg, sp, dm, cw, location, filename)
    [~, ps] = param_selection(theta, p_selection);

    
    fprintf('Running manual calibration for %s\n', filename)
    % Create directory for output
    % folder = fullfile('Calibration_outcomes', filename);  % Create path: 'Calibration_outcomes/filename'
    % if ~isfolder(folder)
    %     mkdir(folder);  % Create the directory if it doesn't exist
    % end

    data_mom=data_moments(dm, moments_selection);
    %Create the combinations
    combinations = create_combinations(cgrid, p_selection,theta);
    nCombinations = size(combinations, 1);
    distances=zeros(nCombinations,1);
    moments=zeros(nCombinations,length(moments_selection));

    parfor idx = 1:nCombinations
        p_vec=[combinations(idx,:)];
        [distances(idx), moments(idx,:)]=wrapperSMM(p_vec, p_selection, theta, ps, tg, sp, dm, moments_selection,cw);
        % Display every 10th iteration
        if mod(idx, 10) == 0
            fprintf('Iteration %d of %d\n', idx, nCombinations)
        end
    end

    % save('calibrationv0.mat','distances','moments','combinations','data_mom')
    % Save data to 'Calibration_outcomes/filename'
   
    if location=="hpc"
        data_mom=data_moments(dm, moments_selection);
        save('/home/gst247/HPC_Model_Matlab/Calibration_outcomes/'+string(filename)+'/manual_calib'+string(filename)+'.mat','distances','moments','combinations','data_mom','p_selection','moments_selection');
        folder = '/home/gst247/HPC_Model_Matlab/Calibration_outcomes/'+string(filename);
        calibration_report(p_selection,ps ,moments_selection, combinations, data_mom, moments, distances, 0.05, folder, filename)
    else
        data_mom=data_moments(dm, moments_selection);
        save('calibrationv0.mat', 'distances', 'moments', 'combinations', 'data_mom','p_selection','moments_selection');
        folder = 'Calibration_outcomes/'+string(filename);
        calibration_report(p_selection,ps ,moments_selection, combinations, data_mom, moments, distances, 0.05, folder, filename)
    end

    fprintf('Calibration for %s done\n', filename)