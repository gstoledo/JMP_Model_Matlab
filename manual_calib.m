function manual_calib(moments_selection, p_selection, cgrid, theta, tg, sp, dm,location,filename)
    [p_vec, ps] = param_selection(theta, p_selection);
    
    fprintf('Running manual calibration for %s\n', filename)
    data_mom=data_moments(dm, moments_selection);
    %Create the combinations
    combinations = create_combinations(cgrid, p_selection,theta);
    nCombinations = size(combinations, 1);
    distances=zeros(nCombinations,1);
    moments=zeros(nCombinations,length(moments_selection));

    parfor idx = 1:nCombinations
        p_vec=[combinations(idx,:)];
        [distances(idx), moments(idx,:)]=wrapperSMM(p_vec, p_selection, theta, ps, tg, sp, dm, moments_selection);
    end

    save('calibrationv0.mat','distances','moments','combinations','data_mom')
    if location=="hpc"
        data_mom=data_moments(dm, moments_selection);
        save('/home/gst247/HPC_Model_Matlab/Calibration_outcomes/manual_calib'+string(filename)+'.mat','distances','moments','combinations','data_mom')
        calibration_report(p_selection, moments_selection, combinations,data_mom ,moments, distances, 0.05,filename)
    end
    data_mom=data_moments(dm, moments_selection);
    calibration_report(p_selection, moments_selection, combinations,data_mom ,moments, distances, 0.05, filename)
    fprintf('Calibration for %s done\n', filename)