%Function that runs, simulate compute moments of the model and get the difference with the data moments
function [distance,model_mom]=f_SMM(p,ps,tg,sp,dm,moments_selection)
    %Solve the model
    [v,e,w]=joint_loop(p,ps,tg,"spec3");
    save('solved_model.mat','v','e','w')
    % Simulations 
    fs=SimulateFirm_cl(p,ps,tg,sp,v,e,w);
    ws=SimulateWorker_cl(p,ps,sp,tg,v,e,w,fs);
    % save('simulations.mat','fs','ws','sp')
    %Data moments
    data_mom=data_moments(dm, moments_selection);

    %Model moments
    model_mom= model_moments(ps,sp,fs,ws,moments_selection);
    numMoments = length(model_mom);

    %Difference
    diff=model_mom-data_mom;
    %Weightning matrix, identity matrix for now
    WMat= eye(numMoments);

    %Distance
    distance=diff*WMat*diff';






    