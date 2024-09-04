%Function that runs, simulate compute moments of the model and get the difference with the data moments
function [distance,mm]=f_SMM(p,ps,tg,sp,dm)
    %Solve the model
    [v,e,w]=joint_loop(p,ps,tg,"spec3");
    save('solved_model.mat','v','e','w')
    % Simulations 
    fs=SimulateFirm_cl(p,ps,tg,sp,v,e,w);
    ws=SimulateWorker_cl(p,ps,sp,tg,v,e,w,fs);
    % save('simulations.mat','fs','ws','sp')
    %Model moments
    mm=model_moments(ps,sp,fs,ws);
    numMoments = length(mm);

    %Difference
    diff=dm-mm;
    %Weightning matrix, identity matrix for now
    WMat= eye(numMoments);

    %Distance
    distance=diff*WMat*diff';






    