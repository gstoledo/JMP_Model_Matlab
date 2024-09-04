%Function to wrap the SMM function
function [dist, mm] = wrapperSMM(paramVec, ps, tg, sp, dm)
    %Unpack the parameters
    %Change this one as we add more parameters to the SMM 
    p = struct();
    p.cost_d = paramVec(1);
    p.cost_p = paramVec(2);
    % p.qup = paramVec(3)*ones(ps.ats,1);
    
    %Run the SMM function
    [dist,mm] = f_SMM(p,ps,tg,sp,dm);
end