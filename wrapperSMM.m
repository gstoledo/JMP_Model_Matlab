%Function to wrap the SMM function
function [distance,model_mom] = wrapperSMM(p_vec, p_selection, theta, ps, tg, sp, dm, moments_selection)
    %Unpack the parameters
    %Change this one as we add more parameters to the SMM 
    p=p_vec_to_struct(p_vec, p_selection, theta);
    
    %Run the SMM function
    [distance,model_mom] = f_SMM(p,ps,tg,sp,dm,moments_selection);
end