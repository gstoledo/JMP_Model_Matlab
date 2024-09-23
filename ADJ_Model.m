location = "local";
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


%Attempt 0 at calibration for very few parameters

% Learning parameters q_up
%Cost of demotion and promotion

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load baseline parameters
run baseline_param.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %Function to wrap the SMM function
p_selection={'cost_p','cost_d','lamu','lam','alpha_m'};
[p_vec, ps] = param_selection(theta, p_selection);
p=p_vec_to_struct(p_vec, p_selection, theta);
p.lamu=0.4; %Set lamu to a fixed value
tic
[v,e,w]=joint_loop(p,ps,tg,"spec3")
toc
save('solved_model.mat','v','e','w')

fs=SimulateFirm_cl(p,ps,tg,sp,v,e,w);
ws=SimulateWorker_cl(p,ps,sp,tg,v,e,w,fs);






load('solved_model.mat')

%% Policy Functions
%Open up the parameters
fieldNames = fieldnames(p);
% Loop over each field and assign it to a variable in the workspace
for i = 1:length(fieldNames)% Dynamically create the variable name
    varName = fieldNames{i};
    % Use eval to assign the value to the variable with the same name
    eval([varName ' = p.' varName ';']);
end
%Opening up the pre set parameters
fieldNames = fieldnames(ps);
% Loop over each field and assign it to a variable in the workspace
for i = 1:length(fieldNames)% Dynamically create the variable name
    varName = fieldNames{i};
    % Use eval to assign the value to the variable with the same name
    eval([varName ' = ps. ' varName ';']);
end
%Open the VFs
fieldNames = fieldnames(v);
% Loop over each field and assign it to a variable in the workspace
for i = 1:length(fieldNames)% Dynamically create the variable name
    varName = fieldNames{i};
    % Use eval to assign the value to the variable with
    eval([varName ' = v.' varName ';']);
end
%Open the distributions
fieldNames = fieldnames(e);
% Loop over each field and assign it to a variable in the workspace
for i = 1:length(fieldNames)% Dynamically create the variable name
    varName = fieldNames{i};
    % Use eval to assign the value to the variable with
    eval([varName ' = e.' varName ';']);
end
%Open the wages
fieldNames = fieldnames(w);
% Loop over each field and assign it to a variable in the workspace
for i = 1:length(fieldNames)% Dynamically create the variable name
    varName = fieldNames{i};
    % Use eval to assign the value to the variable with
    eval([varName ' = w.' varName ';']);
end




%Get the policy functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Policy functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Hiring policies, looking at origins not allocations yet 
[h_e_u, h_e_m, h_e_nm, h_e_t_m, h_e_t_nm, h_m_u, h_m_m, h_m_nm, h_m_t_m, h_m_t_nm, h_nm_u, h_nm_m, h_nm_nm, h_nm_t_m, h_nm_t_nm, h_t_u, h_t_m, h_t_nm, h_t_t_m, h_t_t_nm]...
    =hire_policies(Vmh,Vnh,Veh,U,Vth,ats,tpts,0,cost_d,cost_p);
    
% Allocation Policies after hiring at SM
[p_e_m, p_e_n, p_m_m_u, p_m_m_d, p_m_n, p_n_n_u, p_n_n_p, p_n_m, p_t_m_u, p_t_m_d, p_t_n_u, p_t_n_p]...
    =alloc_policies(Vmh,Vnh,U,Vth,ats,tpts,cost_d,cost_p);


sum(e.eplus_tdist,"all")
sum(e.eplus_ndist,"all")
sum(e.eplus_mdist,"all")
sum(e.eplus_edist,"all")
sum(e.eplus_udist,"all")


