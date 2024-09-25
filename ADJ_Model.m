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




% Learning parameters q_up
%Cost of demotion and promotion

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load baseline parameters
run baseline_param.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Function to wrap the SMM function
tg.display_iter=1;
tg.update_speed=1;
tg.update_speed_v=1;
p_selection={'cost_p','cost_d','lamu','lam','alpha_m'};
[p_vec, ps] = param_selection(theta, p_selection);
ps.phim=0.8;
ps.phin=0.8;
p=p_vec_to_struct(p_vec, p_selection, theta);
[v,e,w]=joint_loopNCC(p,ps,tg,'NCC');
% [v,e,w]=joint_loop(p,ps,tg,'Baseline');




%% Opening up the parameters
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
%Opening up the toggles
fieldNames = fieldnames(tg);
% Loop over each field and assign it to a variable in the workspace
for i = 1:length(fieldNames)% Dynamically create the variable name
    varName = fieldNames{i};
    % Use eval to assign the value to the variable with the same name
    eval([varName ' = tg.' varName ';']);
end

%% Derivated from the main parameters
[a_trans,q_trans,u_trans,fteam,fman,fnman,fe,b,mnew_high,typebirth,wmin,wmax,wgrid] = derivatated_p(p,ps);


Veini = zeros(1, ats);
Vmini = fman / (1 - bt);
Vnini = fnman / (1 - bt);
Vtini = fteam / (1 - bt);
Uini = b / (1 - bt);

phim=0;
phin=0;



[Ve, Vm, Vn, Vt, U, Veh, Vmh, Vnh, Vth, Vetl, Vmtl, Vntl, Vttl, Utl]= vf_iterationNCC(eplus_edist,eplus_mdist,eplus_ndist,eplus_tdist,eplus_udist,Veini,Vmini,Vnini,Vtini,Uini,ats,tpts,phim,phin,cost_d,cost_p,nm_penal,lamu, lam, del, bt, death, bpf, bpw, n, b, fteam,fman,fnman,fe, u_trans, a_trans,q_trans,speed,display_iter_v);
%Test
% e_mdist=eplus_ndist
% e_ndist=eplus_mdist
% e_tdist=permute(eplus_tdist,[1,3,2])

% [e1_udist,e1_edist, e1_mdist, e1_ndist,e1_tdist]=e1_distNCC(Veh,Vmh,Vnh,U,Vth,ats,tpts,phim,phin,cost_d,cost_p,nm_penal,e_udist,e_edist,e_mdist,e_ndist,e_tdist,lamu,lam,n,del);
% [e1_udist,e1_edist, e1_mdist, e1_ndist,e1_tdist]=e1_distV3(Veh,Vmh,Vnh,U,Vth,ats,tpts,cost_d,cost_p,e_udist,e_edist,e_mdist,e_ndist,e_tdist,lamu,lam,n,del);


[e1_udist,e1_edist,e1_mdist,e1_ndist,e1_tdist]=e1_distNCC_phi(Veh,Vmh,Vnh,Uh,Vth,ats,tpts,phim,phin,cost_d,cost_p,e_udist,e_edist,e_mdist,e_ndist,e_tdist,lamu,lam,n,del,nm_penal);

%It is not really getting to one...
sum(e1_edist)+sum(e1_mdist,"all")+sum(e1_ndist,"all")+sum(e1_tdist,"all") % We are probabily missing a flow somewhere, that we will check later 

1-(sum(e1_edist)+sum(e1_mdist,"all")+sum(e1_ndist,"all")+sum(e1_tdist,"all")) % We are probabily missing a flow somewhere, that we will check later

pop1=sum(e1_udist)+sum(e1_mdist,"all")+sum(e1_ndist,"all")+2*sum(e1_tdist,"all")








% % %Function to wrap the SMM function
% p_selection={'cost_p','cost_d','lamu','lam','alpha_m'};
% [p_vec, ps] = param_selection(theta, p_selection);
% p=p_vec_to_struct(p_vec, p_selection, theta);
% p.lamu=0.4; %Set lamu to a fixed value
% tic
% [v,e,w]=joint_loop(p,ps,tg,"spec3")
% toc
% save('solved_model.mat','v','e','w')

% fs=SimulateFirm_cl(p,ps,tg,sp,v,e,w);
% ws=SimulateWorker_cl(p,ps,sp,tg,v,e,w,fs);























% load('solved_model.mat')

% %% Policy Functions
% %Open up the parameters
% fieldNames = fieldnames(p);
% % Loop over each field and assign it to a variable in the workspace
% for i = 1:length(fieldNames)% Dynamically create the variable name
%     varName = fieldNames{i};
%     % Use eval to assign the value to the variable with the same name
%     eval([varName ' = p.' varName ';']);
% end
% %Opening up the pre set parameters
% fieldNames = fieldnames(ps);
% % Loop over each field and assign it to a variable in the workspace
% for i = 1:length(fieldNames)% Dynamically create the variable name
%     varName = fieldNames{i};
%     % Use eval to assign the value to the variable with the same name
%     eval([varName ' = ps. ' varName ';']);
% end
% %Open the VFs
% fieldNames = fieldnames(v);
% % Loop over each field and assign it to a variable in the workspace
% for i = 1:length(fieldNames)% Dynamically create the variable name
%     varName = fieldNames{i};
%     % Use eval to assign the value to the variable with
%     eval([varName ' = v.' varName ';']);
% end
% %Open the distributions
% fieldNames = fieldnames(e);
% % Loop over each field and assign it to a variable in the workspace
% for i = 1:length(fieldNames)% Dynamically create the variable name
%     varName = fieldNames{i};
%     % Use eval to assign the value to the variable with
%     eval([varName ' = e.' varName ';']);
% end
% %Open the wages
% fieldNames = fieldnames(w);
% % Loop over each field and assign it to a variable in the workspace
% for i = 1:length(fieldNames)% Dynamically create the variable name
%     varName = fieldNames{i};
%     % Use eval to assign the value to the variable with
%     eval([varName ' = w.' varName ';']);
% end




% %Get the policy functions
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Policy functions
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Hiring policies, looking at origins not allocations yet 
% [h_e_u, h_e_m, h_e_nm, h_e_t_m, h_e_t_nm, h_m_u, h_m_m, h_m_nm, h_m_t_m, h_m_t_nm, h_nm_u, h_nm_m, h_nm_nm, h_nm_t_m, h_nm_t_nm, h_t_u, h_t_m, h_t_nm, h_t_t_m, h_t_t_nm]...
%     =hire_policies(Vmh,Vnh,Veh,U,Vth,ats,tpts,0,cost_d,cost_p);
    
% % Allocation Policies after hiring at SM
% [p_e_m, p_e_n, p_m_m_u, p_m_m_d, p_m_n, p_n_n_u, p_n_n_p, p_n_m, p_t_m_u, p_t_m_d, p_t_n_u, p_t_n_p]...
%     =alloc_policies(Vmh,Vnh,U,Vth,ats,tpts,cost_d,cost_p);


% sum(e.eplus_tdist,"all")
% sum(e.eplus_ndist,"all")
% sum(e.eplus_mdist,"all")
% sum(e.eplus_edist,"all")
% sum(e.eplus_udist,"all")


