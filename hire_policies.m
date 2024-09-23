%%Hiring policies, looking at origins not allocations yet 
function [h_e_u, h_e_m, h_e_nm, h_e_t_m, h_e_t_nm, h_m_u, h_m_m, h_m_nm, h_m_t_m, h_m_t_nm, h_nm_u, h_nm_m, h_nm_nm, h_nm_t_m, h_nm_t_nm, h_t_u, h_t_m, h_t_nm, h_t_t_m, h_t_t_nm]=hire_policies(Vmh,Vnh,Veh,Uh,Vth,ats,tpts,true,cost_d,cost_p);
    % Gains from trade after convergence
    [gt_ef_meet_u, gt_ef_meet_m, gt_ef_meet_nm, gt_ef_meet_t_m, gt_ef_meet_t_nm, gt_ef_meet_t]=ef_gt(Vmh,Vnh,Veh,Uh,Vth,ats,tpts,true);
    [gt_em_meet_u, gt_em_meet_m, gt_em_meet_nm, gt_em_meet_t_m, gt_em_meet_t_nm, gt_em_meet_t]=mf_gt(Vmh,Vnh,Veh,Uh,Vth,ats,tpts,true,cost_d);
    [gt_en_meet_u, gt_en_meet_m, gt_en_meet_nm, gt_en_meet_t_m, gt_en_meet_t_nm, gt_en_meet_t]=nf_gt(Vmh,Vnh,Veh,Uh,Vth,ats,tpts,true,cost_p);
    [gt_tf_meet_u, gt_tf_meet_m, gt_tf_meet_nm, gt_tf_meet_t_m, gt_tf_meet_t_nm, gt_tf_meet_t]=tf_gt(Vmh,Vnh,Veh,Uh,Vth,ats,tpts,true,cost_p,cost_d);
    tol=1e-8;
    
    %Empty firm hiring policy
    %Hire policy by empty firms from unemployed
    h_e_u=double(gt_ef_meet_u>tol);
    
    %Hire policy by empty firms from firms with manager
    h_e_m=double(gt_ef_meet_m>tol);
    
    %Hire policy by empty firms from firms with no manager
    h_e_nm=double(gt_ef_meet_nm>tol);
    
    %Hire policy by empty firms from teams, hiring their manager 
    h_e_t_m=double(gt_ef_meet_t_m>=gt_ef_meet_t_nm & gt_ef_meet_t_m>tol);
    
    %Hire policy by empty firms from teams, hiring their non manager
    h_e_t_nm=double(gt_ef_meet_t_nm>gt_ef_meet_t_m & gt_ef_meet_t_nm>tol);
    
    %Firm with manager hiring policy
    %Hire policy by firms with manager from unemployed
    h_m_u=double(gt_em_meet_u>tol);
    
    %Hire policy by firms with manager from firms with manager
    h_m_m=double(gt_em_meet_m>tol);
    
    %Hire policy by firms with manager from firms with no manager
    h_m_nm=double(gt_em_meet_nm>tol);
    
    %Hire policy by firms with manager from teams, hiring their manager
    h_m_t_m=double(gt_em_meet_t_m>=gt_em_meet_t_nm & gt_em_meet_t_m>tol);
    
    %Hire policy by firms with manager from teams, hiring their non manager
    h_m_t_nm=double(gt_em_meet_t_nm>gt_em_meet_t_m & gt_em_meet_t_nm>tol);
    
    %Firm with no manager hiring policy
    %Hire policy by firms with no manager from unemployed
    h_nm_u=double(gt_en_meet_u>tol);
    
    %Hire policy by firms with no manager from firms with manager       
    h_nm_m=double(gt_en_meet_m>tol);
    
    %Hire policy by firms with no manager from firms with no manager
    h_nm_nm=double(gt_en_meet_nm>tol);
    
    %Hire policy by firms with no manager from teams, hiring their manager
    h_nm_t_m=double(gt_en_meet_t_m>=gt_en_meet_t_nm & gt_en_meet_t_m>tol);
    
    %Hire policy by firms with no manager from teams, hiring their non manager
    h_nm_t_nm=double(gt_en_meet_t_nm>gt_en_meet_t_m & gt_en_meet_t_nm>tol);
    
    %Team hiring policy
    %Hire policy by teams from unemployed
    h_t_u=double(gt_tf_meet_u>tol);
    
    %Hire policy by teams from firms with manager
    h_t_m=double(gt_tf_meet_m>tol);
    
    %Hire policy by teams from firms with no manager    
    h_t_nm=double(gt_tf_meet_nm>tol);
    
    %Hire policy by teams from teams, hiring their manager
    h_t_t_m=double(gt_tf_meet_t_m>=gt_tf_meet_t_nm & gt_tf_meet_t_m>tol);
    
    %Hire policy by teams from teams, hiring their non manager
    h_t_t_nm=double(gt_tf_meet_t_nm>gt_tf_meet_t_m & gt_tf_meet_t_nm>tol);
    
    
