%%Alloc policies, conditional on hiring
function [p_e_m, p_e_n, p_m_m_u, p_m_m_d, p_m_n, p_n_n_u, p_n_n_p, p_n_m, p_t_m_u, p_t_m_d, p_t_n_u, p_t_n_p]=alloc_policies(Vmh,Vnh,Uh,Vth,ats,tpts,cost_d,cost_p)
    % Allocation Policies after hiring at SM
    %Empty firm allocation policy
    %Probability of empty (a) firm allocating the hire z_tilda as manager 
    %p_e_m(a,z_tilda)
    p_e_m=double(Vmh>=Vnh);
    %Probability of empty (a) firm allocating the hire z_tilda as non manager
    %p_e_nm(a,z_tilda)
    p_e_n=double(Vmh<Vnh);
    
    %Firm with manager allocation policy 
    %Probability of firm with manager (a,z) allocating the hire z_tilda as manager, firing current manager
    %p_m_m_u(a,z,z_tilda)
    p_m_m_u=zeros(ats,tpts,tpts);
    %Probability of firm with manager (a,z) allocating the hire z_tilda as manager, demoting current manager
    %p_m_m_m(a,z,z_tilda)
    p_m_m_d=zeros(ats,tpts,tpts);
    
    %Probability of firm with manager (a,z) allocating the hire z_tilda as non manager
    %p_m_n(a,z,z_tilda)
    p_m_n=zeros(ats,tpts,tpts);
    
    
    %One loop for all of them 
    for a=1:ats
        for z=1:tpts
            for z_tilda=1:tpts
                p_m_m_u(a,z,z_tilda)= double(Vmh(a,z_tilda)+Uh(z)>max(Vth(a,z,z_tilda),Vth(a,z_tilda,z)-cost_d));
                p_m_m_d(a,z,z_tilda)= double(Vth(a,z_tilda,z)-cost_d>max(Vmh(a,z_tilda)+Uh(z),Vth(a,z,z_tilda)));
                p_m_n(a,z,z_tilda)= 1-p_m_m_u(a,z,z_tilda)-p_m_m_d(a,z,z_tilda);
            end
        end
    end
    
    
    %Firm with non-manager allocation policy
    %Probability of firm with non manager (a,q) allocating the hire z_tilda as non-manager, firing current non manager
    %p_n_n_u(a,q,z_tilda)
    p_n_n_u=zeros(ats,tpts,tpts);
    
    %Probability of firm with non manager (a,q) allocating the hire z_tilda as non-manager, promoting current non manager 
    %p_n_n_p(a,q,z_tilda)
    p_n_n_p=zeros(ats,tpts,tpts);
    
    %Probability of firm with non manager (a,q) allocating the hire z_tilda as manager
    %p_n_m(a,q,z_tilda) 
    p_n_m=zeros(ats,tpts,tpts);
    
    %One loop for all of them
    for a=1:ats
        for q=1:tpts
            for z_tilda=1:tpts
                p_n_n_u(a,q,z_tilda)= double(Vnh(a,z_tilda)+Uh(q)>max(Vth(a,q,z_tilda)-cost_p,Vth(a,z_tilda,q)));
                p_n_n_p(a,q,z_tilda)= double(Vth(a,q,z_tilda)-cost_p>max(Vnh(a,z_tilda)+Uh(q),Vth(a,z_tilda,q)));
                p_n_m(a,q,z_tilda)= 1-p_n_n_u(a,q,z_tilda)-p_n_n_p(a,q,z_tilda);
            end
        end
    end
    
    %Team allocation policy
    %Probability of team (a,z,q) allocating the hire z_tilda as manager, firing current manager
    %p_t_m_u(a,z,q,z_tilda)
    p_t_m_u=zeros(ats,tpts,tpts,tpts);
    
    %Probability of team (a,z,q) allocating the hire z_tilda as manager, demoting current manager
    %p_t_m_d(a,z,q,z_tilda)
    p_t_m_d=zeros(ats,tpts,tpts,tpts);
    
    %Probability of team (a,z,q) allocating the hire z_tilda as non manager, firing current non manager
    %p_t_n_u(a,z,q,z_tilda)
    p_t_n_u=zeros(ats,tpts,tpts,tpts);
    
    %Probability of team (a,z,q) allocating the hire z_tilda as non manager, promoting current non manager
    %p_t_n_p(a,z,q,z_tilda)
    p_t_n_p=zeros(ats,tpts,tpts,tpts);
    
    %One loop for all of them
    for z_tilda=1:tpts
        for z=1:tpts
            for q=1:tpts
                for a=1:ats
                    alloc =[ Vth(a,z_tilda,q)+Uh(z), Vth(a,z_tilda,z)-cost_d+Uh(q), Vth(a,z,z_tilda)+Uh(q),Vth(a,q,z_tilda)-cost_p + Uh(z)];
                    %p_t_m_u(a,z,q,z_tilda)= double(alloc(1)>=max(alloc([2 3 4])));
                    p_t_m_d(a,z,q,z_tilda)= double(alloc(2)>max(alloc([1 3 4])));
                    p_t_n_u(a,z,q,z_tilda)= double(alloc(3)>max(alloc([1 2 4])));
                    p_t_n_p(a,z,q,z_tilda)= double(alloc(4)>max(alloc([1 2 3])));
                    p_t_m_u(a,z,q,z_tilda)= 1-p_t_m_d(a,z,q,z_tilda)-p_t_n_u(a,z,q,z_tilda)-p_t_n_p(a,z,q,z_tilda);
                end
            end
        end
    end
    