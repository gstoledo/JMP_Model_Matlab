function [e2_udist,e2_edist,e2_mdist,e2_ndist,e2_tdist]=e2_dist(Ve,Vm,Vn,U,Vt,ats,tpts,cost_d,cost_p,e1_udist,e1_edist,e1_mdist,e1_ndist,e1_tdist)
    %Reallocation Policies, before producing (using Vs not Vh)
    % Notation here follows the paper, a bit wierd with previous notation in the code
    [d_m, r_m, d_n, r_n, d_t_m, d_t_n, r_t_m, r_t_n, d_t_b]...
    =reallocation_policies(Ve,Vm,Vn,Vt,U,ats,tpts,cost_d,cost_p);
    
    
    %Unemployed distribution
    e2_udist=zeros(1,tpts);
    %How to become unemployed z at this stage
    for z=1:tpts
        store_a=zeros(1,ats);
        store_qm=zeros(ats,tpts);
        store_qn=zeros(ats,tpts);
        for a=1:ats
            % Firm with manager or no manager z firing its worker
            store_a(a)=e1_mdist(a,z)*d_m(a,z)+e1_ndist(a,z)*d_n(a,z);
            for q=1:tpts
                % Team firing its manager z
                store_qn(a,q)=e1_tdist(a,z,q)*(d_t_m(a,z,q)+d_t_b(a,z,q));
                % Team firing its non manager z
                store_qm(a,q)=e1_tdist(a,q,z)*(d_t_n(a,q,z)+d_t_b(a,q,z));
            end
        end
        e2_udist(z)=e1_udist(z)+sum(store_a)+sum(store_qm,"all")+sum(store_qn,"all");
    end
    
    
    %Empty firm distribution
    e2_edist=zeros(1,ats);
    for a=1:ats
        %Firm With manager or non-manager firing their only worker
        e2_edist(a)=e1_edist(a) + e1_mdist(a,:)*d_m(a,:)' + e1_ndist(a,:)*d_n(a,:)' + sum(e1_tdist(a,:,:).*d_t_b(a,:,:),"all");
    end
    
    
    %Firm with manager distribution
    e2_mdist=zeros(ats,tpts);
    
    for a=1:ats
        for z=1:tpts
            store_q=zeros(1,tpts);
            for q=1:tpts
                %Team firing its non manager q or reallocating z into manager
                store_q(q)=e1_tdist(a,z,q)*d_t_n(a,z,q)*(1-r_t_m(a,z,q)) + e1_tdist(a,q,z)*d_t_m(a,q,z)*r_t_n(a,q,z);
            end
            e2_mdist(a,z)=e1_mdist(a,z)*(1-d_m(a,z)-r_m(a,z)) + e1_ndist(a,z)*r_n(a,z) + sum(store_q);
        end
    end
    
    
    %Firm with no manager distribution
    e2_ndist=zeros(ats,tpts);
    
    for a=1:ats
        for q=1:tpts
            store_z=zeros(1,tpts);
            for z=1:tpts
                %Team firing its manager z or reallocating q into non manager
                store_z(z)=e1_tdist(a,z,q)*d_t_m(a,z,q)*(1-r_t_n(a,z,q)) + e1_tdist(a,q,z)*d_t_n(a,q,z)*r_t_m(a,q,z);
            end
            e2_ndist(a,q)=e1_ndist(a,q)*(1-d_n(a,q)-r_n(a,q)) + e1_mdist(a,q)*r_m(a,q) + sum(store_z);
        end 
    end
    
    %Firm with team distribution
    e2_tdist=zeros(ats,tpts,tpts);
    
    for a=1:ats
        for z=1:tpts
            for q=1:tpts
                e2_tdist(a,z,q)=e1_tdist(a,z,q)*(1-d_t_m(a,z,q)-d_t_n(a,z,q)-d_t_b(a,z,q)-r_t_m(a,z,q)*r_t_n(a,z,q)) + e1_tdist(a,q,z)*r_t_m(a,q,z)*r_t_n(a,q,z);
            end
        end
    end
    
    
    
