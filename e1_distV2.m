%Distributions after the SM phase
function [e1_udist,e1_edist,e1_mdist,e1_ndist,e1_tdist]=e1_distV2(Veh,Vmh,Vnh,Uh,Vth,ats,tpts,cost_d,cost_p,e_udist,e_edist,e_mdist,e_ndist,e_tdist,lamu,lam,n,del)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Policy functions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Hiring policies, looking at origins not allocations yet 
    [h_e_u, h_e_m, h_e_nm, h_e_t_m, h_e_t_nm, h_m_u, h_m_m, h_m_nm, h_m_t_m, h_m_t_nm, h_nm_u, h_nm_m, h_nm_nm, h_nm_t_m, h_nm_t_nm, h_t_u, h_t_m, h_t_nm, h_t_t_m, h_t_t_nm]...
    =hire_policies(Vmh,Vnh,Veh,Uh,Vth,ats,tpts,true,cost_d,cost_p);
    
    % Allocation Policies after hiring at SM
    [p_e_m, p_e_n, p_m_m_u, p_m_m_d, p_m_n, p_n_n_u, p_n_n_p, p_n_m, p_t_m_u, p_t_m_d, p_t_n_u, p_t_n_p]...
    =alloc_policies(Vmh,Vnh,Uh,Vth,ats,tpts,cost_d,cost_p);
     
 %%%%%%
 % Distributions and LOMs 
 u=sum(e_udist);
 % Distributions after the SM (implied by hiring and allocation policies)       
 % Unemployed distribution after SM
 e1_udist = zeros(1,tpts);
    for z = 1:tpts%Loop over types of a and q where z could have been 
        %%%%Replacement hiring
        %Firm with manager (a,z) replaces z with hires
        store_m_a=zeros(1,ats);
        for a = 1:ats
            store_u=zeros(1,tpts);
            store_m=zeros(ats,tpts);
            store_n=zeros(ats,tpts);
            store_t=zeros(ats,tpts,tpts);
            for z_tilda=1:tpts
                %From unemployed
                store_u(z_tilda)= e_udist(z_tilda)*h_m_u(z_tilda,a,z)*p_m_m_u(a,z,z_tilda);
                for a_tilda=1:ats
                    %From firm with manager
                    store_m(a_tilda,z_tilda)=e_mdist(a_tilda,z_tilda)*h_m_m(a_tilda,z_tilda,a,z)*p_m_m_u(a,z,z_tilda);
                    for q_tilda=1:tpts
                        %From firm with no manager
                        store_n(a_tilda,q_tilda)=e_ndist(a_tilda,q_tilda)*h_m_nm(a_tilda,q_tilda,a,z)*p_m_m_u(a,z,q_tilda);
                        %From team
                        store_t(a_tilda,z_tilda,q_tilda)=e_tdist(a_tilda,z_tilda,q_tilda)*(h_m_t_m(a_tilda,z_tilda,q_tilda,a,z)*p_m_m_u(a,z,z_tilda)+h_m_t_nm(a_tilda,z_tilda,q_tilda,a,z)*p_m_m_u(a,z,q_tilda));
                    end    
                end
            end
            m_RH_u=(lamu/u)*sum(store_u);
            m_RH_m=(lam/n)*sum(store_m,"all");
            m_RH_n=(lam/n)*sum(store_n,"all");
            m_RH_t=(lam/n)*sum(store_t,"all");
            
            store_m_a(a)=e_mdist(a,z)*(del +m_RH_u+ m_RH_m + m_RH_n + m_RH_t);
        end
    
    
                    
                    
                    
                    
        %Firm with team (a,z,q) replaces z with hires
        store_tm_aq=zeros(ats,tpts);
        for a=1:ats
            for q=1:tpts
                store_u=zeros(1,tpts);
                store_m=zeros(ats,tpts);
                store_n=zeros(ats,tpts);
                store_t=zeros(ats,tpts,tpts);
                for z_tilda=1:tpts
                    %From unemployed
                    store_u(z_tilda)= e_udist(z_tilda)*h_t_u(z_tilda,a,z,q)*(p_t_m_u(a,z,q,z_tilda)+p_t_n_p(a,z,q,z_tilda));
                    for a_tilda=1:ats
                        %From firm with manager
                        store_m(a_tilda,z_tilda)=e_mdist(a_tilda,z_tilda)*h_t_m(a_tilda,z_tilda,a,z,q)*(p_t_m_u(a,z,q,z_tilda)+p_t_n_p(a,z,q,z_tilda));
                        for q_tilda=1:tpts
                            %From firm with no manager
                            store_n(a_tilda,q_tilda)=e_ndist(a_tilda,q_tilda)*h_t_nm(a_tilda,q_tilda,a,z,q)*(p_t_m_u(a,z,q,z_tilda)+p_t_n_p(a,z,q,z_tilda));
                            %From team
                            store_t(a_tilda,z_tilda,q_tilda)=e_tdist(a_tilda,z_tilda,q_tilda)*(h_t_t_m(a_tilda,z_tilda,q_tilda,a,z,q)*(p_t_m_u(a,z,q,z_tilda)+p_t_n_p(a,z,q,z_tilda))+h_t_t_nm(a_tilda,z_tilda,q_tilda,a,z,q)*(p_t_m_u(a,z,q,q_tilda)+p_t_n_p(a,z,q,q_tilda)));
                        end
                    end
                end
                
                tm_RH_u=(lamu/u)*sum(store_u);
                tm_RH_m=(lam/n)*sum(store_m,"all");
                tm_RH_n=(lam/n)*sum(store_n,"all");
                tm_RH_t=(lam/n)*sum(store_t,"all");
                
                store_tm_aq(a,q)=e_tdist(a,z,q)*(del+ tm_RH_u + tm_RH_m + tm_RH_n + tm_RH_t);
            end
        end
                    
    
    
            
        %Firm with non-manager (a,z) replaces z with hires
        store_n_a=zeros(1,ats);   
        for a=1:ats
            store_u=zeros(1,tpts);
            store_m=zeros(ats,tpts);
            store_n=zeros(ats,tpts);
            store_t=zeros(ats,tpts,tpts);
            for z_tilda=1:tpts
                %From unemployed
                store_u(z_tilda)= e_udist(z_tilda)*h_nm_u(z_tilda,a,z)*p_n_n_u(a,z,z_tilda);
                for a_tilda=1:ats
                    %From firm with manager
                    store_m(a_tilda,z_tilda)=e_mdist(a_tilda,z_tilda)*h_nm_m(a_tilda,z_tilda,a,z)*p_n_n_u(a,z,z_tilda);
                    for q_tilda=1:tpts
                        %From firm with no manager
                        store_n(a_tilda,q_tilda)=e_ndist(a_tilda,q_tilda)*h_nm_nm(a_tilda,q_tilda,a,z)*p_n_n_u(a,z,q_tilda);
                        %From team
                        store_t(a_tilda,z_tilda,q_tilda)=e_tdist(a_tilda,z_tilda,q_tilda)*(h_nm_t_m(a_tilda,z_tilda,q_tilda,a,z)*p_n_n_u(a,z,z_tilda)+h_nm_t_nm(a_tilda,z_tilda,q_tilda,a,z)*p_n_n_u(a,z,q_tilda));
                    end    
                end
            end
            
            n_RH_u=(lamu/u)*sum(store_u);
            n_RH_m=(lam/n)*sum(store_m,"all");
            n_RH_n=(lam/n)*sum(store_n,"all");
            n_RH_t=(lam/n)*sum(store_t,"all");
            
            store_n_a(a)=e_ndist(a,z)*(del+ n_RH_u + n_RH_m + n_RH_n + n_RH_t);
        end
    
        
        %Team with non-manager (a,q,z) replaces z with hires
        store_tn_aq=zeros(ats,tpts);        
        for a=1:ats
            for q=1:tpts 
                store_u=zeros(1,tpts);
                store_m=zeros(ats,tpts);
                store_n=zeros(ats,tpts);
                store_t=zeros(ats,tpts,tpts);
                
                for z_tilda=1:tpts
                    %From unemployed
                    store_u(z_tilda)= e_udist(z_tilda)*h_t_u(z_tilda,a,q,z)*(p_t_n_u(a,q,z,z_tilda)+p_t_m_d(a,q,z,z_tilda));
                    for a_tilda=1:ats
                        %From firm with manager
                        store_m(a_tilda,z_tilda)=e_mdist(a_tilda,z_tilda)*h_t_m(a_tilda,z_tilda,a,q,z)*(p_t_n_u(a,q,z,z_tilda)+p_t_m_d(a,q,z,z_tilda));
                        for q_tilda=1:tpts
                            %From firm with no manager
                            store_n(a_tilda,q_tilda)=e_ndist(a_tilda,q_tilda)*h_t_nm(a_tilda,q_tilda,a,q,z)*(p_t_n_u(a,q,z,q_tilda)+p_t_m_d(a,q,z,q_tilda));
                            %From team
                            store_t(a_tilda,z_tilda,q_tilda)=e_tdist(a_tilda,z_tilda,q_tilda)*(h_t_t_m(a_tilda,z_tilda,q_tilda,a,q,z)*(p_t_n_u(a,q,z,z_tilda)+p_t_m_d(a,q,z,z_tilda))+h_t_t_nm(a_tilda,z_tilda,q_tilda,a,q,z)*(p_t_n_u(a,q,z,q_tilda)+p_t_m_d(a,q,z,q_tilda)));
                        end
                    end
                end
                
                tn_RH_u=(lamu/u)*sum(store_u);
                tn_RH_m=(lam/n)*sum(store_m,"all");
                tn_RH_n=(lam/n)*sum(store_n,"all");
                tn_RH_t=(lam/n)*sum(store_t,"all");
                
                
                store_tn_aq(a,q)=e_tdist(a,q,z)*(del+ tn_RH_u+ tn_RH_m + tn_RH_n + tn_RH_t);
            end
        end
        
            
        %Outflows from unemployment state z 
            store_hire_empty=zeros(1,ats);
            store_hire_m=zeros(ats,tpts);
            store_hire_n=zeros(ats,tpts);
            store_hire_team=zeros(ats,tpts,tpts);
            for a_tilda=1:ats
                %Empty firm hires from unemployed
                store_hire_empty(a_tilda)=e_edist(a_tilda)*h_e_u(z,a_tilda);
                for z_tilda = 1:tpts
                    %Firm with manager hires from unemployed
                    store_hire_m(a_tilda,z_tilda)=e_mdist(a_tilda,z_tilda)*h_m_u(z,a_tilda,z_tilda);
                    %Firm with no manager hires from unemployed
                    store_hire_n(a_tilda,z_tilda)=e_ndist(a_tilda,z_tilda)*h_nm_u(z,a_tilda,z_tilda);
                    for q_tilda=1:tpts
                        %Team hires from unemployed
                        store_hire_team(a_tilda,z_tilda,q_tilda)=e_tdist(a_tilda,z_tilda,q_tilda)*h_t_u(z,a_tilda,q_tilda,z_tilda);
                    end
                end
            end
            outflows=(lamu/u)*(sum(store_hire_empty)+sum(store_hire_m,"all")+sum(store_hire_n,"all")+sum(store_hire_team,"all"));
                        
                    
            e1_udist(z)= e_udist(z)*(1-outflows) + sum(store_m_a) + sum(store_tm_aq,"all") + sum(store_n_a) + sum(store_tn_aq,"all");
    end
     
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % For empty firms (a)
    e1_edist=zeros(1,ats);
    for a=1:ats
        %Probability of empty firm to hire anyone
        store_u=zeros(1,tpts);
        store_m=zeros(ats,tpts);    
        store_n=zeros(ats,tpts);
        store_t=zeros(ats,tpts,tpts);
        
        for z_tilda=1:tpts
            %From unemployed
            store_u(z_tilda)= e_udist(z_tilda)*h_e_u(z_tilda,a);
            for a_tilda=1:ats
                %From firm with manager
                store_m(a_tilda,z_tilda)=e_mdist(a_tilda,z_tilda)*h_e_m(a_tilda,z_tilda,a);
                for q_tilda=1:tpts
                    %From firm with no manager
                    store_n(a_tilda,q_tilda)=e_ndist(a_tilda,q_tilda)*h_e_nm(a_tilda,q_tilda,a);
                    %From team
                    store_t(a_tilda,z_tilda,q_tilda)=e_tdist(a_tilda,z_tilda,q_tilda)*(h_e_t_m(a_tilda,z_tilda,q_tilda,a)+h_e_t_nm(a_tilda,z_tilda,q_tilda,a));
                end    
            end
        end
        Hire_any=(lamu/u)*sum(store_u)+ (lam/n)*sum(store_m,"all")+ (lam/n)*sum(store_n,"all")+ (lam/n)*sum(store_t,"all");
        
        %Firms w/ manager (a,z) having their manager poached (for all z)
        store_a=zeros(1,ats);
        store_m=zeros(ats,tpts);
        store_n=zeros(ats,tpts);
        store_t=zeros(ats,tpts,tpts);
        Poach_m=zeros(1,tpts);
        for z=1:tpts
            for a_tilda=1:ats
                %Poached from empty firm
                store_a(a_tilda)=e_edist(a_tilda)*h_e_m(a,z,a_tilda);
                for z_hat=1:tpts
                    %Poached from firm with manager
                    store_m(a_tilda,z_hat)=e_mdist(a_tilda,z_hat)*h_m_m(a,z,a_tilda,z_hat);
                    %Poached from firm with no manager
                    store_n(a_tilda,z_hat)=e_ndist(a_tilda,z_hat)*h_nm_m(a,z,a_tilda,z_hat);
                    for q_tilda=1:tpts
                        %Poached from team
                        store_t(a_tilda,z_hat,q_tilda)=e_tdist(a_tilda,z_hat,q_tilda)*h_t_m(a,z,a_tilda,z_hat,q_tilda);
                    end
                end
            end
            %later have to add delta and multiply by e_m(a,z) for every z
            Poach_m(z)=(lam/n)*(sum(store_a)+sum(store_m,"all")+sum(store_n,"all")+sum(store_t,"all"));
        end
        
            
        %Firm w/ non manager (a,q) having their non manager poached (for all q)
        store_a=zeros(1,ats);
        store_m=zeros(ats,tpts);
        store_n=zeros(ats,tpts);
        store_t=zeros(ats,tpts,tpts);
        
        Poach_n=zeros(1,tpts);
        for q=1:tpts
            for a_tilda=1:ats
                %Poached from empty firm
                store_a(a_tilda)=e_edist(a_tilda)*h_e_nm(a,q,a_tilda);
                for z_hat=1:tpts
                    %Poached from firm with manager
                    store_m(a_tilda,z_hat)=e_mdist(a_tilda,z_hat)*h_m_nm(a,q,a_tilda,z_hat);
                    %Poached from firm with no manager
                    store_n(a_tilda,z_hat)=e_ndist(a_tilda,z_hat)*h_nm_nm(a,q,a_tilda,z_hat);
                    for q_tilda=1:tpts
                        %Poached from team
                        store_t(a_tilda,z_hat,q_tilda)=e_tdist(a_tilda,z_hat,q_tilda)*h_t_nm(a,q,a_tilda,z_hat,q_tilda);
                    end
                end
            end
            %later have to multiply by e_n(a,q) for every q
            Poach_n(q)=(lam/n)*(sum(store_a)+sum(store_m,"all")+sum(store_n,"all")+sum(store_t,"all"));
        end
        
        %Finally we have the distribution of empty firms (a)
        e1_edist(a)=e_edist(a)*(1-lamu-lam - Hire_any)...
            +sum(e_mdist(a,:).*(Poach_m+del))...
            +sum(e_ndist(a,:).*(Poach_n+del));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % Distribution of firm with manager (a,z)
    e1_mdist=zeros(ats,tpts);
    for a=1:ats
        for z=1:tpts
            %Fractions that remains on (a,z) after the SM
            %Need to compute several objects
            store_u=zeros(1,tpts);
            store_m=zeros(ats,tpts);
            store_n=zeros(ats,tpts);
            store_t=zeros(ats,tpts,tpts);
            %Prob of (a,z) to hire anyone
            for z_tilda=1:tpts
                %From unemployed
                store_u(z_tilda)= e_udist(z_tilda)*h_m_u(z_tilda,a,z);
                for a_tilda=1:ats
                    %From firm with manager
                    store_m(a_tilda,z_tilda)=e_mdist(a_tilda,z_tilda)*h_m_m(a_tilda,z_tilda,a,z);
                    for q_tilda=1:tpts
                        %From firm with no manager
                        store_n(a_tilda,q_tilda)=e_ndist(a_tilda,q_tilda)*h_m_nm(a_tilda,q_tilda,a,z);
                        %From team
                        store_t(a_tilda,z_tilda,q_tilda)=e_tdist(a_tilda,z_tilda,q_tilda)*(h_m_t_m(a_tilda,z_tilda,q_tilda,a,z)+h_m_t_nm(a_tilda,z_tilda,q_tilda,a,z));
                    end    
                end
            end
            Hire_any=(lamu/u)*sum(store_u)+ (lam/n)*sum(store_m,"all")+ (lam/n)*sum(store_n,"all") +(lam/n)*sum(store_t,"all");
            
            %Prob of (a,z) having its manager poached
            store_a=zeros(1,ats);
            store_m=zeros(ats,tpts);
            store_n=zeros(ats,tpts);
            store_t=zeros(ats,tpts,tpts);
            
            for a_tilda=1:ats
                %Poached from empty firm
                store_a(a_tilda)=e_edist(a_tilda)*h_e_m(a,z,a_tilda);
                for z_tilda=1:tpts
                    %Poached from firm with manager
                    store_m(a_tilda,z_tilda)=e_mdist(a_tilda,z_tilda)*h_m_m(a,z,a_tilda,z_tilda);
                    %Poached from firm with no manager
                    store_n(a_tilda,z_tilda)=e_ndist(a_tilda,z_tilda)*h_nm_m(a,z,a_tilda,z_tilda);
                    for q_tilda=1:tpts
                        %Poached from team
                        store_t(a_tilda,z_tilda,q_tilda)=e_tdist(a_tilda,z_tilda,q_tilda)*h_t_m(a,z,a_tilda,z_tilda,q_tilda);
                    end
                end
            end
            Poach_m=(lam/n)*(sum(store_a)+sum(store_m,"all")+sum(store_n,"all")+sum(store_t,"all"));
            
            %Firm (a,z,q) having its q poached
            Poach_q=zeros(1,tpts);
            for q=1:tpts
                store_a=zeros(1,ats);
                store_m=zeros(ats,tpts);
                store_n=zeros(ats,tpts);
                store_t=zeros(ats,tpts,tpts);
                for a_tilda=1:ats
                    %Poached from empty firm
                    store_a(a_tilda)=e_edist(a_tilda)*h_e_t_nm(a,z,q,a_tilda);
                    for z_tilda=1:tpts
                        %Poached from firm with manager
                        store_m(a_tilda,z_tilda)=e_mdist(a_tilda,z_tilda)*h_m_t_nm(a,z,q,a_tilda,z_tilda);
                        %Poached from firm with no manager
                        store_n(a_tilda,z_tilda)=e_ndist(a_tilda,z_tilda)*h_nm_t_nm(a,z,q,a_tilda,z_tilda);
                        for q_tilda=1:tpts
                            %Poached from team
                            store_t(a_tilda,z_tilda,q_tilda)=e_tdist(a_tilda,z_tilda,q_tilda)*h_t_t_nm(a,z,q,a_tilda,z_tilda,q_tilda);
                        end
                    end
                end
            %later have to add delta and multiply by e(a,z,q) for every q
            Poach_q(q)=(lam/n)*(sum(store_a)+sum(store_m,"all")+sum(store_n,"all")+sum(store_t,"all"));
            end
            
            
            % Inflows from Unemploment ending up in a,z
            store_ztilda=zeros(1,tpts);
            for z_tilda=1:tpts
                store_ztilda(z_tilda)=e_mdist(a,z_tilda)*h_m_u(z,a,z_tilda)*p_m_m_u(a,z_tilda,z);
            end
            inflow_u=(lamu/u)*(e_edist(a)*h_e_u(z,a)*p_e_m(a,z)+sum(store_ztilda));
            
            
            
            %InFlows from z's out there ending up in the right firm (ending up in a,z)
            store_atilda_z=zeros(1,ats); % From manager firms (a_tilde,z)
            for a_tilda=1:ats
                store_ztilda=zeros(1,tpts);
                for z_tilda=1:tpts
                    store_ztilda(z_tilda)=e_mdist(a,z_tilda)*h_m_m(a_tilda,z,a,z_tilda)*p_m_m_u(a,z_tilda,z);
                end
                store_atilda_z(a_tilda)=(lam/n)*(e_edist(a)*h_e_m(a_tilda,z,a)*p_e_m(a,z)+sum(store_ztilda));
            end
            
            store_atilda_z_qtilda=zeros(ats,tpts); % From team firms (a_tilde,z,q_tilde)
            for a_tilda=1:ats
                for q_tilda=1:tpts
                    store_ztilda=zeros(1,tpts);
                    for z_tilda=1:tpts
                        store_ztilda(z_tilda)=e_mdist(a,z_tilda)*h_m_t_m(a_tilda,z,q_tilda,a,z_tilda)*p_m_m_u(a,z_tilda,z);
                    end
                    store_atilda_z_qtilda(a_tilda,q_tilda)=(lam/n)*(e_edist(a)*h_e_t_m(a_tilda,z,q_tilda,a)*p_e_m(a,z)+sum(store_ztilda));
                end
            end
            
            store_atilda_z_nm=zeros(1,ats); % From non manager firms (a_tilde,z)
            for a_tilda=1:ats
                store_ztilda=zeros(1,tpts);
                for z_tilda=1:tpts
                    store_ztilda(z_tilda)=e_mdist(a,z_tilda)*h_m_nm(a_tilda,z,a,z_tilda)*p_m_m_u(a,z_tilda,z);
                end
                store_atilda_z_nm(a_tilda)=(lam/n)*(e_edist(a)*h_e_nm(a_tilda,z,a)*p_e_m(a,z)+sum(store_ztilda));
            end
            
            
            store_atilda_qtilda_z=zeros(ats,tpts); % From team firms (a_tilde,q_tilde,z)
            for a_tilda=1:ats
                for q_tilda=1:tpts
                    store_ztilda=zeros(1,tpts);
                    for z_tilda=1:tpts
                        store_ztilda(z_tilda)=e_mdist(a,z_tilda)*h_m_t_nm(a_tilda,z,q_tilda,a,z_tilda)*p_m_m_u(a,z_tilda,z);
                    end
                    store_atilda_qtilda_z(a_tilda,q_tilda)=(lam/n)*(e_edist(a)*h_e_t_nm(a_tilda,z,q_tilda,a)*p_e_m(a,z)+sum(store_ztilda));
                end
            end 
            
            
            %So finally we have the distribution of firms with manager (a,z)
            e1_mdist(a,z)=e_udist(z)*inflow_u...
            +e_mdist(a,z)*(1-del-lamu-lam-Hire_any-Poach_m)...
            +sum(reshape(e_tdist(a,z,:),1,tpts).*(Poach_q+del))...
            +sum(e_mdist(:,z)'.*store_atilda_z)...
            +sum(reshape(e_tdist(:,z,:),ats,tpts).*store_atilda_z_qtilda,"all")...
            +sum(e_ndist(:,z)'.*store_atilda_z_nm)...
            +sum(reshape(e_tdist(:,:,z),ats,tpts).*store_atilda_qtilda_z,"all");
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Now same for non-manager firms (a,q)

    e1_ndist=zeros(ats,tpts);
    for a=1:ats
        for q=1:tpts
            %Fractions that remains on (a,q) after the SM
            store_u=zeros(1,tpts);
            store_m=zeros(ats,tpts);
            store_n=zeros(ats,tpts);
            store_t=zeros(ats,tpts,tpts);        
            %Prob of (a,q) to hire anyone
            for z_tilda=1:tpts
                %From unemployed
                store_u(z_tilda)= e_udist(z_tilda)*h_nm_u(z_tilda,a,q);
                for a_tilda=1:ats
                    %From firm with manager
                    store_m(a_tilda,z_tilda)=e_mdist(a_tilda,z_tilda)*h_nm_m(a_tilda,z_tilda,a,q);
                    for q_tilda=1:tpts
                        %From firm with no manager
                        store_n(a_tilda,q_tilda)=e_ndist(a_tilda,q_tilda)*h_nm_nm(a_tilda,q_tilda,a,q);
                        %From team
                        store_t(a_tilda,z_tilda,q_tilda)=e_tdist(a_tilda,z_tilda,q_tilda)*(h_nm_t_m(a_tilda,z_tilda,q_tilda,a,q)+h_nm_t_nm(a_tilda,z_tilda,q_tilda,a,q));
                    end    
                end
            end
            Hire_any=(lamu/u)*sum(store_u)+ (lam/n)*sum(store_m,"all")+ (lam/n)*sum(store_n,"all")+ (lam/n)*sum(store_t,"all");
            
            
            %Prob of (a,q) having its non-manager poached
            store_a=zeros(1,ats);
            store_m=zeros(ats,tpts);
            store_n=zeros(ats,tpts);
            store_t=zeros(ats,tpts,tpts);
            
            for a_tilda=1:ats
                %Poached from empty firm
                store_a(a_tilda)=e_edist(a_tilda)*h_e_nm(a,q,a_tilda);
                for z_tilda=1:tpts
                    %Poached from firm with manager
                    store_m(a_tilda,z_tilda)=e_mdist(a_tilda,z_tilda)*h_m_nm(a,q,a_tilda,z_tilda);
                    %Poached from firm with no manager
                    store_n(a_tilda,z_tilda)=e_ndist(a_tilda,z_tilda)*h_nm_nm(a,q,a_tilda,z_tilda);
                    for q_tilda=1:tpts
                        %Poached from team
                        store_t(a_tilda,z_tilda,q_tilda)=e_tdist(a_tilda,z_tilda,q_tilda)*h_t_nm(a,q,a_tilda,z_tilda,q_tilda);
                    end
                end
            end
            Poach_n=(lam/n)*(sum(store_a)+sum(store_m,"all")+sum(store_n,"all")+sum(store_t,"all"));
            
            %Firm (a,z,q) having its manager z poached, for all z
            Poach_z=zeros(1,tpts);
            for z=1:tpts
                store_a=zeros(1,ats);
                store_m=zeros(ats,tpts);
                store_n=zeros(ats,tpts);
                store_t=zeros(ats,tpts,tpts);
                for a_tilda=1:ats
                    %Poached from empty firm
                    store_a(a_tilda)=e_edist(a_tilda)*h_e_t_m(a,z,q,a_tilda);
                    for z_tilda=1:tpts
                        %Poached from firm with manager
                        store_m(a_tilda,z_tilda)=e_mdist(a_tilda,z_tilda)*h_m_t_m(a,z,q,a_tilda,z_tilda);
                        %Poached from firm with no manager
                        store_n(a_tilda,z_tilda)=e_ndist(a_tilda,z_tilda)*h_nm_t_m(a,z,q,a_tilda,z_tilda);
                        for q_tilda=1:tpts
                            %Poached from team
                            store_t(a_tilda,z_tilda,q_tilda)=e_tdist(a_tilda,z_tilda,q_tilda)*h_t_t_m(a,z,q,a_tilda,z_tilda,q_tilda);
                        end
                    end
                end
                %later have to add delta and multiply by e(a,z,q) for every z 
                Poach_z(z)=(lam/n)*(sum(store_a)+sum(store_m,"all")+sum(store_n,"all")+sum(store_t,"all"));
            end 
            
            
            % Inflows from Unemploment ending up in a,q
            store_ztilda=zeros(1,tpts);
            for z_tilda=1:tpts
                store_ztilda(z_tilda)=e_ndist(a,z_tilda)*h_nm_u(q,a,z_tilda)*p_n_n_u(a,z_tilda,q);
            end
            inflow_u=(lamu/u)*(e_edist(a)*h_e_u(q,a)*p_e_n(a,q)+sum(store_ztilda));
            
            %InFlows from q's out there ending up in the right firm (ending up in a,q)
            store_atilda_q=zeros(1,ats); % From manager firms (a_tilde,q)
            for a_tilda=1:ats
                store_ztilda=zeros(1,tpts);
                for z_tilda=1:tpts
                    store_ztilda(z_tilda)=e_ndist(a,z_tilda)*h_nm_m(a_tilda,q,a,z_tilda)*p_n_n_u(a,z_tilda,q);
                end
                store_atilda_q(a_tilda)=(lam/n)*(e_edist(a)*h_e_m(a_tilda,q,a)*p_e_n(a,q)+sum(store_ztilda));
            end
            
            store_atilda_q_ztilda=zeros(ats,tpts); % From team firms (a_tilda,q,z_tilda)
            for a_tilda=1:ats
                for z_tilda=1:tpts
                    store_qtilda=zeros(1,tpts);
                    for q_tilda=1:tpts
                        store_qtilda(q_tilda)=e_ndist(a,q_tilda)*h_nm_t_m(a_tilda,q,z_tilda,a,q_tilda)*p_n_n_u(a,q_tilda,q);
                    end
                    store_atilda_q_ztilda(a_tilda,z_tilda)=(lam/n)*(e_edist(a)*h_e_t_m(a_tilda,q,z_tilda,a)*p_e_n(a,q)+sum(store_qtilda));
                end
            end
            
            store_atilda_q_nm=zeros(1,ats); % From non manager firms (a_tilde,q)
            for a_tilda=1:ats
                store_ztilda=zeros(1,tpts);
                for z_tilda=1:tpts
                    store_ztilda(z_tilda)=e_ndist(a,z_tilda)*h_nm_nm(a_tilda,q,a,z_tilda)*p_n_n_u(a,z_tilda,q);
                end
                store_atilda_q_nm(a_tilda)=(lam/n)*(e_edist(a)*h_e_nm(a_tilda,q,a)*p_e_n(a,q)+sum(store_ztilda));
            end
            
            store_atilda_ztilda_q=zeros(ats,tpts); % From team firms (a_tilde,z_tilda,q)
            for a_tilda=1:ats
                for z_tilda=1:tpts
                    store_qtilda=zeros(1,tpts);
                    for q_tilda=1:tpts
                        store_qtilda(q_tilda)=e_ndist(a,q_tilda)*h_nm_t_nm(a_tilda,z_tilda,q,a,q_tilda)*p_n_n_u(a,q_tilda,q);
                    end
                    store_atilda_ztilda_q(a_tilda,z_tilda)=(lam/n)*(e_edist(a)*h_e_t_nm(a_tilda,z_tilda,q,a)*p_e_n(a,q)+sum(store_qtilda));
                end
            end
            
            %So finally we have the distribution of firms with no manager (a,q) 
            e1_ndist(a,q)=e_udist(q)*inflow_u...
            +e_ndist(a,q)*(1-del-lamu-lam-Hire_any-Poach_n)...
            +sum(e_tdist(a,:,q).*(Poach_z+del))...
            +sum(e_mdist(:,q)'.*store_atilda_q)...
            +sum(reshape(e_tdist(:,q,:),ats,tpts).*store_atilda_q_ztilda,"all")...
            +sum(e_ndist(:,q)'.*store_atilda_q_nm)...
            +sum(reshape(e_tdist(:,:,q),ats,tpts).*store_atilda_ztilda_q,"all");
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% For firms with team (a,z,q)
    e1_tdist=zeros(ats,tpts,tpts);%Now the most complicated one and the logic is a bit different
    for a=1:ats
        for z=1:tpts
            for q=1:tpts
                %%% For z
                %Probability of becoming a team firm (a,z,q) hiring z from unemployed to become manager
                store_zm=zeros(1,tpts);
                store_zn=zeros(1,tpts);
                for z_tilda=1:tpts
                    store_zm(z_tilda)=e_tdist(a,z_tilda,q)*h_t_u(z,a,z_tilda,q)*p_t_m_u(a,z_tilda,q,z);
                    store_zn(z_tilda)=e_tdist(a,q,z_tilda)*h_t_u(z,a,q,z_tilda)*p_t_m_d(a,q,z_tilda,z);
                end
                Gamma_z_u_m=(lamu/u)*(e_mdist(a,q)*h_m_u(z,a,q)*p_m_m_d(a,q,z)+ e_ndist(a,q)*h_nm_u(z,a,q)*p_n_m(a,q,z)+sum(store_zm)+sum(store_zn));
                
                
                %Probability of becoming a team firm (a,z,q) hiring z from manager firm (a_tilda,z) to become manager
                store_atildaz=zeros(1,ats);
                for a_tilda=1:ats
                    store_zm=zeros(1,tpts);
                    store_zn=zeros(1,tpts);
                    for z_tilda=1:tpts
                        store_zm(z_tilda)=e_tdist(a,z_tilda,q)*h_t_m(a_tilda,z,a,z_tilda,q)*p_t_m_u(a,z_tilda,q,z);
                        store_zn(z_tilda)=e_tdist(a,q,z_tilda)*h_t_m(a_tilda,z,a,q,z_tilda)*p_t_m_d(a,q,z_tilda,z);
                    end
                    store_atildaz(a_tilda)=(lam/n)*(e_mdist(a,q)*h_m_m(a_tilda,z,a,q)*p_m_m_d(a,q,z)+ e_ndist(a,q)*h_nm_m(a_tilda,z,a,q)*p_n_m(a,q,z)+sum(store_zm)+sum(store_zn));
                end
                
                
                %Probability of becoming a team firm (a,z,q) hiring z from non-manager firm (a_tilda,z) to become manager
                store_atildaz_nm=zeros(1,ats);
                for a_tilda=1:ats
                    store_zm=zeros(1,tpts);
                    store_zn=zeros(1,tpts);
                    for z_tilda=1:tpts
                        store_zm(z_tilda)=e_tdist(a,z_tilda,q)*h_t_nm(a_tilda,z,a,z_tilda,q)*p_t_m_u(a,z_tilda,q,z);
                        store_zn(z_tilda)=e_tdist(a,q,z_tilda)*h_t_nm(a_tilda,z,a,q,z_tilda)*p_t_m_d(a,q,z_tilda,z);
                    end
                    store_atildaz_nm(a_tilda)=(lam/n)*(e_mdist(a,q)*h_m_nm(a_tilda,z,a,q)*p_m_m_d(a,q,z)+ e_ndist(a,q)*h_nm_nm(a_tilda,z,a,q)*p_n_m(a,q,z)+sum(store_zm)+sum(store_zn));
                end
                
                
                %Probability of becoming a team firm (a,z,q) hiring z from a team firm (a_tilda,z,q_tilda) to become manager
                store_atilda_z_qtilda=zeros(ats,tpts);
                for a_tilda=1:ats
                    for q_tilda=1:tpts
                        store_zm=zeros(1,tpts);
                        store_zn=zeros(1,tpts);
                        for z_tilda=1:tpts
                            store_zm(z_tilda)=e_tdist(a,z_tilda,q)*h_t_t_m(a_tilda,z,q_tilda,a,z_tilda,q)*p_t_m_u(a,z_tilda,q,z);
                            store_zn(z_tilda)=e_tdist(a,q,z_tilda)*h_t_t_m(a_tilda,z,q_tilda,a,q,z_tilda)*p_t_m_d(a,q,z_tilda,z);
                        end
                        store_atilda_z_qtilda(a_tilda,q_tilda)=(lam/n)*(e_mdist(a,q)*h_m_t_m(a_tilda,z,q_tilda,a,q)*p_m_m_d(a,q,z)+ e_ndist(a,q)*h_nm_t_m(a_tilda,z,q_tilda,a,q)*p_n_m(a,q,z)+sum(store_zm)+sum(store_zn));
                    end
                end
                
                
                %Probability of becoming a team firm (a,z,q) hiring z from a team firm (a_tilda,q_tilda,z) to become manager
                store_atilda_qtilda_z=zeros(ats,tpts);
                for a_tilda=1:ats
                    for q_tilda=1:tpts
                        store_zm=zeros(1,tpts);
                        store_zn=zeros(1,tpts);
                        for z_tilda=1:tpts
                            store_zm(z_tilda)=e_tdist(a,z_tilda,q)*h_t_t_nm(a_tilda,q_tilda,z,a,z_tilda,q)*p_t_m_u(a,z_tilda,q,z);
                            store_zn(z_tilda)=e_tdist(a,q,z_tilda)*h_t_t_nm(a_tilda,q_tilda,z,a,q,z_tilda)*p_t_m_d(a,q,z_tilda,z);
                        end
                        store_atilda_qtilda_z(a_tilda,q_tilda)=(lam/n)*(e_mdist(a,q)*h_m_t_nm(a_tilda,q_tilda,z,a,q)*p_m_m_d(a,q,z)+ e_ndist(a,q)*h_nm_t_nm(a_tilda,q_tilda,z,a,q)*p_n_m(a,q,z)+sum(store_zm)+sum(store_zn));
                    end
                end
                
                %%%For q
                
                %Probability of becoming a team firm (a,z,q) hiring q from unemployed to become non-manager
                store_zm=zeros(1,tpts);
                store_zn=zeros(1,tpts);
                for z_tilda=1:tpts
                    store_zm(z_tilda)=e_tdist(a,z_tilda,z)*h_t_u(q,a,z_tilda,z)*p_t_n_p(a,z_tilda,z,q);
                    store_zn(z_tilda)=e_tdist(a,z,z_tilda)*h_t_u(q,a,z,z_tilda)*p_t_n_u(a,z,z_tilda,q);
                end
                Gamma_q_u_n=(lamu/u)*(e_mdist(a,z)*h_m_u(q,a,z)*p_m_n(a,z,q)+ e_ndist(a,z)*h_nm_u(q,a,z)*p_n_n_p(a,z,q)+sum(store_zm)+sum(store_zn));
                
                %Probability of becoming a team firm (a,z,q) hiring q from manager firm (a_tilda,q) to become non-manager
                store_atildaq=zeros(1,ats);
                for a_tilda=1:ats
                    store_zm=zeros(1,tpts);
                    store_zn=zeros(1,tpts);
                    for z_tilda=1:tpts
                        store_zm(z_tilda)=e_tdist(a,z_tilda,z)*h_t_m(a_tilda,q,a,z_tilda,z)*p_t_n_p(a,z_tilda,z,q);
                        store_zn(z_tilda)=e_tdist(a,z,z_tilda)*h_t_m(a_tilda,q,a,z,z_tilda)*p_t_n_u(a,z,z_tilda,q);
                    end
                    store_atildaq(a_tilda)=(lam/n)*(e_mdist(a,z)*h_m_m(a_tilda,q,a,z)*p_m_n(a,z,q)+ e_ndist(a,z)*h_nm_m(a_tilda,q,a,z)*p_n_n_p(a,z,q)+sum(store_zm)+sum(store_zn));
                end
                
                %Probability of becoming a team firm (a,z,q) hiring q from non-manager firm (a_tilda,q) to become non-manager
                store_atildaq_nm=zeros(1,ats);
                for a_tilda=1:ats
                    store_zm=zeros(1,tpts);
                    store_zn=zeros(1,tpts);
                    for z_tilda=1:tpts
                        store_zm(z_tilda)=e_tdist(a,z_tilda,z)*h_t_nm(a_tilda,q,a,z_tilda,z)*p_t_n_p(a,z_tilda,z,q);
                        store_zn(z_tilda)=e_tdist(a,z,z_tilda)*h_t_nm(a_tilda,q,a,z,z_tilda)*p_t_n_u(a,z,z_tilda,q);
                    end
                    store_atildaq_nm(a_tilda)=(lam/n)*(e_mdist(a,z)*h_m_nm(a_tilda,q,a,z)*p_m_n(a,z,q)+ e_ndist(a,z)*h_nm_nm(a_tilda,q,a,z)*p_n_n_p(a,z,q)+sum(store_zm)+sum(store_zn));
                end 
                
                %Probability of becoming a team firm (a,z,q) hiring q from a team firm (a,q,q_tilda) to become non-manager
                store_atildaq_qtilda=zeros(ats,tpts);
                for a_tilda=1:ats
                    for q_tilda=1:tpts
                        store_zm=zeros(1,tpts);
                        store_zn=zeros(1,tpts);
                        for z_tilda=1:tpts
                            store_zm(z_tilda)=e_tdist(a,z_tilda,z)*h_t_t_m(a_tilda,q,q_tilda,a,z_tilda,z)*p_t_n_p(a,z_tilda,z,q);
                            store_zn(z_tilda)=e_tdist(a,z,z_tilda)*h_t_t_m(a_tilda,q,q_tilda,a,z,z_tilda)*p_t_n_u(a,z,z_tilda,q);
                        end
                        store_atildaq_qtilda(a_tilda,q_tilda)=(lam/n)*(e_mdist(a,z)*h_m_t_m(a_tilda,q,q_tilda,a,z)*p_m_n(a,z,q)+ e_ndist(a,z)*h_nm_t_m(a_tilda,q,q_tilda,a,z)*p_n_n_p(a,z,q)+sum(store_zm)+sum(store_zn));
                    end
                end
                
                %Probability of becoming a team firm (a,z,q) hiring q from a team firm (a,q_tilda,q) to become non-manager
                store_atilda_qtilda_q=zeros(ats,tpts);
                for a_tilda=1:ats
                    for q_tilda=1:tpts
                        store_zm=zeros(1,tpts);
                        store_zn=zeros(1,tpts);
                        for z_tilda=1:tpts
                            store_zm(z_tilda)=e_tdist(a,z_tilda,z)*h_t_t_nm(a_tilda,q_tilda,q,a,z_tilda,z)*p_t_n_p(a,z_tilda,z,q);
                            store_zn(z_tilda)=e_tdist(a,z,z_tilda)*h_t_t_nm(a_tilda,q_tilda,q,a,z,z_tilda)*p_t_n_u(a,z,z_tilda,q);
                        end
                        store_atilda_qtilda_q(a_tilda,q_tilda)=(lam/n)*(e_mdist(a,z)*h_m_t_nm(a_tilda,q_tilda,q,a,z)*p_m_n(a,z,q)+ e_ndist(a,z)*h_nm_t_nm(a_tilda,q_tilda,q,a,z)*p_n_n_p(a,z,q)+sum(store_zm)+sum(store_zn));
                    end
                end
                
                %Mass that stay in the team firm (a,z,q)
                %Probability of (a,z,q) to hire anyone
                store_u=zeros(1,tpts);
                store_m=zeros(ats,tpts);
                store_n=zeros(ats,tpts);
                store_t=zeros(ats,tpts,tpts);
                for z_tilda=1:tpts
                    %From unemployed
                    store_u(z_tilda)= e_udist(z_tilda)*h_t_u(z_tilda,a,z,q);
                    for a_tilda=1:ats
                        %From firm with manager
                        store_m(a_tilda,z_tilda)=e_mdist(a_tilda,z_tilda)*h_t_m(a_tilda,z_tilda,a,z,q);
                        for q_tilda=1:tpts
                            %From firm with no manager
                            store_n(a_tilda,q_tilda)=e_ndist(a_tilda,q_tilda)*h_t_nm(a_tilda,q_tilda,a,z,q);
                            %From team
                            store_t(a_tilda,z_tilda,q_tilda)=e_tdist(a_tilda,z_tilda,q_tilda)*(h_t_t_m(a_tilda,z_tilda,q_tilda,a,z,q)+h_t_t_nm(a_tilda,z_tilda,q_tilda,a,z,q));
                        end
                    end
                end
                Hire_any=(lamu/u)*sum(store_u)+ (lam/n)*sum(store_m,"all")+ (lam/n)*sum(store_n,"all")+ (lam/n)*sum(store_t,"all");
                
                %Prob of (a,z,q) having any of its members poached
                store_a=zeros(1,ats);
                store_m=zeros(ats,tpts);
                store_n=zeros(ats,tpts);
                store_t=zeros(ats,tpts,tpts);
                for a_tilda=1:ats
                    %Poached from empty firm
                    store_a(a_tilda)=e_edist(a_tilda)*(h_e_t_m(a,z,q,a_tilda)+h_e_t_nm(a,z,q,a_tilda));
                    for z_tilda=1:tpts
                        %Poached from firm with manager
                        store_m(a_tilda,z_tilda)=e_mdist(a_tilda,z_tilda)*(h_m_t_m(a,z,q,a_tilda,z_tilda)+h_m_t_nm(a,z,q,a_tilda,z_tilda));
                        %Poached from firm with no manager
                        store_n(a_tilda,z_tilda)=e_ndist(a_tilda,z_tilda)*(h_nm_t_m(a,z,q,a_tilda,z_tilda)+h_nm_t_nm(a,z,q,a_tilda,z_tilda));
                        for q_tilda=1:tpts
                            %Poached from team
                            store_t(a_tilda,z_tilda,q_tilda)=e_tdist(a_tilda,z_tilda,q_tilda)*(h_t_t_m(a,z,q,a_tilda,z_tilda,q_tilda)+h_t_t_nm(a,z,q,a_tilda,z_tilda,q_tilda));
                        end
                    end
                end
                Poach=(lam/n)*(sum(store_a)+sum(store_m,"all")+sum(store_n,"all")+sum(store_t,"all"));
                
                
                %Finally we have the distribution of team firms (a,z,q)
                e1_tdist(a,z,q)=e_udist(z)*Gamma_z_u_m + e_udist(q)*Gamma_q_u_n...
                    +e_tdist(a,z,q)*(1-2*del-lamu-lam-Hire_any-Poach)...
                    +sum(e_mdist(:,z)'.*store_atildaz)...
                    +sum(e_ndist(:,z)'.*store_atildaz_nm)...
                    +sum(reshape(e_tdist(:,z,:),ats,tpts).*store_atilda_z_qtilda,"all")...
                    +sum(reshape(e_tdist(:,:,z),ats,tpts).*store_atilda_qtilda_z,"all")...
                    +sum(e_mdist(:,q)'.*store_atildaq)...
                    +sum(e_ndist(:,q)'.*store_atildaq_nm)...
                    +sum(reshape(e_tdist(:,q,:),ats,tpts).*store_atildaq_qtilda,"all")...
                    +sum(reshape(e_tdist(:,:,q),ats,tpts).*store_atilda_qtilda_q,"all");
            end
        end
    end
    



    



    
    
    
    

    