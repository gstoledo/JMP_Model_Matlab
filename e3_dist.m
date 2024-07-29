function [e3_udist,e3_edist,e3_mdist,e3_ndist,e3_tdist]=e3_dist(ats,tpts,e2_udist,e2_edist,e2_mdist,e2_ndist,e2_tdist,u_trans,a_trans,q_trans)
    %Disrtribtuions after shocks 
    %e3 distributions
    
    %Dist of unemployed
    e3_udist=sum(u_trans.*e2_udist',1); 
    
    %Dist of empty firms
    e3_edist=sum(a_trans.*e2_edist',1);
    
    %Dist of firms with manager
    e3_mdist=zeros(ats,tpts);
    
    for z=1:tpts
        e3_mdist(:,z)=sum(a_trans.*e2_mdist(:,z),1)';
    end
    
    
    %Dist of firms with no manager
    e3_ndist=zeros(ats,tpts);
    for a=1:ats
        for q=1:tpts
            store_aprime_qprime=zeros(ats,tpts);
            for a_prime=1:ats
                for q_prime=1:tpts
                    store_aprime_qprime(a_prime,q_prime)=a_trans(a_prime,a)*q_trans(q_prime,q)*e2_ndist(a_prime,q_prime);
                end
            end
            e3_ndist(a,q)=sum(store_aprime_qprime,"all");
        end
    end
    
    %Dist of firms with team
    e3_tdist=zeros(ats,tpts,tpts);
    for a=1:ats
        for z=1:tpts
            for q=1:tpts
                store_aprime_z_qprime=zeros(ats,tpts);
                for a_prime=1:ats
                    for q_prime=1:tpts
                        store_aprime_z_qprime(a_prime,q_prime)=a_trans(a_prime,a)*q_trans(q_prime,q)*e2_tdist(a_prime,z,q_prime);
                    end
                end
                e3_tdist(a,z,q)=sum(store_aprime_z_qprime,"all");
            end
        end
    end