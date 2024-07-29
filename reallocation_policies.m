
%Reallocation Policies, before producing (using Vs not Vh)
function [d_m, r_m, d_n, r_n, d_t_m, d_t_n, r_t_m, r_t_n, d_t_b]=reallocation_policies(Ve,Vm,Vn,Vt,U,ats,tpts,cost_d,cost_p)
    %Notation here follows the paper, a bit wierd with previous notation in the code

    %For firm with manager (a,z)
    %Fire the manager
    %dm(a,z)
    d_m=zeros(ats,tpts);
    
    %Reallocate the manager
    r_m=zeros(ats,tpts);
    
    for z=1:tpts
        for a=1:ats
            realloc=[Ve(a)+U(z),Vn(a,z)-cost_p,Vm(a,z)];
            d_m(a,z)=double(realloc(1)>max(realloc([2,3])));
            r_m(a,z)=double(realloc(2)>max(realloc([1,3])));
        end
    end
    
    %For firm with no manager (a,q)
    %Fire the non manager
    %dn(a,q)
    d_n=zeros(ats,tpts);
    
    %Reallocate the non manager
    r_n=zeros(ats,tpts);
    for q=1:tpts
        for a=1:ats
            realloc=[Ve(a)+U(q),Vm(a,q)-cost_p,Vn(a,q)];
            d_n(a,q)=double(realloc(1)>max(realloc(2:3)));
            r_n(a,q)=double(realloc(2)>max(realloc([1,3])));
        end
    end
    
    %For firm with team (a,z,q)
    %We need 5 objects to decrtive final alloc and we will use combinations of them
    %Fire the manager only
    %d_t_m(a,z,q)
    d_t_m=zeros(ats,tpts,tpts);
    
    %Fire the non manager only
    %d_t_n(a,z,q)
    d_t_n=zeros(ats,tpts,tpts);
    
    %Reallocate the manager
    %r_t_m(a,z,q)
    r_t_m=zeros(ats,tpts,tpts);
    
    %Reallocate the non manager
    %r_t_n(a,z,q)
    r_t_n=zeros(ats,tpts,tpts);
    
    %Fire both
    %d_t_b(a,z,q)
    d_t_b=zeros(ats,tpts,tpts);
    
    %Page 25 of the pdf helps
    for q=1:tpts
        for z=1:tpts
            for a=1:ats
                realloc=[Vt(a,z,q),...%stay idle 1
                Vn(a,q)-cost_p+U(z),...%fire manager and promote 2
                Vn(a,z)-cost_d+U(q),...%fire non manager and demote 3
                Vm(a,z)+U(q),...%Fire non manager 4
                Vn(a,q)+U(z),... %Fire manager 5
                Vt(a,q,z),... %Reallocte both 6
                Ve(a)+U(q)+U(z)]; %Fire both 7
                d_t_m(a,z,q)=double(max(realloc([2,5]))>max(realloc([1,3,4,6,7])));
                d_t_n(a,z,q)=double(max(realloc([3,4]))>max(realloc([1,2,5,6,7])));
                r_t_m(a,z,q)=double(max(realloc([3,6]))>max(realloc([1,2,4,5,7])));
                r_t_n(a,z,q)=double(max(realloc([2,6]))>max(realloc([1,3,4,5,7])));
                d_t_b(a,z,q)=double(realloc(7)>max(realloc([1,2,3,4,5,6])));
            end
        end
    end
    