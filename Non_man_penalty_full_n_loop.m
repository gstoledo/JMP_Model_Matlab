%Hiring from non mager penalty matrix from allocation into non-manager
function [x]=Non_man_penalty_full_n_loop(Vth,tpts,ats,nm_penal)
    for i = 1:ats
        for z = 1:tpts
            x(i,z,:)=Vth(i,z,:);
            if nm_penal==1
                x(i,z, 2:end) = Vth(i,z,1:end-1); 
                x(i,z,1)=Vth(i,z,1);
            end
        end
    end         

    