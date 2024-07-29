%Hiring from non mager penalty matrix from allocation into manager
function [x]=Non_man_penalty_full_m_loop(Vth,tpts,ats,true)
    for i = 1:ats
        for q = 1:tpts
            x(i,:,q)=Vth(i,:,q);
            if true==1
                x(i, 2:end,q) = Vth(i, 1:end-1,q); 
                x(i,1,q)=Vth(i,1,q);
            end
        end
    end         

    