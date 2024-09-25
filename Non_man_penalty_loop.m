function [x]=Non_man_penalty_loop(Vmh,tpts,ats,nm_penal)
    for i = 1:ats
    x(i,:)=Vmh(i,:);
    if nm_penal==1
        x(i, 2:end) = Vmh(i, 1:end-1); % Shift values to the right
        x(i,1)=Vmh(i,1);
    end
end
end