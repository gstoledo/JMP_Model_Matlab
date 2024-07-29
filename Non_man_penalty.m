
function [x]=Non_man_penalty(Vmh,tpts,ats,true)
    for i = 1:ats
    % Store the first row of the current layer
    x(:,:,i)=Vmh(:,:,i);
    if true==1
        % Assuming Vmh(:,:,i) is a 2D matrix and you want to shift values to the right in each row
        x(:, 2:end, i) = Vmh(:, 1:end-1, i); % Shift values to the right
        x(:,1,i)=Vmh(:,1, i);
    end
end
end