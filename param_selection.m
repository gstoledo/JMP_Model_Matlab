function [p_vec, ps] = param_selection(theta, selection_param)
    % Initialize an empty array for selected parameters
    p_vec = [];
    
    % Loop over the selection_param and extract the values
    for i = 1:length(selection_param)
        field_name = selection_param{i};
        
        % Check if the field contains a scalar or a vector
        field_value = theta.(field_name);
        
        if isvector(field_value)
            % Concatenate the entire vector into p_vec
            p_vec = [p_vec, field_value(:)'];  % Ensure row vector
        else
            % Concatenate scalar values
            p_vec = [p_vec, field_value];
        end
    end
    
    % Remove the selected parameters from the structure
    ps = rmfield(theta, selection_param);
end