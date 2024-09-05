function p=p_vec_to_struct(p_vec, p_selection, theta)
    % Initialize an empty struct to hold the parameters
    p = struct();
    % Keep track of the index in paramVec
    paramIdx = 1;

    % Loop through the parameters in p_selection
    for i = 1:numel(p_selection)
        param_name = p_selection{i};

        % Check if the parameter exists in theta and if it is a vector
        if isfield(theta, param_name) && isvector(theta.(param_name))
            % Get the size of the vector from the original theta parameter
            vector_size = numel(theta.(param_name));

            % Allocate the corresponding elements from paramVec as a column vector
            p.(param_name) = p_vec(paramIdx:paramIdx + vector_size - 1)';

            % Update the paramIdx to skip the allocated vector elements
            paramIdx = paramIdx + vector_size;
        else
            % If it's a scalar parameter, just assign the current value from paramVec
            p.(param_name) = p_vec(paramIdx);

            % Increment the index to the next parameter
            paramIdx = paramIdx + 1;
        end
    end