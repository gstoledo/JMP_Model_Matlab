%Function to recall data moments and use to import in the calibration routine, given a moments selection

function data_moments = data_moments(dm, moments_selection)
    % Initialize an empty vector to hold the selected data moments
    data_moments = zeros(1, length(moments_selection));
    
    % Loop over the fields specified in moments_selection
    for i = 1:length(moments_selection)
        field_name = moments_selection{i};  % Get the field name from the selection
        data_moments(i) = dm.(field_name);  % Access the corresponding field in the struct
    end
end
