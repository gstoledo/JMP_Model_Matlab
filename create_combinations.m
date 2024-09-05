function combinations = create_combinations(cgrid, p_selection, theta)
    % Initialize the grids based on the selection of parameters
    grids = cell(1, numel(p_selection));

    % Loop through the selected parameters
    for i = 1:numel(p_selection)
        param = p_selection{i};
        
        % For 'qup', handle cell arrays, as each element has a separate grid
        if strcmp(param, 'qup')
            % Create grids for each element of theta.qup
            for j = 1:numel(theta.qup)
                grids{end+1} = cgrid.qup{j};  % Add each grid element to the grids cell
            end
        elseif isfield(cgrid, param)
            % For other parameters, directly use the grid from cgrid
            grids{i} = cgrid.(param);
        else
            error('Parameter %s not found in cgrid', param);
        end
    end

    % Remove empty cells from grids (caused by qup's insertion)
    grids = grids(~cellfun('isempty', grids));

    % Generate the grid combinations using ndgrid
    [gridArrays{1:numel(grids)}] = ndgrid(grids{:});

    % Reshape the grid arrays to column vectors and concatenate them
    combinationCells = cellfun(@(x) x(:), gridArrays, 'UniformOutput', false);
    combinations = [combinationCells{:}];
end
