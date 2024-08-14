function plot_realloc_matcombine(H_matrices_set, titles, mainTitle, xlabelText, ylabelText, fontSize, markerShapes, markerColors, legendLabels, filename)
    % plot_hire_mat_grid: Function to overlay plots of multiple binary matrices
    % in a grid with 3 columns and as many rows as needed, with a unified legend below the grid.
    
    % Determine the number of matrices
    numMatrices = length(H_matrices_set);
    numCols = 3; % Fixed number of columns
    numRows = ceil(numMatrices / numCols); % Calculate the number of rows needed
    
    % Calculate positions for central x and y labels
    centerRow = ceil(numRows / 2); % Middle row
    centerCol = ceil(numCols / 2); % Middle column
    lastRow = numRows; % Last row
    
    % Create a tiled layout for subplots
    figure;
    t = tiledlayout(numRows + 1, numCols, 'Padding', 'compact', 'TileSpacing', 'compact');
    
    % Add the main title
    sg = sgtitle(mainTitle, 'FontSize', fontSize + 7, 'FontWeight', 'bold');
    sg.FontName = 'Helvetica';
    
    % Initialize arrays to store overall legend information
    overallLegendEntries = {};
    overallLegendHandles = [];
    
    % Loop over each set of matrices and create a subplot
    for k = 1:numMatrices
        % Extract the current set of matrices
        H_matrices = H_matrices_set{k};
        
        % Create a new subplot
        nexttile;
        
        % Loop over each matrix and overlay the scatter plots
        hold on;
        markerSize = 120; % Marker size for the scatter plot
        for i = 1:length(H_matrices)
            H_u = H_matrices{i};
            if any(H_u(:)) % Only plot if the matrix has at least one `1`
                [row, col] = find(H_u);
                h = scatter(col, row, markerSize, markerShapes{i}, 'filled', 'MarkerFaceColor', markerColors{i});
                
                % Only add to legend if this shape/color hasn't been added yet
                if ~ismember(legendLabels{i}, overallLegendEntries)
                    overallLegendEntries{end+1} = legendLabels{i}; % Add to overall legend
                    overallLegendHandles(end+1) = h; % Store handle for overall legend
                end
            end
        end
        hold off;
        
        % Set axis properties to achieve the desired orientation
        set(gca, 'YDir', 'normal');
        set(gca, 'XDir', 'normal');
        
        % Adjust the axis limits to match the matrix dimensions
        axis([0.5, size(H_matrices{1}, 2) + 0.5, 0.5, size(H_matrices{1}, 1) + 0.5]);
        
        % Explicitly set the tick labels to match the matrix indices
        xticks(1:size(H_matrices{1}, 2));
        yticks(1:size(H_matrices{1}, 1));
        
        % Set tick labels to show integers starting from 0
        xticklabels(0:size(H_matrices{1}, 2)-1);
        yticklabels(0:size(H_matrices{1}, 1)-1);
        
        % Increase the font size of tick labels
        set(gca, 'FontSize', fontSize-1);
        
        % Set the font to Helvetica
        set(gca, 'FontName', 'Helvetica');
        
        % Add titles with specified font sizes
        title(titles{k}, 'FontSize', fontSize+2, 'FontWeight', 'normal');
        
        % Conditionally add y and x labels only in the center row/column
        if mod(k-1, numCols) == 0 && floor((k-1)/numCols)+1 == centerRow % Middle row, first column
            ylabel(ylabelText, 'FontSize', fontSize+4, 'FontWeight', 'bold');
        end
        if mod(k-1, numCols) == centerCol-1 && ceil(k/numCols) == lastRow % Center column, middle row
            xlabel(xlabelText, 'FontSize', fontSize+4, 'FontWeight', 'bold');
        end
        
        % Add grid lines
        grid on;
    end
    
    % Create a tile for the unified legend
    legendTile = nexttile([1, numCols]); % Span the legend across all columns
    axis(legendTile, 'off'); % Turn off the axis for the legend tile
    
    % % Add the unified legend and set it to have 2 columns
    lgd = legend(legendTile, overallLegendHandles, overallLegendEntries, 'Orientation', 'horizontal', 'NumColumns', 2);
    set(lgd, 'FontSize', fontSize + 2, 'Location', 'southoutside'); % Adjust legend font size
    
    % Adjust the legend's position to move it up and center it
    lgd.Position(1) = 0.65 - lgd.Position(3) / 2; % Center horizontally
    lgd.Position(2) = lgd.Position(2) + 0.17; % Move up slightly
    
    % Adjust marker sizes in the legend
    legendIcons = findobj(lgd, 'type', 'patch');
    for i = 1:length(legendIcons)
        set(legendIcons(i), 'MarkerSize', 20); % Adjust marker size in the legend
    end
    
    % Adjust the figure layout and save
    figHeight = 500; % Adjusted for 16:9 ratio
    figWidth = 889;  % Typical width for 16:9 ratio
    set(gcf, 'Position', [100, 100, figWidth, figHeight]); % Set the figure size for 16:9 ratio
    
    % Use DataAspectRatio to maintain correct aspect ratios
    daspect([1 1 1]);
    pbaspect([16 9 1]);
    
    if ~isempty(filename)
        saveas(gcf, filename);
    end
end
