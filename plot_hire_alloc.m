function plot_hire_alloc(H_matrices, titleText, xlabelText, ylabelText, fontSize, markerShapes, markerColors, legendLabels, filename)
    % plot_hire_mat_overlay: Function to overlay plots of multiple binary matrices
    % with a legend positioned outside the plot at the bottom in a 2-column grid.
    %
    % Inputs:
    % H_matrices - Cell array containing exactly 6 binary matrices to be plotted
    % titleText  - Title of the plot
    % xlabelText - X-axis label
    % ylabelText - Y-axis label
    % fontSize   - Font size for the title and labels
    % markerShapes - Cell array of shapes for the scatter plot markers (e.g., {'o', 's', '^', '*'})
    % markerColors - Cell array of colors for the scatter plot markers (e.g., {'k', 'r', 'b'})
    % legendLabels - Cell array of custom legend labels (e.g., {'Label1', 'Label2', ...})
    % filename   - Filename to save the plot (e.g., 'plot.png')

    % % Ensure there are exactly 6 matrices
    % assert(length(H_matrices) == 6, 'You must provide exactly 6 matrices.');
    % assert(length(legendLabels) == 6, 'You must provide exactly 6 legend labels.');

    % Create the figure
    figure;

    % Initialize an array to store legend entries
    legendEntries = {};
    legendHandles = [];

    % Loop over each matrix and overlay the scatter plots
    hold on;
    markerSize = 150; % Increased marker size for the plot
    for i = 1:length(H_matrices)
        H_u = H_matrices{i};
        if any(H_u(:)) % Only plot if the matrix has at least one `1`
            [row, col] = find(H_u);
            h = scatter(col, row, markerSize, markerShapes{i}, 'filled', 'MarkerFaceColor', markerColors{i});
            legendEntries{end+1} = legendLabels{i}; % Use provided legend label
            legendHandles(end+1) = h; % Store handle for legend
        end
    end
    hold off;

    % Set axis properties to achieve the desired orientation
    set(gca, 'YDir', 'normal'); % Ensures (1,1) is in the lower-left corner
    set(gca, 'XDir', 'normal'); % Ensures X axis goes from left to right

    % Adjust the axis limits to match the matrix dimensions (assuming all matrices are the same size)
    axis([0.5, size(H_matrices{1}, 2) + 0.5, 0.5, size(H_matrices{1}, 1) + 0.5]);

    % Explicitly set the tick labels to match the matrix indices
    xticks(1:size(H_matrices{1}, 2));
    yticks(1:size(H_matrices{1}, 1));

    % Set tick labels to show integers starting from 0
    xticklabels(0:size(H_matrices{1}, 2)-1);
    yticklabels(0:size(H_matrices{1}, 1)-1);

    % Increase the font size of tick labels
    set(gca, 'FontSize', fontSize);

    % Set the font to Helvetica
    set(gca, 'FontName', 'Helvetica');

    % Add titles and axis labels with specified font sizes
    title(titleText, 'FontSize', fontSize + 4, 'FontWeight', 'bold'); % Slightly larger title
    xlabel(xlabelText, 'FontSize', fontSize, 'FontWeight', 'normal');
    ylabel(ylabelText, 'FontSize', fontSize, 'FontWeight', 'normal');

    % Add grid lines
    grid on;

    % Add a legend only if there are entries
    if ~isempty(legendEntries)
        % Create the legend with a 2-column layout
        legendHandle = legend(legendHandles, legendEntries, 'Orientation', 'horizontal', 'Location', 'southoutside', 'NumColumns', 2);
        set(legendHandle, 'FontSize', fontSize); % Set the legend font size

        % Increase the marker size in the legend
        % Directly modify the MarkerSize of the handles
        for i = 1:length(legendHandles)
            set(findobj(legendHandle, 'type', 'line', '-and', 'Marker', markerShapes{i}), 'MarkerSize', 20); % Increase marker size
        end

        % Adjust the figure layout to accommodate the legend outside
        set(gca, 'Position', [0.1, 0.35, 0.8, 0.55]); % Adjust the plot to make room for the legend
    end

    % Save the figure if a filename is provided
    if ~isempty(filename)
        saveas(gcf, filename);
    end
end
