function plot_hire_matcombine(H_u, titleText, xlabelText, ylabelText, fontSize, markerShape, markerColor, filename, subplotPos)
    % plot_hire_mat: Function to plot a binary matrix with customizable options
    % with support for subplotting
    %
    % Inputs:
    % H_u        - Binary matrix to be plotted
    % titleText  - Title of the plot
    % xlabelText - X-axis label
    % ylabelText - Y-axis label
    % fontSize   - Font size for the title and labels
    % markerShape- Shape of the scatter plot marker (e.g., 'o', 's', '^')
    % markerColor- Color of the scatter plot marker (e.g., 'k' for black, 'r' for red)
    % filename   - Filename to save the plot (e.g., 'plot.png')
    % subplotPos - Position for subplot (e.g., [2,2,1] for 2x2 grid, first plot)

    % Create the subplot
    subplot(subplotPos(1), subplotPos(2), subplotPos(3));
    
    % Find the indices of the ones in the matrix
    [row, col] = find(H_u);

    % Create the scatter plot
    scatter(col, row, 100, markerShape, 'filled', 'MarkerFaceColor', markerColor); 

    % Set axis properties to achieve the desired orientation
    set(gca, 'YDir', 'normal'); % Ensures (1,1) is in the lower-left corner
    set(gca, 'XDir', 'normal'); % Ensures X axis goes from left to right

    % Adjust the axis limits to match the matrix dimensions
    axis([0.5, size(H_u, 2) + 0.5, 0.5, size(H_u, 1) + 0.5]);

    % Explicitly set the tick labels to match the matrix indices
    xticks(1:size(H_u, 2));
    yticks(1:size(H_u, 1));

    % Set tick labels to show integers starting from 0
    xticklabels(0:size(H_u, 2)-1);
    yticklabels(0:size(H_u, 1)-1);

    % Set the font to Helvetica
    set(gca, 'FontName', 'Helvetica');

    % Add titles and axis labels with specified font sizes
    title(titleText, 'FontSize', fontSize, 'FontWeight', 'bold');
    xlabel(xlabelText, 'FontSize', fontSize, 'FontWeight', 'normal');
    ylabel(ylabelText, 'FontSize', fontSize, 'FontWeight', 'normal');

    % Add grid lines
    grid on;

    % Save the figure if a filename is provided
    if ~isempty(filename)
        saveas(gcf, filename);
    end
end
