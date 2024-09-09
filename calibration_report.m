function calibration_report(p_selection, moments_selection, combinations, data_mom ,moments, distances, percentageBottom, filename)
    % This function generates a PDF report for the calibration results.
    % Inputs:
    % - p_selection: A cell array with column headers (e.g., {'cost_p', 'cost_d'})
    % - moments_selection: A cell array with moment names (e.g., {'NiMi', 'MiNi'})
    % - combinations: A matrix containing the combinations of parameters
    % - dm: A vector of data moments
    % - moments: A matrix containing the moments to be plotted
    % - distances: A vector of distances used to find the bottom X percent
    % - percentageBottom: A number representing the percentage of bottom elements to include (e.g., 0.1 for bottom 10%)

    % Calculate the number of elements corresponding to the bottom percentage
    numBottom = floor(percentageBottom * numel(distances));

    % Sort the vector and get the corresponding indices
    [~, sortedIdx] = sort(distances);
    bottomIndices = sortedIdx(1:numBottom);


    filename = 'CalibReport_'+string(filename)+'.pdf';
    %% --- First Page: Display unique values ---
    figure;
    axis off;  % Hide axis for clean presentation

    % Main title at the top of the figure
    text(0.5, 1.0, 'Manual Calibration', 'FontSize', 24, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    text(0.5, 0.8, datestr(now, 'mmmm dd, yyyy'), 'FontSize', 14, 'HorizontalAlignment', 'center');  % Current date

    % Loop over columns in 'combinations' to get unique values and display
    maxLen = 0;  % Track the max length of unique values for padding
    uniqueTable = {};  % Store unique values for each column

    for i = 1:length(p_selection)
        % Get unique values for column 'i'
        uniqueVals = unique(combinations(:, i));
        uniqueTable{i} = uniqueVals;  % Store unique values
        maxLen = max(maxLen, length(uniqueVals));  % Update max length
    end

    % Create a padded table for the unique values
    paddedTable = NaN(maxLen, length(p_selection));  % Initialize NaN matrix for padding
    for i = 1:length(uniqueTable)
        paddedTable(1:length(uniqueTable{i}), i) = uniqueTable{i};
    end

    % Format the padded table for display
    formattedUniqueTable = num2str(paddedTable, '%.4f    ');

    % Display column titles for the unique table
    text(0.1, 0.65, 'Grids Used', 'FontSize', 12, 'FontWeight', 'bold');
    columnHeaders = strjoin(p_selection, '    ');  % Create a string for headers
    text(0.1, 0.6, columnHeaders, 'FontSize', 10, 'FontName', 'Courier', 'FontWeight', 'bold');
    text(0.1, 0.45, formattedUniqueTable, 'FontSize', 10, 'FontName', 'Courier', 'VerticalAlignment', 'top');

    % Export the unique values table to the PDF
    exportgraphics(gcf, filename, 'Append', false);  % Overwrite the PDF file

    %% --- Second Page: Display the combinations table ---
    figure;
    axis off;  % Hide axis for clean presentation

    % Prepare the text for the combinations table (bottom X%)
    titleText = strjoin(p_selection, '    ');  % Create headers for the table
    formattedText = num2str(combinations(bottomIndices, :), '%.4f    ');  % Format table values

    % Display the column titles for the combinations table
    text(0.1, 0.35, 'Combinations Parameters in the Bottom Selection', 'FontSize', 12, 'FontWeight', 'bold');
    text(0.1, 0.3, titleText, 'FontSize', 10, 'FontName', 'Courier', 'FontWeight', 'bold');
    text(0.1, 0.2, formattedText, 'FontSize', 10, 'FontName', 'Courier', 'VerticalAlignment', 'top');

    % Export the combinations table figure to the PDF
    exportgraphics(gcf, filename, 'Append', true);  % Append to the PDF

    %% --- Loop over each element in moments_selection and plot corresponding moments ---
    for i = 1:length(moments_selection)
        % Convert combinations to strings for x-tick labels
        xLabels = arrayfun(@(j) sprintf('(%.2f, %.2f)', combinations(bottomIndices(j), 1), combinations(bottomIndices(j), 2)), ...
                           1:length(bottomIndices), 'UniformOutput', false);

        % Create a figure for each moment
        figure;
        plot(moments(bottomIndices, i), 'o-', 'Color', [0 0.4470 0.7410], 'DisplayName', 'Model', 'LineWidth', 2, 'MarkerSize', 6);
        hold on;
        yline(data_mom(i), '--', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2, 'DisplayName', 'Data');
        xlabel('Combinations', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('Moments', 'FontSize', 12, 'FontWeight', 'bold');

        % Customize x-axis to display combinations values as tick labels
        xticks(1:length(bottomIndices));
        xticklabels(xLabels);  % Set custom x-tick labels with combination values
        xtickangle(45);  % Rotate labels for better visibility

        % Customize title and subtitle dynamically based on moments_selection
        title(sprintf('Moments and Data for %s', moments_selection{i}), 'FontSize', 14, 'FontWeight', 'bold');
        % subtitle(sprintf('Analysis of Bottom Selection - %s', moments_selection{i}), 'FontSize', 12);

        % Customize legend
        legend('Model', 'Data', 'Location', 'best', 'FontSize', 12);

        % Customize grid, box, and other aesthetics
        grid on;
        set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6);  % Thinner grid lines
        box off;  % Remove box lines

        % Export each figure to the PDF
        exportgraphics(gcf,filename, 'Append', true);  % Append to PDF
    end

    % Save the PDF and close all figures
    close all;
end
