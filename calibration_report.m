function calibration_report(p_selection,ps ,moments_selection, combinations, data_mom, moments, distances, percentageBottom, folder, filename)
    % This function generates a PDF report for the calibration results.
    % Inputs:
    % - p_selection: A cell array with column headers (e.g., {'cost_p', 'cost_d'})
    % - moments_selection: A cell array with moment names (e.g., {'NiMi', 'MiNi'})
    % - combinations: A matrix containing the combinations of parameters
    % - data_mom: A vector of data moments
    % - moments: A matrix containing the moments to be plotted
    % - distances: A vector of distances used to find the bottom X percent
    % - percentageBottom: A number representing the percentage of bottom elements to include (e.g., 0.1 for bottom 10%)
    % - folder: The folder where the report will be saved
    % - filename: The base name for the output file (without extension)

    % Calculate the number of elements corresponding to the bottom percentage
    numBottom = floor(percentageBottom * numel(distances));
    if numBottom==0
        numBottom=numel(distances);
    end


    % Sort the vector and get the corresponding indices
    [~, sortedIdx] = sort(distances);
    bottomIndices = sortedIdx(1:numBottom);

    % Ensure the folder exists
    if ~isfolder(folder)
        mkdir(folder);
    end

    % Create the full file path with the PDF extension
    fullFilePath = fullfile(folder, [filename, '.pdf']);

     %% --- First Page: Display calibration variables from the struct 'ps' ---
    figure;
    axis off;  % Hide axis for clean presentation

    % Main title at the top of the figure
    text(0.5, 1.0, 'Manual Calibration Report', 'FontSize', 24, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    text(0.5, 0.9, datestr(now, 'mmmm dd, yyyy'), 'FontSize', 14, 'HorizontalAlignment', 'center');  % Current date

    % Display the preset calibration variables from the struct
    text(0.1, 0.85, 'Preset Calibration Variables', 'FontSize', 12, 'FontWeight', 'bold');

    % Extract field names and values from the struct 'ps'
    fields = fieldnames(ps);
    numFields = length(fields);

    % Define vertical space based on the number of fields to avoid overlap
    startPosY = 0.85;    % Starting Y position
    lineSpacing = min(0.05, 0.85/numFields);  % Dynamically adjust line spacing based on field count

    % Loop through each field in the struct and display its value
    for i = 1:numFields
        fieldName = fields{i};
        fieldValue = getfield(ps, fieldName);
        % Display each field name and value
        text(0.1, startPosY - i*lineSpacing, sprintf('%s: %s', fieldName, mat2str(fieldValue)), 'FontSize', 10, 'FontName', 'Courier');
    end

    % Export this page to the PDF
    exportgraphics(gcf, fullFilePath, 'Append', false);  % Overwrite the PDF file


    %% --- Second Page: Display unique values ---
    figure;
    axis off;  % Hide axis for clean presentation

    % Main title at the top of the figure
    text(0.5, 1.0, 'Manual Calibration', 'FontSize', 24, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    text(0.5, 0.8, datestr(now, 'mmmm dd, yyyy'), 'FontSize', 14, 'HorizontalAlignment', 'center');  % Current date

    % Loop over columns in 'combinations' to get unique values and display
    maxLen = 0;  % Track the max length of unique values for padding
    uniqueTable = {};  % Store unique values for each column

    for i = 1:length(p_selection)
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
    exportgraphics(gcf, fullFilePath, 'Append', true);  % Overwrite the PDF file

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
    exportgraphics(gcf, fullFilePath, 'Append', true);  % Append to the PDF

   %% --- Third Page: Display the best combinations for each moment ---
    figure;
    axis off;

    % Title
    text(0.5, 1.0, 'Best Combinations by Moment', 'FontSize', 24, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

    % Find the best indices for each moment
    bestIndices = zeros(length(moments_selection), 1);
    for i = 1:length(moments_selection)
        [~, bestIndices(i)] = min(abs(moments(:, i) - data_mom(i)));
    end

    % Extract the best combinations
    bestCombinations = combinations(bestIndices, :)';

    % Display moment names (column headers)
    columnHeaders = strjoin(moments_selection, '    ');
    text(0.5, 0.85, columnHeaders, 'FontSize', 12, 'FontName', 'Courier', 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

    % Display parameter names (row headers) and the corresponding best combinations
    rowHeaders = p_selection;  % use p_selection for row labels
    numParams = length(p_selection);
    startY = 0.75;  % Start a little below the headers
    lineSpacing = 0.05;  % Control spacing between lines

    for i = 1:numParams
        % Display row headers (parameter names)
        text(0.1, startY - (i-1)*lineSpacing, rowHeaders{i}, 'FontSize', 10, 'FontName', 'Courier', 'FontWeight', 'bold');
        
        % Display the corresponding best combination values in each row
        formattedRow = num2str(bestCombinations(i, :), '%.4f    ');
        text(0.3, startY - (i-1)*lineSpacing, formattedRow, 'FontSize', 10, 'FontName', 'Courier');
    end

    % Export the third page to the PDF
    exportgraphics(gcf, fullFilePath, 'Append', true);


    %% --- Check for rows with NaN in moments and display the failed combinations ---
    if any(isnan(moments), 'all')  % Check if there are NaN values in the moments matrix
        % Find the indices of rows where at least one moment is NaN
        failedIndices = find(any(isnan(moments), 2));
        
        if ~isempty(failedIndices)
            % Create a new figure for the "Parameters that Failed to Converge" page
            figure;
            axis off;

            % Title
            text(0.5, 1.0, 'Parameters that Failed to Converge', 'FontSize', 24, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

            % Loop over each failed combination and display its values
            startY = 0.85;  % Starting Y position
            lineSpacing = 0.05;  % Control spacing between lines
            for i = 1:length(failedIndices)
                rowIdx = failedIndices(i);
                failedCombination = combinations(rowIdx, :);

                % Display the combination label and values
                text(0.1, startY - (i-1)*lineSpacing, sprintf('Cb %d:', rowIdx), 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Courier');
                formattedCombination = num2str(failedCombination, '%.4f    ');

                % Check for NaN in the combination and highlight it
                isNaNInCombination = isnan(failedCombination);
                if any(isNaNInCombination)
                    % Convert the combination to a string and replace NaN values with 'NaN'
                    combinationStr = cellstr(num2str(failedCombination, '%.4f    '));
                    combinationStr(isNaNInCombination) = {'NaN'};
                    text(0.3, startY - (i-1)*lineSpacing, strjoin(combinationStr, '    '), 'FontSize', 10, 'FontName', 'Courier');
                else
                    % Display the combination normally if no NaNs (backup safety)
                    text(0.3, startY - (i-1)*lineSpacing, formattedCombination, 'FontSize', 10, 'FontName', 'Courier');
                end
            end

            % Export this page to the PDF
            exportgraphics(gcf, fullFilePath, 'Append', true);
        end
    end
        
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

        % Customize legend
        legend('Model', 'Data', 'Location', 'best', 'FontSize', 12);

        % Customize grid, box, and other aesthetics
        grid on;
        set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6);  % Thinner grid lines
        box off;  % Remove box lines

        % Export each figure to the PDF
        exportgraphics(gcf, fullFilePath, 'Append', true);  % Append to PDF
    end

    if length(p_selection) > 1
        num_grids = size(combinations, 2);
        for i = 1:num_grids
            for j = i + 1:num_grids
                % Heatmap for distances
                grid1 = combinations(:, i);
                grid2 = combinations(:, j);
                unique_grid1 = unique(grid1);
                unique_grid2 = unique(grid2);
                [X, Y] = meshgrid(unique_grid1, unique_grid2);
                heatmapData = nan(size(X));
                
                % Populate the heatmap data
                for idx = 1:numel(X)
                    rows = find(grid1 == X(idx) & grid2 == Y(idx));
                    if ~isempty(rows)
                        % Handle multiple rows by taking the mean of distances (or use sum, min, max, etc.)
                        heatmapData(idx) = mean(distances(rows));
                    end
                end
                
                % Plot the heatmap for the current grid pair
                figure;
                h = heatmap(unique_grid1, unique_grid2, heatmapData, 'Colormap', cool);
                h.YDisplayData = flip(unique_grid2);
                h.CellLabelColor = 'none';
                h.XLabel = [p_selection{i} ' Axis'];
                h.YLabel = [p_selection{j} ' Axis'];
                h.Title = ['Heatmap for ' p_selection{i} ' vs ' p_selection{j} ' for Distances'];
                exportgraphics(gcf, fullFilePath, 'Append', true);

                % Heatmap for moments
                for mm = 1:length(moments_selection)
                    heatmapData = nan(size(X));
                    for idx = 1:numel(X)
                        rows = find(grid1 == X(idx) & grid2 == Y(idx));
                        if ~isempty(rows)
                            % Handle multiple rows by taking the mean of moments
                            heatmapData(idx) = mean(moments(rows, mm));
                        end
                    end
                    figure;
                    h = heatmap(unique_grid1, unique_grid2, heatmapData, 'Colormap', cool);
                    h.YDisplayData = flip(unique_grid2);
                    h.XLabel = [p_selection{i} ' Axis'];
                    h.CellLabelColor = 'none';  % Hide cell labels for cleaner heatmap
                    h.YLabel = [p_selection{j} ' Axis'];
                    h.Title = ['Heatmap for ' p_selection{i} ' vs ' p_selection{j} ' for ' moments_selection{mm}, ' Target= ' num2str(data_mom(mm))];
                    exportgraphics(gcf, fullFilePath, 'Append', true);
                end
            end
        end
    end
    
    % Close all figures
    close all;
