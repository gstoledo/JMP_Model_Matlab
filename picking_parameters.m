%Local file to vizualize the calibration outcomes and pick best parameters
clear all
close all
cd('/Users/gabrieltoledo/Library/CloudStorage/GoogleDrive-gabrielstoledo.gt@gmail.com/My Drive/PHD NYU/Labor Firm Structure/JMP_Model_Matlab')
addpath('/Users/gabrieltoledo/Library/CloudStorage/GoogleDrive-gabrielstoledo.gt@gmail.com/My Drive/PHD NYU/Labor Firm Structure/JMP_Model_Matlab')
addpath('/Users/gabrieltoledo/Library/CloudStorage/GoogleDrive-gabrielstoledo.gt@gmail.com/My Drive/PHD NYU/Labor Firm Structure/Python_Firm Structure/replication_HLMP/matlab/run_A4_revisit')


%Load the results
load('Calibration_outcomes_HPC/manual_calibthe5_qup.mat')
%Pick index if bottom 1% in distance
% Assume 'vec' is your long column vector


% Determine the number of elements corresponding to the bottom 1%
numBottom = floor(0.05 * numel(distances));

% Sort the vector and get the corresponding indices
[sortedVec, sortedIdx] = sort(distances);

% Get the indices of the bottom 1% values
bottomIndices = sortedIdx(1:numBottom);


combinations(bottomIndices,:) % Want to print this in the pdf
% This suggestes smaller values of costs are fittting a bit better

plot(moments(bottomIndices,1),dm(1))


figure;
plot(moments(bottomIndices, 1), 'o-', 'DisplayName', 'Moments');
% Add a horizontal line at the value of dm(1)
hold on;
yline(dm(1), '--r', 'DisplayName', sprintf('dm(1) = %.2f', dm(1)));
% Label the axes
xlabel('Index (Bottom 1% Rows)');
% Set the custom x-tick labels using combinations(bottomIndices)
xticks(1:length(bottomIndices));  % Set x-ticks to match the number of bottom indices
xticklabels(combinations(bottomIndices,:));  % Use the values of combinations(bottomIndices) as labels
ylabel('Moments');
% Add a legend
legend;
% Add title
title('Moments and Horizontal Line at dm(1)');

figure
plot(moments(bottomIndices, 2), 'o-', 'DisplayName', 'Moments');
% Add a horizontal line at the value of dm(1)
hold on;
yline(dm(2), '--r', 'DisplayName', sprintf('dm(2) = %.2f', dm(2)));
% Label the axes
xlabel('Index (Bottom 1% Rows)');
ylabel('Moments');
% Add a legend
legend;
% Add title
title('Moments and Horizontal Line at dm(2)');
% Show the plot
hold off;

%DM2 which is the demotion rate is pretty good



clear all
close all
load('Calibration_outcomes_HPC/qup/manual_calibqup.mat')
numBottom = floor(1 * numel(distances));
[sortedVec, sortedIdx] = sort(distances);
bottomIndices = sortedIdx(1:numBottom);
combinations(bottomIndices,:)

plot(moments(:,1))

moments(bottomIndices,1)

%Find indeces such that moments are close to the data moments
% Around 5% of data_mom(1)
bestInd_data_mom1=find(moments(:,1)>data_mom(1)*0.95 & moments(:,1)<data_mom(1)*1.05);
bestComb_data_mom1=combinations(bestInd_data_mom1,:);

%Single best index
thebestInd_data_mom=zeros(1,length(data_mom));
best_Comb_data_mom=zeros(length(data_mom),length(combinations));
best_moments_data_mom=zeros(length(data_mom),length(moments));
for i=1:length(data_mom)
    [~,thebestInd_data_mom(i)]=min(abs(moments(:,i) - data_mom(i)));
    best_Comb_data_mom(i,:)=combinations(thebestInd_data_mom(i),:);
    best_moments_data_mom(i,:)=moments(thebestInd_data_mom(i),:);
end

combinations(thebestInd_data_mom1,:)
moments(thebestInd_data_mom1,:)

figure;
plot(moments(bestInd_data_mom1,1), 'o-', 'DisplayName', 'Moments');
% plot(moments(bottomIndices, 1), 'o-', 'DisplayName', 'Moments');
% Add a horizontal line at the value of data_mom(1)
hold on;
yline(data_mom(1), '--r', 'DisplayName', sprintf('data_mom(1) = %.2f', data_mom(1)));
% Label the axes
xlabel('Index (Bottom 1% Rows)');
% Set the custom x-tick labels using combinations(bottomIndices)
xticks(1:length(bottomIndices));  % Set x-ticks to match the number of bottom indices
xticklabels(combinations(bottomIndices,:));  % Use the values of combinations(bottomIndices) as labels
ylabel('Moments');
% Add a legend
legend;
% Add title
title('Moments and Horizontal Line at data_mom(1)');
saveas(gcf,'Calibration_outcomes_HPC/qup/Post_NiMi.png')

figure
plot(moments(bottomIndices, 2), 'o-', 'DisplayName', 'Moments');
% Add a horizontal line at the value of data_mom(1)
hold on;
yline(data_mom(2), '--r', 'DisplayName', sprintf('data_mom(2) = %.4f', data_mom(2)));
% Label the axes
xlabel('Index (Bottom 1% Rows)');
ylabel('Moments');
% Add a legend
legend;
% Add title
title('Moments and Horizontal Line at data_mom(2)');
saveas(gcf,'Calibration_outcomes_HPC/qup/Post_MiNi.png')
% Show the plot
hold off;


% load('Calibration_outcomes/calibrationv0_mincosts.mat')




clear all;

% Load the data


% load('Calibration_outcomes_HPC/manual_calibalpha_m.mat')
% Determine the number of elements corresponding to the bottom 1%

numBottom = floor(0.05 * numel(distances));


% Sort the vector and get the corresponding indices
[sortedVec, sortedIdx] = sort(distances);
bottomIndices = sortedIdx(1:numBottom);
% moments_selection = {'ManWorkerRatio','b_cross'};

% moments
% combinations
% distances


% % Example grids (replace with your actual grid data)
% gridx = unique(combinations(:, 1));  % Unique values of grid x
% gridy = unique(combinations(:, 2));  % Unique values of grid y

% % Create a grid for the heatmap based on gridx and gridy dimensions
% [X, Y] = meshgrid(gridx, gridy);  % Create the full grid

% % Initialize the heatmap data grid
% heatmapData = nan(size(X));  % same size as grid

% % Loop through all grid points and match with combinations
% for i = 1:numel(X)
%     % Find which row in 'combinations' matches the (gridx, gridy) pair
%     row = find(combinations(:, 1) == X(i) & combinations(:, 2) == Y(i));
    
%     % If a match is found, assign the corresponding distance value
%     if ~isempty(row)
%         heatmapData(i) = moments(row,1);
%     end
% end

% % Plot the heatmap using gridx and gridy as the axes
% h = heatmap(gridx, gridy, heatmapData, 'Colormap', parula);

% % Flip the y-axis by reversing the order of gridy values
% h.YDisplayData = flip(gridy);
% % Add labels and title
% h.XLabel = 'Grid X';
% h.YLabel = 'Grid Y';
% h.Title = 'Distance Heatmap';


clear all
load('Calibration_outcomes/manual_calibcosts_lowcp.mat')
p_selection = {'cost_p', 'cost_d'};  % Column headers, modify if more columns are added
moments_selection = {'NiMi', 'MiNi'};  % Moment names, modify as needed
% Get the number of grids/columns
num_grids = size(combinations, 2);
% Loop through all unique pairs of grids
for i = 1:num_grids
    for j = i+1:num_grids
        
        % Select grid i and grid j
        grid1 = combinations(:, i);
        grid2 = combinations(:, j);
        
        % Get unique values for the current grid pair
        unique_grid1 = unique(grid1);
        unique_grid2 = unique(grid2);
        
        % Create a grid for the heatmap based on grid1 and grid2 dimensions
        [X, Y] = meshgrid(unique_grid1, unique_grid2);
        
        % Initialize the heatmap data grid
        heatmapData = nan(size(X));
        
        % Loop through all grid points and match with combinations
        for idx = 1:numel(X)
            row = find(grid1 == X(idx) & grid2 == Y(idx));
            if ~isempty(row)
                heatmapData(idx) = distances(row);
            end
        end
        
        % Plot the heatmap for the current grid pair
        h = heatmap(unique_grid1, unique_grid2, heatmapData, 'Colormap', parula);
        
        % Flip the y-axis for correct display
        h.YDisplayData = flip(unique_grid2);
        
        % Add axis labels and title dynamically based on p_selection
        h.XLabel = [p_selection{i} ' Axis'];
        h.YLabel = [p_selection{j} ' Axis'];
        h.Title = ['Heatmap for ' p_selection{i} ' vs ' p_selection{j}];
        
        % Optionally, pause between each heatmap or save the figure
        pause(1);  % Allows viewing each heatmap before the next one

        %Same but plotting the actual moments
        for mm=1:length(moments_selection)
            heatmapData = nan(size(X));
            for idx = 1:numel(X)
                row = find(grid1 == X(idx) & grid2 == Y(idx));
                if ~isempty(row)
                    heatmapData(idx) = moments(row,mm);
                end
            end
            h = heatmap(unique_grid1, unique_grid2, heatmapData, 'Colormap', parula);
            h.YDisplayData = flip(unique_grid2);
            h.XLabel = [p_selection{i} ' Axis'];
            h.YLabel = [p_selection{j} ' Axis'];
            h.Title = ['Heatmap for ' p_selection{i} ' vs ' p_selection{j} ' for ' moments_selection{mm}, ' Target= ' num2str(data_mom(mm))];
            pause(1);
        end
    end
end




% % Vector for column titles (p_selection) and moments selection
% p_selection = {'cost_p', 'cost_d'};  % Column headers, modify if more columns are added
% [~, ps] = param_selection(theta, p_selection);
% moments_selection = {'NiMi', 'MiNi'};  % Moment names, modify as needed

% Open the PDF file to save figures (overwrite any existing file)
filename = 'ResultsReport.pdf';
folder = 'Calibration_outcomes';
calibration_report(p_selection, moments_selection, combinations,data_mom ,moments, distances, 0.1, folder, filename)
calibration_report(p_selection,ps ,moments_selection, combinations, data_mom, moments, distances, 0.05, folder, filename)





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

% Prepare the text for the combinations table (bottom 1%)
titleText = strjoin(p_selection, '    ');  % Create headers for the table
formattedText = num2str(combinations(bottomIndices, :), '%.4f    ');  % Format table values

% Display the column titles for the combinations table
text(0.1, 0.35, 'Combinations Parameters in the Bottom 1%', 'FontSize', 12, 'FontWeight', 'bold');
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
    yline(dm(i), '--', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2, 'DisplayName', 'Data');
    xlabel('Combinations', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Moments', 'FontSize', 12, 'FontWeight', 'bold');
    
    % Customize x-axis to display combinations values as tick labels
    xticks(1:length(bottomIndices));
    xticklabels(xLabels);  % Set custom x-tick labels with combination values
    xtickangle(45);  % Rotate labels for better visibility
    
    % Customize title and subtitle dynamically based on moments_selection
    title(sprintf('Moments and Data for %s', moments_selection{i}), 'FontSize', 14, 'FontWeight', 'bold');
    % subtitle(sprintf('Analysis of Bottom 1%% of Distances - %s', moments_selection{i}), 'FontSize', 12);

    % Customize legend
    legend('Model', 'Data', 'Location', 'best', 'FontSize', 12);

    % Customize grid, box, and other aesthetics
    grid on;
    set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6);  % Thinner grid lines
    box off;  % Remove box lines

    % Export each figure to the PDF
    exportgraphics(gcf, filename, 'Append', true);  % Append to PDF
end

% Save the PDF and close all figures
close all;


