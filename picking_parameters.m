%Local file to vizualize the calibration outcomes and pick best parameters
clear all
close all
cd('/Users/gabrieltoledo/Library/CloudStorage/GoogleDrive-gabrielstoledo.gt@gmail.com/My Drive/PHD NYU/Labor Firm Structure/JMP_Model_Matlab')
addpath('/Users/gabrieltoledo/Library/CloudStorage/GoogleDrive-gabrielstoledo.gt@gmail.com/My Drive/PHD NYU/Labor Firm Structure/JMP_Model_Matlab')
addpath('/Users/gabrieltoledo/Library/CloudStorage/GoogleDrive-gabrielstoledo.gt@gmail.com/My Drive/PHD NYU/Labor Firm Structure/Python_Firm Structure/replication_HLMP/matlab/run_A4_revisit')


%Load the results
load('Calibration_outcomes/calibrationv0.mat')
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


%Load the results
load('Calibration_outcomes/manual_calibcosts.mat')
numBottom = floor(0.05 * numel(distances));
[sortedVec, sortedIdx] = sort(distances);
bottomIndices = sortedIdx(1:numBottom);
combinations(bottomIndices,:)


plot(moments(:,1))

figure;
plot(moments(bottomIndices, 1), 'o-', 'DisplayName', 'Moments');
% Add a horizontal line at the value of data_mom(1)
hold on;
yline(data_mom(1), '--r', 'DisplayName', sprintf('dm(1) = %.2f', data_mom(1)));
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

figure
plot(moments(bottomIndices, 2), 'o-', 'DisplayName', 'Moments');
% Add a horizontal line at the value of data_mom(1)
hold on;
yline(data_mom(2), '--r', 'DisplayName', sprintf('dm(2) = %.2f', data_mom(2)));
% Label the axes
xlabel('Index (Bottom 1% Rows)');
ylabel('Moments');
% Add a legend
legend;
% Add title
title('Moments and Horizontal Line at data_mom(2)');
% Show the plot
hold off;


% load('Calibration_outcomes/calibrationv0_mincosts.mat')




clear all;

% Load the data
load('Calibration_outcomes/calibrationv0.mat')

% Determine the number of elements corresponding to the bottom 1%
numBottom = floor(0.1 * numel(distances));

% Sort the vector and get the corresponding indices
[sortedVec, sortedIdx] = sort(distances);
bottomIndices = sortedIdx(1:numBottom);

% Vector for column titles (p_selection) and moments selection
p_selection = {'cost_p', 'cost_d'};  % Column headers, modify if more columns are added
moments_selection = {'NiMi', 'MiNi'};  % Moment names, modify as needed

% Open the PDF file to save figures (overwrite any existing file)
filename = 'ResultsReport.pdf';
calibration_report(p_selection, moments_selection, combinations,dm ,moments, distances, 0.05, filename)

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
