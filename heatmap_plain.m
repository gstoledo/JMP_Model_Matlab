function [h]=heatmap_plain(E1, bluePalette_rgb, titleStr, axisFontSize, titleFontSize,savename)
    % Function to create and customize a heatmap

    % Create the heatmap
    h = heatmap(E1);

    % Customize the heatmap
    h.Title = titleStr;
    h.XLabel = 'Worker type';
    h.YLabel = 'Manager type';
    h.XDisplayLabels = num2cell(0:(size(E1, 2) - 1));
    h.YDisplayLabels = num2cell(0:(size(E1, 1) - 1));

    % Flip the vertical axis
    h.NodeChildren(3).YDir = 'normal';

    % Remove cell labels
    h.CellLabelColor = 'none';

    % Remove grid lines
    h.GridVisible = 'off';

    % Customize font size and font name
    h.FontSize = axisFontSize;
    h.FontName = 'Helvetica';
    

    % Suppress the warning temporarily
    warning('off', 'MATLAB:structOnObject')

    % Access the title and change its font weight
    hs = struct(h);
    hs.Axes.Title.FontWeight = 'normal';
    hs.Axes.Title.FontSize = titleFontSize;

    % Reduce the size of the tick labels
    hs.Axes.XAxis.FontSize = axisFontSize; % Adjust the font size as needed
    hs.Axes.YAxis.FontSize = axisFontSize; % Adjust the font size as needed

    % Access the colorbar and remove ticks and font size
    hs.Colorbar.TickLength = 0;
    hs.Colorbar.FontSize = 12;
    minValue = min(E1(:));
    maxValue = max(E1(:));
    hs.Colorbar.Ticks = linspace(minValue, maxValue, 5); % Adjust the number of ticks as needed

    % Re-enable the warning
    warning('on', 'MATLAB:structOnObject')

    % Apply the custom colormap
    colormap(bluePalette_rgb);

    % Print figure
    print('-dpdf', savename);
    % Save the figure
    if ~isempty(savename)
        print('-dpdf', savename);
    end
end
