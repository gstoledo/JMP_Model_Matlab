function heatmap_combined(E, normalize, pallete, ats, tpts, Title, xlabel, ylabel, filename)
    % Normalize the matrix E and calculate min and max values
    if normalize == 1    
        En = zeros(tpts+1, tpts+1, ats);
        for a = 1:ats
            s = sum(E(:,:,a), 'all');
            En(:,:,a) = E(:,:,a) / s;
        end
    else
        En = E;
    end 

    minValue = min(En(:));
    maxValue = max(En(:));

    % Determine number of rows and columns for the layout
    numCols = min(ats, 3); % Max 3 columns per row
    numRows = ceil(ats / numCols);

    % Create a tiled layout for better control over the layout
    t = tiledlayout(numRows, numCols, 'TileSpacing', 'Compact', 'Padding', 'Compact');

    % Create a finer palette by interpolating the existing palette
    numShades = 1000; % Increase this to get more shades
    finerPalette = interp1(linspace(0, 1, size(pallete, 1)), pallete, linspace(0, 1, numShades));

    for a = 1:ats
        nexttile;
        h1 = heatmap(En(:,:,a));
        h1.Title = 'a=' + string(a);
        h1.XLabel = xlabel;
        h1.YLabel = ylabel;
        h1.NodeChildren(3).YDir = 'normal';
        h1.CellLabelColor = 'none';
        h1.GridVisible = 'off';
        h1.FontSize = 16;
        h1.FontName = 'Helvetica';
        h1.ColorLimits = [minValue, maxValue];

        if a < ats
            h1.ColorbarVisible = 'off';
        end

        % Relabel the axes to start at 0
        h1.XDisplayLabels = num2cell(0:(size(En(:,:,a), 2) - 1));
        h1.YDisplayLabels = num2cell(0:(size(En(:,:,a), 1) - 1));

        % Suppress the warning temporarily
        warning('off', 'MATLAB:structOnObject')
        hs1 = struct(h1);
        hs1.Axes.Title.FontWeight = 'normal';
        hs1.Axes.Title.FontSize = 12;
        hs1.Axes.XAxis.FontSize = 12;
        hs1.Axes.YAxis.FontSize = 12;
        warning('on', 'MATLAB:structOnObject')

        % Adjust the heatmap to be square
        originalUnits = h1.Units;
        h1.Units = 'centimeters';
        sz1 = size(h1.ColorData);
        h1.Position(3:4) = min(h1.Position(3:4)) * [1, 1];
        if sz1(1) > sz1(2)
            h1.Position(3) = h1.Position(3) * (sz1(2) / sz1(1));
        else
            h1.Position(4) = h1.Position(4) * (sz1(1) / sz1(2));
        end
        h1.Units = originalUnits;
        
        % Apply the finer colormap
        colormap(finerPalette);
    end
    
    % Add a title to the entire tiled layout
    title(t, Title,'FontSize', 20, 'FontName', 'Helvetica');
    
    % Adjust the figure's PaperPosition property for printing
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig.PaperUnits = 'centimeters';
    fig.PaperPosition = [0 0 30 20]; % Adjust these values as needed for better fit

    % Save the figure to a PDF file using the -bestfit option
    print(fig, filename, '-dpdf', '-bestfit');
end
