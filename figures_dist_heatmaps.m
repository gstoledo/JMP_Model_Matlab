clear all
close all
location="local";
% location="hpc";
cd('/Users/gabrieltoledo/Library/CloudStorage/GoogleDrive-gabrielstoledo.gt@gmail.com/My Drive/PHD NYU/Labor Firm Structure/JMP_Model_Matlab')
addpath('/Users/gabrieltoledo/Library/CloudStorage/GoogleDrive-gabrielstoledo.gt@gmail.com/My Drive/PHD NYU/Labor Firm Structure/JMP_Model_Matlab')
addpath('/Users/gabrieltoledo/Library/CloudStorage/GoogleDrive-gabrielstoledo.gt@gmail.com/My Drive/PHD NYU/Labor Firm Structure/Python_Firm Structure/replication_HLMP/matlab/run_A4_revisit')
% load('model_solutions_sp4local.mat')
load('model_solutions_sp3hpc.mat')
load('xmin_results_combined.mat')


%Check sum of the eplus final distributions
nplus=sum(eplus_edist)+sum(eplus_mdist,"all")+sum(eplus_ndist,"all")+sum(eplus_tdist,"all");
popplus=sum(eplus_udist)+sum(eplus_mdist,"all")+sum(eplus_ndist,"all")+2*sum(eplus_tdist,"all");


x=xmin;

tpts=size(Vm,2);
ats=size(Vm,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Predefine qualitative color palettes
load('color_style.mat'); % Ran in color_style.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Using imagesc
 % % Create a figure and axes
 % % figure;
 % ax = axes;

 % % Display heatmap using imagesc
 % imagesc(ax, E);

 % % Customize the colorbar
 % cb = colorbar(ax, 'southoutside');
 % cb.TickLength = 0;

 % % cb.Position(3)=0.5*cb.Position(3);
 % % cb.Position(1) = 0.25 + 0.5 * cb.Position(1);  % Center the colorbar
 % % cb.Position(2) = cb.Position(2) - 0.05;  % Move it further down


 % % Invert the direction of the vertical axis
 % set(ax, 'YDir', 'normal');

 % % Remove all gridlines
 % ax.XGrid = 'off';
 % ax.YGrid = 'off';
 % ax.TickLength = [0, 0];  % Remove tick marks
 % % Ensure all cells are labeled
 % ax.XTick = 1:size(E, 2);  % Label all columns
 % ax.YTick = 1:size(E, 1);  % Label all rows
 % % Reduce the size of the tick labels
 % ax.FontSize = 14;  % Adjust the font size as needed

 % % Customize the title and labels
 % title(ax, 'Firm distribution conditional on Productivity level a', 'FontSize', 20, 'FontName', 'Helvetica', 'FontWeight', 'normal');
 % xlabel(ax, 'Worker type', 'FontSize', 16, 'FontName', 'Helvetica');
 % ylabel(ax, 'Manager type', 'FontSize', 16, 'FontName', 'Helvetica');
 % % Apply the custom colormap
 % colormap(ax, bluePalette_rgb / 255);

 % % Print figure
 % print('-dpdf', 'basic4.pdf');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Using heatmap
 % clear h
 % % Create the heatmap
 % h = heatmap(E);% 
 % % Customize the heatmap
 % h.Title = 'Firm distribution conditional on Productivity level a';
 % h.XLabel = 'Worker type';
 % h.YLabel = 'Manager type';

 % % Flip the vertical axis
 % h.NodeChildren(3).YDir = 'normal';

 % % Remove cell labels
 % h.CellLabelColor = 'none';

 % % Remove grid lines
 % h.GridVisible = 'off';

 % % Customize font size and font name
 % h.FontSize = 20;
 % h.FontName = 'Helvetica';
 % % Suppress the warning temporarily
 % warning('off', 'MATLAB:structOnObject')

 % % Access the title and change its font weight
 % hs = struct(h);
 % hs.Axes.Title.FontWeight = 'normal';

 % % Reduce the size of the tick labels
 % hs.Axes.XAxis.FontSize = 16; % Adjust the font size as needed
 % hs.Axes.YAxis.FontSize = 16; % Adjust the font size as needed


 % % Access the colorbar and remove ticks and font size
 % hs.Colorbar.TickLength = 0;
 % hs.Colorbar.FontSize = 12;


 % % Set custom ticks for the colorbar at larger intervals
 % % Assume the range of the data in E is [minValue, maxValue]
 % minValue = min(E(:));
 % maxValue = max(E(:));
 % hs.Colorbar.Ticks = linspace(minValue, maxValue, 5); % Adjust the number of ticks as needed


 % % Re-enable the warning
 % warning('on', 'MATLAB:structOnObject')

 % % Apply the custom colormap
 % colormap(bluePalette_rgb);

 % % Print figure
 % print('-dpdf', 'basic2.pdf');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %Lets build our E matrix for each type of firm
E=zeros(tpts+1,tpts+1,ats);
for a=1:ats
    E(1,1,a)=eplus_edist(1,a);
    for z=2:tpts+1
        for q=2:tpts+1
            E(z,1,a)=eplus_mdist(a,z-1);
            E(1,q,a)=eplus_ndist(a,q-1);
            E(z,q,a)=eplus_tdist(a,z-1,q-1);
        end
    end
end

E1=E(:,:,1);
E2=E(:,:,2);
 
% %% Using heatmap (plain one)
%  clear h
%  % Create the heatmap
%  h = heatmap(E1);
%  % Customize the heatmap
%  h.Title = 'Firm distribution conditional on Productivity level a';
%  h.XLabel = 'Worker type';
%  h.YLabel = 'Manager type';

%  % Flip the vertical axis
%  h.NodeChildren(3).YDir = 'normal';

%  % Remove cell labels
%  h.CellLabelColor = 'none';

%  % Remove grid lines
%  h.GridVisible = 'off';

%  % Customize font size and font name
%  h.FontSize = 20;
%  h.FontName = 'Helvetica';

%  % Suppress the warning temporarily
%  warning('off', 'MATLAB:structOnObject')

%  % Access the title and change its font weight
%  hs = struct(h);
%  hs.Axes.Title.FontWeight = 'normal';

%  % Reduce the size of the tick labels
%  hs.Axes.XAxis.FontSize = 16; % Adjust the font size as needed
%  hs.Axes.YAxis.FontSize = 16; % Adjust the font size as needed

%  % Access the colorbar and remove ticks and font size
%  hs.Colorbar.TickLength = 0;
%  hs.Colorbar.FontSize = 12;
%  minValue = min(E(:));
%  maxValue = max(E(:));
%  hs.Colorbar.Ticks = linspace(minValue, maxValue, 5); % Adjust the number of ticks as needed

%  % Re-enable the warning
%  warning('on', 'MATLAB:structOnObject')

%  % Apply the custom colormap
%  colormap(bluepurplePalette_rgb);

%  % Print figure  
%  print('-dpdf', 'basicE1.pdf');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %Lets build our E matrix for each type of firm
%% Using heatmap (plain one) for each a type as a function
% clear h1 h2
% h1=heatmap_plain(E1, bluepurplePalette_rgb, 'Firm distribution conditional on Productivity level a', 16, 20,'basicE1.pdf');
% h2=heatmap_plain(E2, bluepurplePalette_rgb, 'Firm distribution conditional on Productivity level a', 16, 20,'basicE2.pdf');

% % In a loop for each a type
% for a=1:ats
%     clear h e s normalized
%     e=E(:,:,a);
%     s=sum(E(:,:,a),'all');
%     normalized=e/s;
%     h=heatmap_plain(normalized, bluepurplePalette_rgb, 'Firm distribution conditional on Productivity level a='+string(a)+'', 16, 20,'basicE_'+string(a)+'.pdf');
% end
    


%% Using heatmap (skeewed colorbar)
% h = heatmap(E1);

% % Customize the heatmap
% h.Title = 'Firm distribution conditional on Productivity level a';
% h.XLabel = 'Worker type';
% h.YLabel = 'Manager type';
% % Relabel the axes to start at 0
% h.XDisplayLabels = num2cell(0:(size(E1, 2) - 1));
% h.YDisplayLabels = num2cell(0:(size(E1, 1) - 1));

% % Flip the vertical axis
% h.NodeChildren(3).YDir = 'normal';

% % Remove cell labels
% h.CellLabelColor = 'none';

% % Remove grid lines
% h.GridVisible = 'off';

% % Customize font size and font name
% h.FontSize = 20;
% h.FontName = 'Helvetica';

% % Suppress the warning temporarily
% warning('off', 'MATLAB:structOnObject')

% % Access the title and change its font weight
% hs = struct(h);
% hs.Axes.Title.FontWeight = 'normal';

% % Reduce the size of the tick labels
% hs.Axes.XAxis.FontSize = 16; % Adjust the font size as needed
% hs.Axes.YAxis.FontSize = 16; % Adjust the font size as needed


% % Interpolate to create a smoother colormap with more granularity at lower values
% numColors = 100;
% lowPercentile = 0.85; % 90th percentile
% lowCutoff = prctile(E1(:), 100 * lowPercentile);
% customCmap = interp1(linspace(0, lowPercentile, size(bluePalette_rgb, 1)), bluePalette_rgb, linspace(0, lowPercentile, numColors));

% % Apply the custom colormap
% colormap(customCmap);

% % Adjust the color axis limits to be more sensitive to lower values
% caxis([min(E1(:)), lowCutoff]);
% % Access the colorbar and remove ticks and font size
% hs.Colorbar.TickLength = 0;
% hs.Colorbar.FontSize = 12;
% minValue = min(E1(:));
% maxValue = max(E1(:));
% % hs.Colorbar.Ticks = linspace(minValue, maxValue, 5); % Adjust the number of ticks as needed

% % Re-enable the warning
% warning('on', 'MATLAB:structOnObject')


% % Print figure  
% print('-dpdf', 'basicE1.pdf');



% %% Combined heatmaps
% %Normalized E1 and E2
% sum1=sum(E1,'all');
% sum2=sum(E2,'all');
% En1=E1/sum1;
% En2=E2/sum2;


% clear h1 h2
% % Determine the common color limits for both heatmaps
% minValue = min(min(En1(:)), min(En2(:)));
% maxValue = max(max(En1(:)), max(En2(:)));

% % Create a tiled layout for better control over the layout
% t = tiledlayout(1, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');

% % Plot the first heatmap
% nexttile;
% h1 = heatmap(En1);
% h1.Title = 'Heatmap of E1';
% %h1.XLabel = 'Worker type';
% h1.YLabel = 'Manager type';
% h1.NodeChildren(3).YDir = 'normal';
% h1.CellLabelColor = 'none';
% h1.GridVisible = 'off';
% h1.FontSize = 16;
% h1.FontName = 'Helvetica';
% h1.ColorLimits = [minValue, maxValue];
% h1.ColorbarVisible = 'off';
% % Relabel the axes to start at 0
% h1.XDisplayLabels = num2cell(0:(size(E1, 2) - 1));
% h1.YDisplayLabels = num2cell(0:(size(E1, 1) - 1));

% % Suppress the warning temporarily
% warning('off', 'MATLAB:structOnObject')
% hs1 = struct(h1);
% hs1.Axes.Title.FontWeight = 'normal';
% hs1.Axes.Title.FontSize = 20;
% hs1.Axes.XAxis.FontSize = 16;
% hs1.Axes.YAxis.FontSize = 16;
% warning('on', 'MATLAB:structOnObject')

% % Adjust the heatmap to be square
% originalUnits = h1.Units;
% h1.Units = 'centimeters';
% sz1 = size(h1.ColorData);
% h1.Position(3:4) = min(h1.Position(3:4)) * [1, 1];
% if sz1(1) > sz1(2)
%     h1.Position(3) = h1.Position(3) * (sz1(2) / sz1(1));
% else
%     h1.Position(4) = h1.Position(4) * (sz1(1) / sz1(2));
% end
% h1.Units = originalUnits;

% % Plot the second heatmap
% nexttile;
% h2 = heatmap(En2);
% h2.Title = 'Heatmap of E2';
% h2.NodeChildren(3).YDir = 'normal';
% h2.CellLabelColor = 'none';
% h2.GridVisible = 'off';
% h2.FontSize = 16;
% h2.FontName = 'Helvetica';
% h2.ColorLimits = [minValue, maxValue];
% % Relabel the axes to start at 0
% h2.XDisplayLabels = num2cell(0:(size(E2, 2) - 1));
% h2.YDisplayLabels = num2cell(0:(size(E2, 1) - 1));

% % Suppress the warning temporarily
% warning('off', 'MATLAB:structOnObject')
% hs2 = struct(h2);
% hs2.Colorbar.TickLength = 0;
% hs2.Axes.Title.FontWeight = 'normal';
% hs2.Axes.Title.FontSize = 20;
% hs2.Axes.XAxis.FontSize = 16;
% hs2.Axes.YAxis.FontSize = 16;
% warning('on', 'MATLAB:structOnObject')

% % Adjust the second heatmap to be square
% originalUnits = h2.Units;
% h2.Units = 'centimeters';
% sz2 = size(h2.ColorData);
% h2.Position(3:4) = min(h2.Position(3:4)) * [1, 1];
% if sz2(1) > sz2(2)
%     h2.Position(3) = h2.Position(3) * (sz2(2) / sz2(1));
% else
%     h2.Position(4) = h2.Position(4) * (sz2(1) / sz2(2));
% end
% h2.Units = originalUnits;

% % Apply the custom colormap to both heatmaps
% colormap(bluepurplePalette_rgb);

% % Add a centralized xlabel for the entire layout
% xlabel(t, 'Worker type', 'FontSize', 16, 'FontName', 'Helvetica');

% % Print the figure
% print('-dpdf', 'sideBySideHeatmaps.pdf');


heatmap_combined(E,1,bluepurplePalette_rgb, ats,tpts, 'Firm Distribution Conditional on a' , 'Worker type', 'Manager type', 'Figures_heatmaps/dist_cond_a.pdf');


%% Same can be dine with the value functions
V=zeros(tpts+1,tpts+1,ats);

for a=1:ats
    V(1,1,a)=Ve(1,a);
    for z=2:tpts+1
        for q=2:tpts+1
            V(z,1,a)=Vm(a,z-1);
            V(1,q,a)=Vn(a,q-1);
            V(z,q,a)=Vt(a,z-1,q-1);
        end
    end
end

heatmap_combined(V,0,bluePalette_rgb, ats,tpts, 'Value Functions Conditional on a' , 'Worker type', 'Manager type', 'Figures_heatmaps/VFs_cond_a.pdf');