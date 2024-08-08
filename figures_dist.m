clear all
close all
cd('/Users/gabrieltoledo/Library/CloudStorage/GoogleDrive-gabrielstoledo.gt@gmail.com/My Drive/PHD NYU/Labor Firm Structure/JMP_Model_Matlab')
addpath('/Users/gabrieltoledo/Library/CloudStorage/GoogleDrive-gabrielstoledo.gt@gmail.com/My Drive/PHD NYU/Labor Firm Structure/JMP_Model_Matlab')
addpath('/Users/gabrieltoledo/Library/CloudStorage/GoogleDrive-gabrielstoledo.gt@gmail.com/My Drive/PHD NYU/Labor Firm Structure/Python_Firm Structure/replication_HLMP/matlab/run_A4_revisit')
load('model_solutions.mat')
load('xmin_results_combined.mat')


x=xmin;

tpts=size(Vm,2);
ats=size(Vm,1);

lamu                    =x(1);
% lamu=0.1;
lam                    =x(2);
% lam=0.1;
fcomp                   =x(3); %Complementarity in F
mubar                   =x(4);
del                     =x(5); %Delta 
ulose                   =x(6);
sologain                =x(7);
team_learn_param_up  	=x(8);
mnew	                =x(9);
team_learn_param_up_sym	=x(10);
bpf                     =x(11); %Firm bargaining power (1-gamma)
homeprod                =x(12);
sololose                =0;
A                       =x(14);
team_learn_param_down   =0;
mnew_high               =floor(x(15)); 



%Fundamentals
bt   =1/(1.1^(1/12))   ; %beta, discount factor
death=1/(35*12)        ; %Probability agent dies  
bpw  =1-bpf            ; %bpw is bargaining power of worker
nfirm=1                ; %Measure of firms 
% tpts =3                ; %type point space
% ats=2                  ; %Productivity space 
spts =tpts+2           ; %state S for updating -- {u,0,{j}} -- dimension of that space is type+2
cost_p=0                   ; %cost of promoting a non manager to manager
cost_d=0                   ; %cost of demoting a manager to non manager
alpha_m=1              ; %Manager returns
alpha_n=0.1              ; %Non manager returns

typemin=1; %Lowest type
typemax=mubar ; %Highest type

amin=1; %Lowest productivity
amax=2; %Highest productivity

type=[typemin:(typemax-typemin)/(tpts-1):typemax]; %Types
a_type=[amin:(amax-amin)/(ats-1):amax]; %Productivity


%Need to ajdust this later, pay attention to normalization to measure of workers 
mnew_high=tpts; %Number of types
typebirth=zeros(1,tpts);
for i=1:mnew_high
    typebirth(i)= mnew*exp(-type(i)*mnew)./sum(mnew*exp(-type(1:mnew_high)*mnew)); %Birth rate
end


soloprod=1    ;
pnew    =fcomp;
pold    =1    ;
 
fteam=zeros(ats,tpts,tpts);
fman=zeros(ats,tpts);
fnman=zeros(ats,tpts);
fe=zeros(1,ats);
for a=1:ats
    for z=1:tpts
        for q=1:tpts
            fteam(a,z,q)=A*a_type(a)*(type(z).^alpha_m)*(type(q).^alpha_n);
        end
        fman(a,z)=A*a_type(a)*type(z).^alpha_m;
        fnman(a,z)=0;
    end
    fe(a)=0;
end

b=homeprod*fman(1,:) ; %Type home production vector
n=1; %mass of firms 


%Transtion matrices (They are not perfect yet just to write it all down)
%A transition
aup=0.1;
adown=0.3;
astay=1-aup-adown;
a_trans=create_trans(adown,astay,aup,ats);

%Q transition
qup=0.2;
qdown=0.1;
qstay=1-qup-qdown;
q_trans=create_trans(qdown,qstay,qup,tpts);


%Unemp transition
ugain=0.00           ; %Probability unemployed move up
ustay=1-ulose-ugain  ; %Probability unemployed stay
u_trans=create_trans(ulose,ustay,ugain,tpts);


% Astethcics settings


%% Set default properties for 4:3 figure

figureWidth = 8.5;
figureHeight = 6.375;
set(groot, 'defaultFigureUnits', 'inches')
set(groot, 'defaultFigurePosition', [1,1,figureWidth,figureHeight]);
set(groot, 'defaultFigurePaperPosition', [0, 0, figureWidth,figureHeight]);
set(groot, 'defaultFigurePaperSize', [figureWidth,figureHeight]);

%% Set default properties for axes

set(groot, 'defaultAxesFontName', 'Helvetica')
set(groot, 'defaultAxesFontSize', 24)
set(groot, 'defaultAxesLabelFontSizeMultiplier', 1)
set(groot, 'defaultAxesTitleFontSizeMultiplier', 1)
set(groot, 'defaultAxesTitleFontWeight', 'normal')
set(groot, 'defaultAxesXColor', 'k')
set(groot, 'defaultAxesYColor', 'k')
set(groot, 'defaultAxesGridColor', 'k')
set(groot, 'defaultAxesLineWidth', 1)
set(groot, 'defaultAxesYGrid', 'on')
set(groot, 'defaultAxesXGrid', 'off')
set(groot, 'defaultAxesTickDirMode', 'manual')
set(groot, 'defaultAxesTickDir', 'out')
set(groot, 'defaultAxesTickLength',[0.005, 0.005])
set(groot, 'defaultAxesBox','off')

%% Predefine qualitative color palettes

% Dark colors

darkPalette = ['#1b9e77';'#d95f02';'#7570b3';'#e7298a';'#66a61e';'#e6ab02';'#a6761d';'#666666'];

greenColor = darkPalette(1,:);
orangeColor = darkPalette(2,:);
purpleColor = darkPalette(3,:);
pinkColor = darkPalette(4,:);
appleColor = darkPalette(5,:);
yellowColor = darkPalette(6,:);
brownColor = darkPalette(7,:);
grayColor = darkPalette(8,:);

% Paired colors

pairedPalette = ['#a6cee3';'#1f78b4';'#b2df8a';'#33a02c';'#fb9a99';'#e31a1c';'#fdbf6f';'#ff7f00';'#cab2d6';'#6a3d9a';'#ffff99';'#b15928'];

blueLight = pairedPalette(1,:);
blueDark = pairedPalette(2,:);
greenLight = pairedPalette(3,:);
greenDark = pairedPalette(4,:);
redLight = pairedPalette(5,:);
redDark = pairedPalette(6,:);
orangeLight = pairedPalette(7,:);
orangeDark = pairedPalette(8,:);
purpleLight = pairedPalette(9,:);
purpleDark = pairedPalette(10,:);
yellowLight = pairedPalette(11,:);
yellowDark = pairedPalette(12,:);

%% Predefine sequential color palettes

% Orange

orangePalette = ['#fff5eb';'#fee6ce';'#fdd0a2';'#fdae6b';'#fd8d3c';'#f16913';'#d94801';'#a63603';'#7f2704'];

orange1 = orangePalette(1,:);
orange2 = orangePalette(2,:);
orange3 = orangePalette(3,:);
orange4 = orangePalette(4,:);
orange5 = orangePalette(5,:);
orange6 = orangePalette(6,:);
orange7 = orangePalette(7,:);
orange8 = orangePalette(8,:);

% Blue

bluePalette = ['#f7fbff';'#deebf7';'#c6dbef';'#9ecae1';'#6baed6';'#4292c6';'#2171b5';'#08519c';'#08306b'];
bluePalette_rgb=[247 251 255; 222 235 247; 198 219 239; 158 202 225; 107 174 214; 66 146 198; 33 113 181;8 81 156; 8 48 107]/255;

bluetoorange=[16 91 143; 70 114 159; 107 138 175; 141 162 191; 174 188 208; 207 214 224; 241 241 241; 243 222 208; 242 203 176; 240 185 144; 235 167 113; 229 149 82; 222 131 49];



blue1 = bluePalette(1,:);
blue2 = bluePalette(2,:);
blue3 = bluePalette(3,:);
blue4 = bluePalette(4,:);
blue5 = bluePalette(5,:);
blue6 = bluePalette(6,:);
blue7 = bluePalette(7,:);
blue8 = bluePalette(8,:);

% Purple

purplePalette = ['#fcfbfd';'#efedf5';'#dadaeb';'#bcbddc';'#9e9ac8';'#807dba';'#6a51a3';'#54278f';'#3f007d'];

purple1 = purplePalette(1,:);
purple2 = purplePalette(2,:);
purple3 = purplePalette(3,:);
purple4 = purplePalette(4,:);
purple5 = purplePalette(5,:);
purple6 = purplePalette(6,:);
purple7 = purplePalette(7,:);
purple8 = purplePalette(8,:);

% Gray

grayPalette = ['#ffffff';'#f0f0f0';'#d9d9d9';'#bdbdbd';'#969696';'#737373';'#525252';'#252525';'#000000'];

gray1 = grayPalette(1,:);
gray2 = grayPalette(2,:);
gray3 = grayPalette(3,:);
gray4 = grayPalette(4,:);
gray5 = grayPalette(5,:);
gray6 = grayPalette(6,:);
gray7 = grayPalette(7,:);
gray8 = grayPalette(8,:);



% % Define the matrix E
% E = randi(100, 10, 10); % Example: a 10x10 matrix with random integers between 1 and 100



%% Using imagesc
% % Create a figure and axes
% figure;
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


%% Using heatmap
% clear h
% % Create the heatmap
% h = heatmap(E);

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

%% Using heatmap (plain one)
clear h
% Create the heatmap
h = heatmap(E1);
% Customize the heatmap
h.Title = 'Firm distribution conditional on Productivity level a';
h.XLabel = 'Worker type';
h.YLabel = 'Manager type';

% Flip the vertical axis
h.NodeChildren(3).YDir = 'normal';

% Remove cell labels
h.CellLabelColor = 'none';

% Remove grid lines
h.GridVisible = 'off';

% Customize font size and font name
h.FontSize = 20;
h.FontName = 'Helvetica';

% Suppress the warning temporarily
warning('off', 'MATLAB:structOnObject')

% Access the title and change its font weight
hs = struct(h);
hs.Axes.Title.FontWeight = 'normal';

% Reduce the size of the tick labels
hs.Axes.XAxis.FontSize = 16; % Adjust the font size as needed
hs.Axes.YAxis.FontSize = 16; % Adjust the font size as needed

% Access the colorbar and remove ticks and font size
hs.Colorbar.TickLength = 0;
hs.Colorbar.FontSize = 12;
minValue = min(E(:));
maxValue = max(E(:));
hs.Colorbar.Ticks = linspace(minValue, maxValue, 5); % Adjust the number of ticks as needed

% Re-enable the warning
warning('on', 'MATLAB:structOnObject')

% Apply the custom colormap
colormap(bluePalette_rgb);

% Print figure  
print('-dpdf', 'basicE1.pdf');


%% Using heatmap (plain one) for each a type as a function
clear h1 h2
h1=heatmap_plain(E1, bluePalette_rgb, 'Firm distribution conditional on Productivity level a', 16, 20,'basicE1.pdf');
h2=heatmap_plain(E2, bluePalette_rgb, 'Firm distribution conditional on Productivity level a', 16, 20,'basicE2.pdf');

% In a loop for each a type
for a=1:ats
    clear h
    e=E(:,:,a);
    h=heatmap_plain(e, bluePalette_rgb, 'Firm distribution conditional on Productivity level a='+string(a)+'', 16, 20,'basicE_'+string(a)+'.pdf');
end
    
eplus_mdist

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



%% Combined heatmaps
% clear h1 h2
% % Determine the common color limits for both heatmaps
% minValue = min(min(E1(:)), min(E2(:)));
% maxValue = max(max(E1(:)), max(E2(:)));

% % Create a tiled layout for better control over the layout
% t = tiledlayout(1, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');

% % Plot the first heatmap
% nexttile;
% h1 = heatmap(E1);
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
% h2 = heatmap(E2);
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
% colormap(bluePalette_rgb);

% % Add a centralized xlabel for the entire layout
% xlabel(t, 'Worker type', 'FontSize', 16, 'FontName', 'Helvetica');

% % Print the figure
% print('-dpdf', 'sideBySideHeatmaps.pdf');




