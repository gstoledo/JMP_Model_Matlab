clear all
close all
location="local";
% location="hpc";
cd('/Users/gabrieltoledo/Library/CloudStorage/GoogleDrive-gabrielstoledo.gt@gmail.com/My Drive/PHD NYU/Labor Firm Structure/JMP_Model_Matlab')
addpath('/Users/gabrieltoledo/Library/CloudStorage/GoogleDrive-gabrielstoledo.gt@gmail.com/My Drive/PHD NYU/Labor Firm Structure/JMP_Model_Matlab')
addpath('/Users/gabrieltoledo/Library/CloudStorage/GoogleDrive-gabrielstoledo.gt@gmail.com/My Drive/PHD NYU/Labor Firm Structure/Python_Firm Structure/replication_HLMP/matlab/run_A4_revisit')
% load('model_solutions.mat')
load('model_solutions_sp1hpc.mat')
load('xmin_results_combined.mat')


%Check sum of the eplus final distributions
nplus=sum(eplus_edist)+sum(eplus_mdist,"all")+sum(eplus_ndist,"all")+sum(eplus_tdist,"all");
popplus=sum(eplus_udist)+sum(eplus_mdist,"all")+sum(eplus_ndist,"all")+2*sum(eplus_tdist,"all");


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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Predefine qualitative color palettes
 bluePalette_rgb=[247 251 255; 222 235 247; 198 219 239; 158 202 225; 107 174 214; 66 146 198; 33 113 181;8 81 156; 8 48 107]/255;
 bluetoorange=[16 91 143; 70 114 159; 107 138 175; 141 162 191; 174 188 208; 207 214 224; 241 241 241; 243 222 208; 242 203 176; 240 185 144; 235 167 113; 229 149 82; 222 131 49];
 bluepurplePalette_rgb=[247 252 253; 224 236 244; 191 211 230; 158 188 218; 140 150 198; 140 107 177; 136 65 157; 129 15 124; 77 0 75]/255;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Define a tiled layout with 3 rows and 2 columns
% Distribution of firms by type a 
e_a=zeros(1,ats);
for i=1:ats
    e_a(i)=eplus_edist(1,i) + sum(eplus_mdist(i,:),'all') + sum(eplus_ndist(i,:),'all')+ sum(eplus_tdist(i,:,:),'all');
end 

%PLot e_a 
bar(e_a, 'FaceColor', [0, 0.4470, 0.7410]); % MATLAB default blue
% Customize appearance
ax = gca;
ax.FontName = 'Helvetica';            % Set font to Helvetica
ax.FontSize = 14;                     % Increase font size
ax.FontWeight = 'bold';               % Set font weight to bold
ax.XGrid = 'off';                      % Add grid lines for better readability
ax.YGrid = 'off';
ax.GridLineStyle = '--';              % Dashed grid lines for a softer look

% Set labels and title
xlabel('Type a', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Density', 'FontSize', 16, 'FontWeight', 'bold');
title('Distribution of firms by type a', 'FontSize', 18, 'FontWeight', 'bold');
% Alternatively, for higher resolution:
exportgraphics(gcf, 'Figures_dist_uni/figures_dist_a.png', 'Resolution', 300);


% Distribution of talen in manager positions
e_m=zeros(1,tpts);
for i=1:tpts
    e_m(i)=sum(eplus_mdist(:,i),'all') + sum(eplus_tdist(:,i,:),'all');
end

%PLot e_m
bar(e_m, 'FaceColor', [0.8500, 0.3250, 0.0980]); % MATLAB default orange
% Customize appearance
ax = gca;
ax.FontName = 'Helvetica';            % Set font to Helvetica
ax.FontSize = 14;                     % Increase font size
ax.FontWeight = 'bold';               % Set font weight to bold
ax.XGrid = 'off';                      % Add grid lines for better readability
ax.YGrid = 'off';
ax.GridLineStyle = '--';              % Dashed grid lines for a softer look

% Set labels and title
xlabel('Type m', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Density', 'FontSize', 16, 'FontWeight', 'bold');
title('Distribution of talent in manager positions', 'FontSize', 18, 'FontWeight', 'bold');
% Alternatively, for higher resolution:
exportgraphics(gcf, 'Figures_dist_uni/figures_dist_m.png', 'Resolution', 300);


% Distribution of talent in non-manager positions
e_n=zeros(1,tpts);
for i=1:tpts
    e_n(i)=sum(eplus_ndist(:,i),'all') + sum(eplus_tdist(:,:,i),'all');
end

%PLot e_n
bar(e_n, 'FaceColor', [0.9290, 0.6940, 0.1250]); % MATLAB default yellow
% Customize appearance
ax = gca;
ax.FontName = 'Helvetica';            % Set font to Helvetica
ax.FontSize = 14;                     % Increase font size
ax.FontWeight = 'bold';               % Set font weight to bold
ax.XGrid = 'off';                      % Add grid lines for better readability
ax.YGrid = 'off';
ax.GridLineStyle = '--';              % Dashed grid lines for a softer look

% Set labels and title
xlabel('Type n', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Density', 'FontSize', 16, 'FontWeight', 'bold');
title('Distribution of talent in non-manager positions', 'FontSize', 18, 'FontWeight', 'bold');
% Alternatively, for higher resolution:
exportgraphics(gcf, 'Figures_dist_uni/figures_dist_n.png', 'Resolution', 300);



%Total final distribution of talent
typefinal=zeros(1,tpts);
for i=1:tpts
    typefinal(i)=sum(eplus_mdist(:,i),'all') + sum(eplus_ndist(:,i),'all') + sum(eplus_tdist(:,i,:),'all')+ sum(eplus_tdist(:,:,i),'all');
end
figure
bar(typefinal, 'FaceColor', [0.4660, 0.6740, 0.1880]); % MATLAB default green
% Customize appearance
ax = gca;
ax.FontName = 'Helvetica';            % Set font to Helvetica
ax.FontSize = 14;                     % Increase font size
ax.FontWeight = 'bold';               % Set font weight to bold
ax.XGrid = 'off';                      % Add grid lines for better readability
ax.YGrid = 'off';
ax.GridLineStyle = '--';              % Dashed grid lines for a softer look

% Set labels and title
xlabel('Type', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Density', 'FontSize', 16, 'FontWeight', 'bold');
title('Total final distribution of talent', 'FontSize', 18, 'FontWeight', 'bold');
% Alternatively, for higher resolution:
exportgraphics(gcf, 'Figures_dist_uni/figures_dist_final.png', 'Resolution', 300);



%To compare the "neutral"distributions, the one from birth 
%Need to ajdust this later, pay attention to normalization to measure of workers 
mnew	                =x(9);
mnew_high=tpts; %Number of types
typebirth=zeros(1,tpts);
for i=1:mnew_high
    typebirth(i)= mnew*exp(-type(i)*mnew)./sum(mnew*exp(-type(1:mnew_high)*mnew)); %Birth rate
end
bar(typebirth,'FaceColor', [0.4940, 0.1840, 0.5560]); % MATLAB default purple
% Customize appearance
ax = gca;
ax.FontName = 'Helvetica';            % Set font to Helvetica
ax.FontSize = 14;                     % Increase font size
ax.FontWeight = 'bold';               % Set font weight to bold
ax.XGrid = 'off';                      % Add grid lines for better readability
ax.YGrid = 'off';
ax.GridLineStyle = '--';              % Dashed grid lines for a softer look

% Set labels and title
xlabel('Type', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Density', 'FontSize', 16, 'FontWeight', 'bold');
title('Distribution of talent at birth', 'FontSize', 18, 'FontWeight', 'bold');
% Alternatively, for higher resolution:
exportgraphics(gcf, 'Figures_dist_uni/figures_dist_birth.png', 'Resolution', 300);



