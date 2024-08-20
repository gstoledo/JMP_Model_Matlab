clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Predefine qualitative color palettes

bluePalette_rgb=[247 251 255; 222 235 247; 198 219 239; 158 202 225; 107 174 214; 66 146 198; 33 113 181;8 81 156; 8 48 107]/255;
bluetoorange=[16 91 143; 70 114 159; 107 138 175; 141 162 191; 174 188 208; 207 214 224; 241 241 241; 243 222 208; 242 203 176; 240 185 144; 235 167 113; 229 149 82; 222 131 49];
bluepurplePalette_rgb=[247 252 253; 224 236 244; 191 211 230; 158 188 218; 140 150 198; 140 107 177; 136 65 157; 129 15 124; 77 0 75]/255;


%% Predefine qualitative color palettes


% Dark colors

darkPalette = ['#1b9e77';'#d95f02';'#7570b3';'#e7298a';'#66a61e';'#e6ab02';'#a6761d';'#666666';'#FF00FF'];

greenColor = darkPalette(1,:);
orangeColor = darkPalette(2,:);
purpleColor = darkPalette(3,:);
pinkColor = darkPalette(4,:);
appleColor = darkPalette(5,:);
yellowColor = darkPalette(6,:);
brownColor = darkPalette(7,:);
grayColor = darkPalette(8,:);
magentaColor = darkPalette(9,:);

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

%Qualitative color palettes


qualitative=['#e41a1c';'#377eb8';'#4daf4a';'#984ea3';'#ff7f00';'#ffff33';'#a65628';'#f781bf';'#999999'];
red_qual = qualitative(1,:);
blue_qual = qualitative(2,:);
green_qual = qualitative(3,:);
purple_qual = qualitative(4,:);
orange_qual = qualitative(5,:);
yellow_qual = qualitative(6,:);
brown_qual = qualitative(7,:);
pink_qual = qualitative(8,:);
gray_qual = qualitative(9,:);




save('color_style.mat')