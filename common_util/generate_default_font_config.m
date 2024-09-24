%This script generate the default font configuration file. This config file
%is used by the figure plotting scripts

clear;clc

LabelFontSize = 22; %Font size for x and y labels
GenFontSize = 18; %Font size for legend, ticks, etc
PanelFontSize = 26; %Font size for panel label letters
FontName = 'Arial'; %Font to be used

save('font_config.mat')