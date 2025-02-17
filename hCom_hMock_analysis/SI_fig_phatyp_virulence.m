%% This script compares the VTR of hCom2 and the Yachida et al. cohort, as 
%measured by PhaTYP instead of BACPHLIP/Phanta. It generates the associated SI
%figure for this analysis

clear;clc

%Read data with phatyp virulence calls
hCom_manifest = ...
    wrapped_hCom_phanta_import('phanta_output/hCom_UHGV_final_merged_outputs/','phatyp');
hCom_manifest = sortrows(hCom_manifest,'mouse');
read_cutoff = 1e5;
hCom_manifest = hCom_manifest(hCom_manifest.total_reads > read_cutoff,:);

%Compute metrics of non-challenged fecal samples
in_vivo_ind = strcmp(hCom_manifest.site,'feces') & strcmp(hCom_manifest.challenged,'n');
hCom_VTR = hCom_manifest.vir_temp_ratio(in_vivo_ind);


%% Read Yachida et al. data, get PhaTYP VTR

cd('../cohort_analyses/process_yachida_2019/')
yachida_manifest= ...
    wrapped_yachida_2019_phanta_import('yachida_2019_UHGV_final_merged_outputs/','phatyp');
cd('../../hCom_hMock_analysis/')
yachida_manifest = yachida_manifest(yachida_manifest.total_reads > read_cutoff,:);

hum_VTR = yachida_manifest.vir_temp_ratio;

%% Plot the distribution of PhaTYP VTRs

load('../common_util/font_config.mat');
colors = [225 190 106; 64 176 166; 220 58 32]/255;
colors(1,:) = colors(1,:)*0.6;
colors(3,:) = colors(3,:)*0.7;
offset = 0;
MarkerSize = 12;

figure
hold on
violinplot_modified(log10(hum_VTR),1,ones(size(hum_VTR)),...
    'ViolinColor',colors(3,:),'MarkerSize',MarkerSize);
violinplot_modified(log10(hCom_VTR),2,ones(size(hCom_VTR)),...
    'ViolinColor',colors(3,:),'MarkerSize',MarkerSize);

ylim(log10([1e-1,1e1]))
yticks([-1,0,1])
yticklabels({'10^{-1}','10^0','10^1'})
set(gca,'YMinorTick','on')
xlim([0-0.5,1+0.5])
ylabel({'PhaTYP estimated','virulent to temperate ratio'},'Interpreter','tex',...
    'FontSize',LabelFontSize,'FontName',FontName)
set(gca,'FontSize',GenFontSize,'FontName',FontName)
xticks(0:1);
xticklabels({'Stool','hCom2'})
set(gca,'TickLabelInterpreter','tex')

name_vector = 'plots/SI_fig_phatyp_VTR.pdf';
exportgraphics(gcf,name_vector,'ContentType','vector')

