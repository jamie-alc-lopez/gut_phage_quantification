%% This script generates Figure 4 and associated SI figures,
% which compares the hCom microbiome to stool microbiota

%Import hCom phanta data, compute metrics

clear;clc
close all

%Read data
hCom_manifest = wrapped_hCom_phanta_import(...
    'phanta_output/hCom_UHGV_final_merged_outputs/');
hCom_manifest = sortrows(hCom_manifest,'mouse');
read_cutoff = 1e5;
hCom_manifest = hCom_manifest(hCom_manifest.total_reads > read_cutoff,:);

%Compute metrics of non-challenged fecal samples
in_vivo_ind = strcmp(hCom_manifest.site,'feces') & strcmp(hCom_manifest.challenged,'n');
hCom_bact_H = hCom_manifest.bacteria_H(in_vivo_ind);
hCom_phage_H = hCom_manifest.phage_H(in_vivo_ind);
hCom_bact_rich = hCom_manifest.bacteria_corr_richness(in_vivo_ind);
hCom_phage_rich = hCom_manifest.phage_corr_richness(in_vivo_ind);
hCom_VMR = hCom_manifest.phage_to_microbe_ratio(in_vivo_ind);
hCom_VTR = hCom_manifest.vir_temp_ratio(in_vivo_ind);

%Compute metrics of challenged fecal samples
in_vivo_chal_ind = strcmp(hCom_manifest.site,'feces') & strcmp(hCom_manifest.challenged,'y');
hCom_chal_bact_H = hCom_manifest.bacteria_H(in_vivo_chal_ind);
hCom_chal_phage_H = hCom_manifest.phage_H(in_vivo_chal_ind);
hCom_chal_bact_rich = hCom_manifest.bacteria_corr_richness(in_vivo_chal_ind);
hCom_chal_phage_rich = hCom_manifest.phage_corr_richness(in_vivo_chal_ind);
hCom_chal_VMR = hCom_manifest.phage_to_microbe_ratio(in_vivo_chal_ind);
hCom_chal_VTR = hCom_manifest.vir_temp_ratio(in_vivo_chal_ind);

%Compute metrics of ex vivo samples
ex_vivo_ind = strcmp(hCom_manifest.site,'inoc');
hCom_ex_vivo_bact_H = hCom_manifest.bacteria_H(ex_vivo_ind);
hCom_ex_vivo_phage_H = hCom_manifest.phage_H(ex_vivo_ind);
hCom_ex_vivo_bact_rich = hCom_manifest.bacteria_corr_richness(ex_vivo_ind);
hCom_ex_vivo_phage_rich = hCom_manifest.phage_corr_richness(ex_vivo_ind);
hCom_ex_vivo_VMR = hCom_manifest.phage_to_microbe_ratio(ex_vivo_ind);
hCom_ex_vivo_VTR = hCom_manifest.vir_temp_ratio(ex_vivo_ind);


%% Load in Yachida et al. 2019 Phanta data

%This initially uses the scripts in the cohort analyses folder, requires a 
%temporary change of working directory
cd('../cohort_analyses/process_yachida_2019/')
yachida_manifest= wrapped_yachida_2019_phanta_import(...
    'yachida_2019_UHGV_final_merged_outputs/');
%yachida_manifest= wrapped_yachida_2019_phanta_import(...
%    'yachida_2019_UHGV_high_cov_final_merged_outputs/');
cd('../../hCom_hMock_analysis/')
yachida_manifest = yachida_manifest(yachida_manifest.total_reads > read_cutoff,:);

%Summary statistics
hum_bact_H = yachida_manifest.bacteria_H;
hum_phage_H = yachida_manifest.phage_H;
hum_bact_rich = yachida_manifest.bacteria_corr_richness;
hum_phage_rich = yachida_manifest.phage_corr_richness;
hum_VMR = yachida_manifest.phage_to_microbe_ratio;
hum_VTR = yachida_manifest.vir_temp_ratio;

%Compute centiles and standard deviation distances of hCom wrt human data
[phage_H_cent,phage_H_d_std] = compute_relative_stats(hum_phage_H,mean(hCom_phage_H));
[phage_rich_cent,phage_rich_d_std] = compute_relative_stats(hum_phage_rich,mean(hCom_phage_rich));
[VMR_cent,VMR_d_std] = compute_relative_stats(hum_VMR,mean(hCom_VMR));
max_phage_H_cent = compute_relative_stats(hum_phage_H,max(hCom_phage_H)); 
min_phage_H_VMR_cent = compute_relative_stats(hum_phage_H,min(hCom_phage_H)); 

%% Plot the main human comparison figure
close all
figure
pos = get(gcf,'Position');
pos(3) = 1.5*pos(3);
pos(4) = 2*0.85*pos(4);
set(gcf,'Position',pos)

load('../common_util/font_config.mat');
colors = [225 190 106; 64 176 166; 220 58 32]/255;
colors(1,:) = colors(1,:)*0.6;
colors(3,:) = colors(3,:)*0.7;
offset = 0;
MarkerSize = 12;
remove_outlier = 0;

% Make figure 4A - Shannon entropy
subplot(2,2,1)
hold on
violinplot_modified(hum_bact_H,1,ones(size(hum_bact_H)),...
    'ViolinColor',colors(1,:),'MarkerSize',MarkerSize);
violinplot_modified(hCom_bact_H,2,ones(size(hCom_bact_H)),...
    'ViolinColor',colors(1,:),'MarkerSize',MarkerSize);
violinplot_modified(hum_phage_H,3+offset,ones(size(hum_phage_H)),...
    'ViolinColor',colors(2,:),'MarkerSize',MarkerSize);
violinplot_modified(hCom_phage_H,4+offset,ones(size(hCom_phage_H)),...
    'ViolinColor',colors(2,:),'MarkerSize',MarkerSize);
ylim([0,8])
xlim([0-0.5,3+0.5])
ylabel('Shannon diversity, H','Interpreter','tex','FontSize',LabelFontSize,...
    'FontName',FontName)
set(gca,'FontSize',GenFontSize,'FontName',FontName)
xticks(0:3);
xticklabels({'Stool bact.','hCom2 bact.','Stool phage','hCom2 phage'})
set(gca,'TickLabelInterpreter','tex')
text(-0.2,1.06,'A','Interpreter','tex','Units','normalized',...
    'FontSize',PanelFontSize,'FontName',FontName)


% Make figure 4B - Weighted richness
subplot(2,2,2)
hold on
violinplot_modified(hum_bact_rich,1,ones(size(hum_bact_rich)),...
    'ViolinColor',colors(1,:),'MarkerSize',MarkerSize);
violinplot_modified(hCom_bact_rich,2,ones(size(hCom_bact_rich)),...
    'ViolinColor',colors(1,:),'MarkerSize',MarkerSize);
violinplot_modified(hum_phage_rich,3+offset,ones(size(hum_phage_rich)),...
    'ViolinColor',colors(2,:),'MarkerSize',MarkerSize);
violinplot_modified(hCom_phage_rich,4+offset,ones(size(hCom_phage_rich)),...
    'ViolinColor',colors(2,:),'MarkerSize',MarkerSize);
xlim([0-0.5,3+0.5])
ylabel('Weighted richness','Interpreter','tex','FontSize',LabelFontSize,...
    'FontName',FontName)
set(gca,'FontSize',GenFontSize,'FontName',FontName)
xticks(0:3);
xticklabels({'Stool bact.','hCom2 bact.','Stool phage','hCom2 phage'})
set(gca,'TickLabelInterpreter','tex')
text(-0.2,1.06,'B','Interpreter','tex','Units','normalized',...
    'FontSize',PanelFontSize,'FontName',FontName)


% Make Figure 4C - VMR
subplot(2,2,3)
hold on
violinplot_modified(log10(hum_VMR),1,ones(size(hum_VMR)),...
    'ViolinColor',colors(3,:),'MarkerSize',MarkerSize);
violinplot_modified(log10(hCom_VMR),2,ones(size(hCom_VMR)),...
    'ViolinColor',colors(3,:),'MarkerSize',MarkerSize);
ylim(log10([5e-1,2e1]))
yticks([0,1])
yticklabels({'10^0','10^1'})
set(gca,'YScale','linear')
xlim([0-0.5,1+0.5])
ylabel('Virus to microbe ratio','Interpreter','tex','FontSize',LabelFontSize,...
    'FontName',FontName)
set(gca,'FontSize',GenFontSize,'FontName',FontName)
xticks(0:1);
xticklabels({'Stool','hCom2'})
set(gca,'YMinorTick','on')
set(gca,'TickLabelInterpreter','tex')
text(-0.2,1.06,'C','Interpreter','tex','Units','normalized',...
    'FontSize',PanelFontSize,'FontName',FontName)

% Make Figure 4D - VTR
subplot(2,2,4)
hold on
violinplot_modified(log10(hum_VTR),1,ones(size(hum_VTR)),...
    'ViolinColor',colors(3,:),'MarkerSize',MarkerSize);
violinplot_modified(log10(hCom_VTR),2,ones(size(hCom_VTR)),...
    'ViolinColor',colors(3,:),'MarkerSize',MarkerSize);
%plot([1-0.2,1+0.2],[0,0],'k--','LineWidth',2);

%ylim([-0.5,6])
ylim(log10([1e-2,3e1]))
yticks([-2,-1,0,1])
yticklabels({'10^{-2}','10^{-1}','10^0','10^1'})
set(gca,'YScale','linear')
xlim([0-0.5,1+0.5])
ylabel({'Estimated virulent', 'to temperate ratio'},'Interpreter','tex')
set(gca,'FontSize',GenFontSize,'FontName',FontName)
xticks(0:1);
xticklabels({'Stool','hCom2'})
set(gca,'YMinorTick','on')
set(gca,'TickLabelInterpreter','tex')
text(-0.2,1.06,'D','Interpreter','tex','Units','normalized',...
    'FontSize',PanelFontSize,'FontName',FontName)

name_vector = 'plots/Fig_hCom_human_comparison.pdf';
exportgraphics(gcf,name_vector,'ContentType','vector')


%% Plot the comparison between different hCom types
close all
figure
pos = get(gcf,'Position');
pos(3) = 1.5*pos(3);
pos(4) = 2*0.85*pos(4);
set(gcf,'Position',pos)

colors = [225 190 106; 64 176 166; 220 58 32]/255;
colors(1,:) = colors(1,:)*0.6;
colors(3,:) = colors(3,:)*0.7;
offset = 0;
MarkerSize = 12;


% Shannon entropy
subplot(2,2,1)
hold on
violinplot_modified(hCom_bact_H,1,ones(size(hCom_bact_H)),...
    'ViolinColor',colors(1,:),'MarkerSize',MarkerSize);
violinplot_modified(hCom_chal_bact_H,2,ones(size(hCom_chal_bact_H)),...
    'ViolinColor',colors(1,:),'MarkerSize',MarkerSize);
violinplot_modified(hCom_ex_vivo_bact_H,3,ones(size(hCom_ex_vivo_bact_H)),...
    'ViolinColor',colors(1,:),'MarkerSize',MarkerSize);

violinplot_modified(hCom_phage_H,4+offset,ones(size(hCom_bact_H)),...
    'ViolinColor',colors(2,:),'MarkerSize',MarkerSize);
violinplot_modified(hCom_chal_phage_H,5+offset,ones(size(hCom_chal_bact_H)),...
    'ViolinColor',colors(2,:),'MarkerSize',MarkerSize);
violinplot_modified(hCom_ex_vivo_phage_H,6+offset,ones(size(hCom_ex_vivo_bact_H)),...
    'ViolinColor',colors(2,:),'MarkerSize',MarkerSize);

ylim([0,8])
xlim([0-0.5,5+0.5])
ylabel('Shannon diversity, H','Interpreter','tex',...
    'FontSize',LabelFontSize,'FontName',FontName)
set(gca,'FontSize',GenFontSize,'FontName',FontName)
xticks(0:6);
xticklabels({'Unchallenged bact.','Challenged bact.','in vitro bact.',...
    'Unchallenged phage','Challenged phage','in vitro phage'})
set(gca,'TickLabelInterpreter','tex')
text(-0.2,1.06,'A','Interpreter','tex','Units','normalized',...
    'FontSize',PanelFontSize,'FontName',FontName)


% Weighted richness
subplot(2,2,2)
hold on
violinplot_modified(hCom_bact_rich,1,ones(size(hCom_bact_H)),...
    'ViolinColor',colors(1,:),'MarkerSize',MarkerSize);
violinplot_modified(hCom_chal_bact_rich,2,ones(size(hCom_chal_bact_H)),...
    'ViolinColor',colors(1,:),'MarkerSize',MarkerSize);
violinplot_modified(hCom_ex_vivo_bact_rich,3,ones(size(hCom_ex_vivo_bact_H)),...
    'ViolinColor',colors(1,:),'MarkerSize',MarkerSize);

violinplot_modified(hCom_phage_rich,4+offset,ones(size(hCom_bact_H)),...
    'ViolinColor',colors(2,:),'MarkerSize',MarkerSize);
violinplot_modified(hCom_chal_phage_rich,5+offset,ones(size(hCom_chal_bact_H)),...
    'ViolinColor',colors(2,:),'MarkerSize',MarkerSize);
violinplot_modified(hCom_ex_vivo_phage_rich,6+offset,ones(size(hCom_ex_vivo_bact_H)),...
    'ViolinColor',colors(2,:),'MarkerSize',MarkerSize);

ylim([0,200])
xlim([0-0.5,5+0.5])
ylabel('Weighted richness','Interpreter','tex',...
    'FontSize',LabelFontSize,'FontName',FontName)
set(gca,'FontSize',GenFontSize,'FontName',FontName)
xticks(0:5);
xticklabels({'Unchallenged bact.','Challenged bact.','in vitro bact.',...
    'Unchallenged phage','Challenged phage','in vitro phage'})
set(gca,'TickLabelInterpreter','tex')
text(-0.2,1.06,'B','Interpreter','tex','Units','normalized',...
    'FontSize',PanelFontSize,'FontName',FontName)


% VMR
subplot(2,2,3)
hold on
violinplot_modified(log10(hCom_VMR),1,ones(size(hCom_bact_H)),...
    'ViolinColor',colors(3,:),'MarkerSize',MarkerSize);
violinplot_modified(log10(hCom_chal_VMR),2,ones(size(hCom_chal_bact_H)),...
    'ViolinColor',colors(3,:),'MarkerSize',MarkerSize);
violinplot_modified(log10(hCom_ex_vivo_VMR),3,ones(size(hCom_ex_vivo_bact_H)),...
    'ViolinColor',colors(3,:),'MarkerSize',MarkerSize);
ylim(log10([5e-1,2e1]))
yticks([0,1])
yticklabels({'10^0','10^1'})
set(gca,'YMinorTick','on')
xlim([0-0.5,2+0.5])
ylabel('Virus to microbe ratio','Interpreter','tex',...
    'FontSize',LabelFontSize,'FontName',FontName)
set(gca,'FontSize',GenFontSize,'FontName',FontName)
xticks(0:2);
xticklabels({'Unchallenged','Challenged','in vitro'})
set(gca,'TickLabelInterpreter','tex')
text(-0.2,1.06,'C','Interpreter','tex','Units','normalized',...
    'FontSize',PanelFontSize,'FontName',FontName)

% VTR
subplot(2,2,4)
hold on
violinplot_modified(log10(hCom_VTR),1,ones(size(hCom_bact_H)),...
    'ViolinColor',colors(3,:),'MarkerSize',MarkerSize);
violinplot_modified(log10(hCom_chal_VTR),2,ones(size(hCom_chal_bact_H)),...
    'ViolinColor',colors(3,:),'MarkerSize',MarkerSize);
violinplot_modified(log10(hCom_ex_vivo_VTR),3,ones(size(hCom_ex_vivo_bact_H)),...
    'ViolinColor',colors(3,:),'MarkerSize',MarkerSize);

ylim(log10([1e-1,3e1]))
yticks([-1,0,1])
yticklabels({'10^{-1}','10^0','10^1'})
set(gca,'YMinorTick','on')
xlim([0-0.5,2+0.5])
ylabel({'Estimated virulent','to temperate ratio'},'Interpreter','tex',...
    'FontSize',LabelFontSize,'FontName',FontName)
set(gca,'FontSize',GenFontSize,'FontName',FontName)
xticks(0:2);
xticklabels({'Unchallenged','Challenged','in vitro'})
set(gca,'TickLabelInterpreter','tex')
text(-0.2,1.06,'D','Interpreter','tex','Units','normalized',...
    'FontSize',PanelFontSize,'FontName',FontName)


name_vector='plots/SI_fig_hCom_internal_comparison.pdf';
exportgraphics(gcf,name_vector,'ContentType','vector')



