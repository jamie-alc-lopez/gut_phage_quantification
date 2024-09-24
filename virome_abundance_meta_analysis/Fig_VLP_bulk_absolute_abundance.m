%% This script generates the bulk vs. VLP abundance comparisons found in 
%Fig 1CD and relevant supplementary figures

clear;clc

%Load data
dataset_paths = {'../cohort_analyses/process_shkoporov_2019/shkoporov_2019_absolute_data.mat',...
    '../cohort_analyses/process_liang_2020/liang_2020_absolute_data.mat'};
dataset_names = {'Adult','Infant'};
for i = 1:length(dataset_paths)
    data_structs{i} = load(dataset_paths{i});
end

%Loop through datasets and process absolute abundances
for i = 1:length(data_structs)
    bulk{i} = data_structs{i}.data.stool_abs;
    VLP{i} = data_structs{i}.data.VLP_abs;

    bulk_rel{i} = data_structs{i}.data.stool_rel_phage;
    VLP_rel{i} = data_structs{i}.data.VLP_rel_phage;

    n_subjects(i) = length(unique(data_structs{i}.data.subject));

    %Identify which entries are shared, bulk only, and VLP only
    shared_ind{i} = (bulk{i} > 0) & (VLP{i} > 0);
    bulk_only_ind{i} = (bulk{i} > 0) & (VLP{i} == 0);
    VLP_only_ind{i} = (bulk{i} == 0) & (VLP{i} > 0);

    %Compute shared stats
    shared_mat{i} = [bulk{i}(shared_ind{i}),VLP{i}(shared_ind{i})];
    shared_rel_mat{i} = [bulk_rel{i}(shared_ind{i}),VLP_rel{i}(shared_ind{i})];
    shared_subjects{i} = data_structs{i}.data.subject(shared_ind{i});
    shared_virulence{i} = data_structs{i}.data.virscore(shared_ind{i});
    shared_g = findgroups(shared_subjects{i});
    total_shared_bulk{i} = splitapply(@sum,shared_mat{i}(:,1),shared_g);
    mean_shared_bulk(i) = mean(total_shared_bulk{i});
    total_shared_bulk_rel{i} = splitapply(@sum,shared_rel_mat{i}(:,1),shared_g);
    mean_shared_bulk_rel(i) = mean(total_shared_bulk_rel{i});

    total_shared_VLP{i} = splitapply(@sum,shared_mat{i}(:,2),shared_g);
    mean_shared_VLP(i) = mean(total_shared_VLP{i});
    total_shared_VLP_rel{i} = splitapply(@sum,shared_rel_mat{i}(:,2),shared_g);
    mean_shared_VLP_rel(i) = mean(total_shared_VLP_rel{i});
    frac_VMR_under_1_to_100(i) = sum((shared_mat{i}(:,2)./shared_mat{i}(:,1))<0.01)/length(shared_mat{i}(:,2));

    %Compute shared spearman correlation, to use splitapply need to employ
    %a slightly odd cell2mat constructions
    shared_mat_cell{i} = mat2cell(shared_mat{i},ones(size(shared_mat{i},1),1));
    total_shared_spearman{i} = splitapply(@(x) {corr(cell2mat(x),'type','Spearman')},shared_mat_cell{i},shared_g);
    total_shared_spearman{i} = cellfun(@(x) x(2,1),total_shared_spearman{i});
    mean_shared_spearman(i) = mean(total_shared_spearman{i});

    %Compute bulk-only stats
    bulk_only{i} = bulk{i}(bulk_only_ind{i});
    bulk_only_rel{i} = bulk_rel{i}(bulk_only_ind{i});
    bulk_only_subjects{i} = data_structs{i}.data.subject(bulk_only_ind{i});
    bulk_only_g = findgroups(bulk_only_subjects{i});
    total_bulk_only{i} = splitapply(@sum,bulk_only{i},bulk_only_g);
    total_bulk_only_rel{i} = splitapply(@sum,bulk_only_rel{i},bulk_only_g);
    mean_bulk_only(i) = mean(total_bulk_only{i});
    mean_bulk_only_rel(i) = mean(total_bulk_only_rel{i});

    %Compute VLP-only stats
    VLP_only{i} = VLP{i}(VLP_only_ind{i});
    VLP_only_rel{i} = VLP_rel{i}(VLP_only_ind{i});
    VLP_only_subjects{i} = data_structs{i}.data.subject(VLP_only_ind{i});
    VLP_only_g = findgroups(VLP_only_subjects{i});
    total_VLP_only{i} = splitapply(@sum,VLP_only{i},VLP_only_g);
    total_VLP_only_rel{i} = splitapply(@sum,VLP_only_rel{i},VLP_only_g);
    mean_VLP_only(i) = mean(total_VLP_only{i});
    mean_VLP_only_rel(i) = mean(total_VLP_only_rel{i});

    %Compute per taxa stats across subjects
    bulk_only_taxa_frac(i) = sum(bulk_only_ind{i})/sum(bulk{i} > 0 | VLP{i} > 0);
    VLP_only_taxa_frac(i) = sum(VLP_only_ind{i})/sum(bulk{i} > 0 | VLP{i} > 0);
    shared_taxa_frac(i) = sum(shared_ind{i})/sum(bulk{i} > 0 | VLP{i} > 0);

end

df = @(x) [num2str(round(x(1),2)),'/',num2str(round(x(2),2))];

%For convenience, print out the key overlap results
disp(['Bulk-only phage account for on average '...
    ,df(mean_bulk_only_rel),' in adult/infant stool samples'])
disp(['VLP-only phage account for on average '...
    ,df(mean_VLP_only_rel),' in adult/infant VLP samples'])
disp(['Shared phage account for on average ',...
    df(mean_shared_bulk_rel),' in adult/infant stool samples'])
disp(['Shared phage account for on average ',...
    df(mean_shared_VLP_rel),' in adult/infant VLP samples'])
disp(['In total, ',df(frac_VMR_under_1_to_100),...
    ' of shared phage have a VMR < 0.01 in adults/infants'])

%% Plot absolute abundance scatter and histograms Fig. 1C
close all

data_colors = [100 143 255; 254 97 0]/255;
ref_color = [0.5 0.5 0.5];
alpha = 0.2;
xlim_vec = [1e3,1e12];
ylim_vec = [1e3,1e12];
sub_lim_vec = [0,0.77];
tick_vec = [1e3 1e6 1e9 1e12];
n_microbe = 0.92e11;
load('../common_util/font_config.mat');
LineWidth = 1.5;
vmr_label_angle = 45;
NumBins = 50;

%Make figure and primary axes
figure
fig = gcf;

%Alter figure size
pos = get(gcf,'Position');
pos(4) = pos(3);
set(gcf,'Position',pos)

%Alter main scatter plot position
main_xy = 0.3;
main_dim = 0.65;
VLP_xy = 0.07;
bulk_xy = 0.09;
sub_dim2 = 0.1;
shared_ax = axes('Position',[main_xy main_xy main_dim main_dim]);

%Make subaxes
bulk_only_ax = axes('Position',[main_xy,bulk_xy,main_dim,sub_dim2]);
VLP_only_ax = axes('Position',[VLP_xy,main_xy,sub_dim2,main_dim]);


for i = 1:length(data_structs)

    %Make the shared scatter plot
    axes(shared_ax)
    hold on
    g(i) = scatter(shared_mat{i}(shared_virulence{i}>=0.5,1),...
        shared_mat{i}(shared_virulence{i}>=0.5,2),'^',...
        'MarkerEdgeColor',data_colors(i,:),'MarkerFaceColor',...
        data_colors(i,:),'MarkerEdgeAlpha',alpha,'MarkerFaceAlpha',...
        alpha,'DisplayName',dataset_names{i});
    scatter(shared_mat{i}(shared_virulence{i}<0.5,1),...
        shared_mat{i}(shared_virulence{i}<0.5,2),'o','MarkerEdgeColor',...
        data_colors(i,:),'MarkerFaceColor',data_colors(i,:),...
        'MarkerEdgeAlpha',alpha,'MarkerFaceAlpha',alpha,...
        'DisplayName',dataset_names{i})
    axes(bulk_only_ax)
    hold on
    histogram(log10(bulk_only{i}),'BinLimits',log10(xlim_vec),...
        'Normalization','pdf','NumBins',NumBins,'EdgeColor',data_colors(i,:),...
        'FaceColor',data_colors(i,:),'FaceAlpha',alpha)

    axes(VLP_only_ax)
    hold on
    histogram(log10(VLP_only{i}),'BinLimits',log10(ylim_vec),...
        'Normalization','pdf','NumBins',NumBins,'EdgeColor',data_colors(i,:),...
        'FaceColor',data_colors(i,:),'FaceAlpha',alpha,'Orientation','horizontal')

end

%Alter bulk axes
axes(bulk_only_ax)
xticks([])
yticks([])
ylim(sub_lim_vec)
xlim(log10(xlim_vec));
box on
text(4,0.5*sub_lim_vec(2),'Stool only','FontSize',GenFontSize,...
    'VerticalAlignment','bottom','Interpreter','none','FontName',FontName)
ylabel('Prob.','FontSize',LabelFontSize,'Interpreter','none')

%Alter VLP axes
axes(VLP_only_ax)
xticks([])
yticks([])
xlim(sub_lim_vec)
ylim(log10(ylim_vec));
box on
text(0.25*sub_lim_vec(2),9,{'VLP', 'only'},'FontSize',GenFontSize,...
    'VerticalAlignment','bottom','Interpreter','none','FontName',FontName)
xlabel('Prob.','FontSize',LabelFontSize,'Interpreter','none','FontName',FontName)


%Make the reference lines for the shared scatter plot
axes(shared_ax)
box on
ref_vec = [1e1,1e11];
h(1) = plot(10000*ref_vec,ref_vec,'-.','Color',ref_color,...
    'LineWidth',LineWidth);
h(2) = plot(1000*ref_vec,ref_vec,'-.','Color',ref_color,...
    'LineWidth',LineWidth);
h(3) = plot(100*ref_vec,ref_vec,'-.','Color',ref_color,...
    'LineWidth',LineWidth);
h(4) = plot([n_microbe, n_microbe],ylim_vec,'--','Color',ref_color,...
    'LineWidth',LineWidth);
h(5) = plot(xlim_vec,[n_microbe, n_microbe],'--','Color',ref_color,...
    'LineWidth',LineWidth);

%Make reference text
text(1e4,1e9,{'Shared between', 'VLP and stool'},'FontSize',GenFontSize,...
    'VerticalAlignment','bottom','Interpreter','tex','FontName',FontName)
text(1e7,n_microbe,'Microbial density','FontSize',GenFontSize,...
    'VerticalAlignment','bottom','Interpreter','tex','Color',ref_color,...
    'FontName',FontName)
text(1e5,1e3,'1:10^2','FontSize',GenFontSize,'Rotation',vmr_label_angle,...
    'VerticalAlignment','bottom','Interpreter','tex','Color',ref_color,...
    'FontName',FontName)
text(1e6,1e3,'1:10^3','FontSize',GenFontSize,'Rotation',vmr_label_angle,...
    'VerticalAlignment','bottom','Interpreter','tex','Color',ref_color,...
    'FontName',FontName)
text(1e7,1e3,'1:10^4','FontSize',GenFontSize,'Rotation',vmr_label_angle,...
    'VerticalAlignment','bottom','Interpreter','tex','Color',ref_color,...
    'FontName',FontName)

%Disable legend representation of ref lines
for k = 1:length(h)
    set( get( get( h(k), 'Annotation'), 'LegendInformation' ),...
        'IconDisplayStyle', 'off' );
end
for k = 1:length(g)
    set( get( get( g(k), 'Annotation'), 'LegendInformation' ),...
        'IconDisplayStyle', 'off' );
end

%Alter main axis axis properties
xlim(xlim_vec)
xticks(tick_vec)
ylim(ylim_vec)
yticks(tick_vec)
set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'TickLabelInterpreter','none')
xlabel('Stool phage genome density (per g)','Interpreter','none',...
    'FontName',FontName,'FontSize',LabelFontSize)
ylabel('VLP phage genome or particle density (per g)','Interpreter','none',...
    'FontName',FontName,'FontSize',LabelFontSize)
set(gca,'FontSize',GenFontSize)
text(-0.45,1.03,'C','Interpreter','none','Units','normalized',...
    'FontSize',PanelFontSize,'FontName',FontName,'FontSize',PanelFontSize)
legend('Location','west','FontSize',GenFontSize,'Interpreter','none',...
    'FontName',FontName)

name_vector = ['plots/Fig_VLP_bulk_absolute_abundance.pdf'];
exportgraphics(gcf,name_vector,'ContentType','vector')

%% Plot distribution of particle numbers Fig. 1D

figure
MarkerSize = 12;
os = 1;
xpos_vec = [1 2 3+os 4+os 5+2*os 6+2*os 7+3*os 8+3*os];
tick_xpos_vec = [1.5 3.5+os 5.5+2*os 7.5+3*os] -1;

violinplot_modified(total_VLP_only{1},xpos_vec(1),ones(size(total_VLP_only{1})),...
    'ViolinColor',data_colors(1,:),'MarkerSize',MarkerSize);
violinplot_modified(total_VLP_only{2},xpos_vec(2),ones(size(total_VLP_only{2})),...
    'ViolinColor',data_colors(2,:),'MarkerSize',MarkerSize);

violinplot_modified(total_shared_VLP{1},xpos_vec(3),ones(size(total_shared_VLP{1})),...
    'ViolinColor',data_colors(1,:),'MarkerSize',MarkerSize);
violinplot_modified(total_shared_VLP{2},xpos_vec(4),ones(size(total_shared_VLP{2})),...
    'ViolinColor',data_colors(2,:),'MarkerSize',MarkerSize);

violinplot_modified(total_bulk_only{1},xpos_vec(5),ones(size(total_bulk_only{1})),...
    'ViolinColor',data_colors(1,:),'MarkerSize',MarkerSize);
violinplot_modified(total_bulk_only{2},xpos_vec(6),ones(size(total_bulk_only{2})),...
    'ViolinColor',data_colors(2,:),'MarkerSize',MarkerSize);

violinplot_modified(total_shared_bulk{1},xpos_vec(7),ones(size(total_shared_bulk{1})),...
    'ViolinColor',data_colors(1,:),'MarkerSize',MarkerSize);
violinplot_modified(total_shared_bulk{2},xpos_vec(8),ones(size(total_shared_bulk{2})),...
    'ViolinColor',data_colors(2,:),'MarkerSize',MarkerSize);

set(gca,'YScale','log')
ylim_vec = [1e5,1e13];
ytick_vec = [1e5 1e7 1e9 1e11 1e13];
ylim(ylim_vec);
ylabel({'Phage genome or', 'particle density (per g)'},'Interpreter','none',...
    'FontSize',LabelFontSize,'FontName',FontName)
set(gca,'FontSize',GenFontSize)

xticks(tick_xpos_vec);
xticklabels([]);
xlim([xpos_vec(1)-1 xpos_vec(end)+1]-1)

pos = get(gcf,'Position');
pos(4) = 1*pos(4);
set(gcf,'Position',pos)

pos = get(gca,'Position');
pos(4) = 0.69;
pos(2) = 0.24;
set(gca,'Position',pos)

%Make tick labels manually
ticklabels = {{'VLP','genomes', 'or particles', '(only in VLP)'},...
    {'VLP', 'genomes', 'or particles', '(in both)'},...
    {'Stool','genomes', '(only in stool)'},{'Stool', 'genomes', '(in both)'}};
for i = 1:length(ticklabels)
    text(tick_xpos_vec(i),0.05*ylim_vec(1),ticklabels{i},...
        'FontSize',GenFontSize,'Interpreter','none','FontName',FontName,...
        'HorizontalAlignment','center')
end

yticks(ytick_vec)
set(gca,'TickLabelInterpreter','tex')
text(-0.18,1.03,'D','Interpreter','none','Units','normalized',...
    'FontSize',PanelFontSize)
box off

name_vector = ['plots/Fig_VLP_bulk_frac_comparison.pdf'];
exportgraphics(gcf,name_vector,'ContentType','vector')

%% Make the virulence-colored scatter plot SI Figure

VT_colors = [40 160 40; 75 0 146]/255;
dataset_markers = {'d','sq'};
VT_names = {'Virulent','Temperate'};
ref_color = [0.5 0.5 0.5];
alpha = 0.2;
xlim_vec = [1e3,1e12];
ylim_vec = [1e3,1e12];
tick_vec = [1e3 1e6 1e9 1e12];
n_microbe = 0.92e11;
LineWidth = 1.5;

%Make figure and primary axes
figure
for i = 1:length(data_structs)

    %Make the shared scatter plot
    hold on
    h1(i) = scatter(shared_mat{i}(shared_virulence{i}>=0.5,1),...
        shared_mat{i}(shared_virulence{i}>=0.5,2),dataset_markers{i},...
        'MarkerEdgeColor',VT_colors(1,:),'MarkerFaceColor',...
        VT_colors(1,:),'MarkerEdgeAlpha',alpha,...
        'MarkerFaceAlpha',alpha,'DisplayName',VT_names{1});
    h2(i) = scatter(shared_mat{i}(shared_virulence{i}<0.5,1),...
        shared_mat{i}(shared_virulence{i}<0.5,2),dataset_markers{i},...
        'MarkerEdgeColor',VT_colors(2,:),'MarkerFaceColor',...
        VT_colors(2,:),'MarkerEdgeAlpha',alpha,...
        'MarkerFaceAlpha',alpha,'DisplayName',VT_names{2});
end

set( get( get( h1(1), 'Annotation'), 'LegendInformation' ),...
    'IconDisplayStyle', 'off' );
set( get( get( h2(1), 'Annotation'), 'LegendInformation' ),...
    'IconDisplayStyle', 'off' );

box on
ref_vec = [1e1,1e11];
h(1) = plot(10000*ref_vec,ref_vec,'-.','Color',ref_color,...
    'LineWidth',LineWidth);
h(2) = plot(1000*ref_vec,ref_vec,'-.','Color',ref_color,...
    'LineWidth',LineWidth);
h(3) = plot(100*ref_vec,ref_vec,'-.','Color',ref_color,...
    'LineWidth',LineWidth);
h(4) = plot([n_microbe, n_microbe],ylim_vec,'--','Color',ref_color,...
    'LineWidth',LineWidth);
h(5) = plot(xlim_vec,[n_microbe, n_microbe],'--','Color',ref_color,...
    'LineWidth',LineWidth);

%Make reference text
text(1e7,n_microbe,'Microbial density','FontSize',GenFontSize,...
    'VerticalAlignment','bottom','Interpreter','tex','Color',ref_color,...
    'FontName',FontName)
text(1e5,1e3,'1:10^2','FontSize',GenFontSize,'Rotation',vmr_label_angle,...
    'VerticalAlignment','bottom','Interpreter','tex','Color',ref_color,...
    'FontName',FontName)
text(1e6,1e3,'1:10^3','FontSize',GenFontSize,'Rotation',vmr_label_angle,...
    'VerticalAlignment','bottom','Interpreter','tex','Color',ref_color,...
    'FontName',FontName)
text(1e7,1e3,'1:10^4','FontSize',GenFontSize,'Rotation',vmr_label_angle,...
    'VerticalAlignment','bottom','Interpreter','tex','Color',ref_color,...
    'FontName',FontName)

%Disable legend representation of ref lines
for k = 1:length(h)
    set( get( get( h(k), 'Annotation'), 'LegendInformation' ),...
        'IconDisplayStyle', 'off' );
end

xlim(xlim_vec)
xticks(tick_vec)
ylim(ylim_vec)
yticks(tick_vec)
set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'TickLabelInterpreter','none')
xlabel('Stool phage genome density (per g)','Interpreter','none',...
    'FontSize',LabelFontSize,'FontName',FontName)
ylabel('VLP phage genome or particle density (per g)','Interpreter','none',...
    'FontSize',LabelFontSize,'FontName',FontName)
set(gca,'FontSize',GenFontSize,'FontName',FontName)
legend('Location','west','FontSize',GenFontSize,'Interpreter','none',...
    'FontName',FontName)

name_vector = ['plots/SI_fig_vir_temp_VLP_bulk_abundance.pdf'];
exportgraphics(gcf,name_vector,'ContentType','vector')



%% Make the relative abundance scatter plot SI Figure

data_colors = [100 143 255; 254 97 0]/255;
ref_color = [0.5 0.5 0.5];
alpha = 0.2;
xlim_vec = [1e-6,1];
ylim_vec = [1e-6,1];
sub_lim_vec = [0,0.77];
tick_vec = [1e-6 1e-4 1e-2 1e0];
LineWidth = 1.5;

%Make figure and primary axes
figure
for i = 1:length(data_structs)

    %Make the shared scatter plot
    hold on
    g(i) = scatter(shared_rel_mat{i}(shared_virulence{i}>=0.5,1),...
        shared_rel_mat{i}(shared_virulence{i}>=0.5,2),'^',...
        'MarkerEdgeColor',data_colors(i,:),'MarkerFaceColor',...
        data_colors(i,:),'MarkerEdgeAlpha',alpha,'MarkerFaceAlpha',...
        alpha,'DisplayName',dataset_names{i});
    scatter(shared_rel_mat{i}(shared_virulence{i}<0.5,1),...
        shared_rel_mat{i}(shared_virulence{i}<0.5,2),'o',...
        'MarkerEdgeColor',data_colors(i,:),'MarkerFaceColor',...
        data_colors(i,:),'MarkerEdgeAlpha',alpha,...
        'MarkerFaceAlpha',alpha,'DisplayName',dataset_names{i})
end

for k = 1:length(g)
    set( get( get( g(k), 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
end

xlim(xlim_vec)
xticks(tick_vec)
ylim(ylim_vec)
yticks(tick_vec)
set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'TickLabelInterpreter','none')
xlabel('Stool phage genome rel. abundance','Interpreter','none',...
    'FontSize',LabelFontSize,'FontName',FontName)
ylabel('VLP phage genome or particle rel. abundace','Interpreter','none',...
    'FontSize',LabelFontSize,'FontName',FontName)
set(gca,'FontSize',GenFontSize,'FontName',FontName)
legend('Location','west','FontSize',GenFontSize,'Interpreter','none',...
    'FontName',FontName)

name_vector = ['plots/SI_fig_VLP_bulk_relative_abundance.pdf'];
exportgraphics(gcf,name_vector,'ContentType','vector')

