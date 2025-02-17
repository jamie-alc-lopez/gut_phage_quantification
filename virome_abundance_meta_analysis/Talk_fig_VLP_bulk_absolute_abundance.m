%% This script generates the bulk vs. VLP abundance comparisons figure used
% one of our talks about this work. These figures have different
% colorschemes than those in the figure

clear;clc

%Load data
dataset_paths = {'../cohort_analyses/process_shkoporov_2019/shkoporov_2019_absolute_data.mat'};
dataset_names = {'Adult'};
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
    bulk_only_virulence{i} =  data_structs{i}.data.virscore(bulk_only_ind{i});
    bulk_only_rel{i} = bulk_rel{i}(bulk_only_ind{i});
    bulk_only_subjects{i} = data_structs{i}.data.subject(bulk_only_ind{i});
    bulk_only_g = findgroups(bulk_only_subjects{i});
    total_bulk_only{i} = splitapply(@sum,bulk_only{i},bulk_only_g);
    total_bulk_only_rel{i} = splitapply(@sum,bulk_only_rel{i},bulk_only_g);
    mean_bulk_only(i) = mean(total_bulk_only{i});
    mean_bulk_only_rel(i) = mean(total_bulk_only_rel{i});

    %Compute VLP-only stats
    VLP_only{i} = VLP{i}(VLP_only_ind{i});
    VLP_only_virulence{i} =  data_structs{i}.data.virscore(VLP_only_ind{i});
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


%% Plot absolute abundance scatter and histograms with different colors 
% and content
close all

fig_type = 'gray_jitter';
VT_colors = 0.7*[125 125 125; 125 125 125]/255;

%fig_type = 'color_jitter';
%VT_colors = [230 97 0; 93 58 155]/255;
%VT_colors = [40 160 40; 75 0 146]/255;

%fig_type = 'empty';
%VT_colors = [1 1 1; 1 1 1];

%fig_type = 'no_lines';
%VT_colors = [1 1 1; 1 1 1];

dataset_markers = {'o','o'};
VT_names = {'Virulent','Temperate'};

ref_color = [0.5 0.5 0.5];
alpha = 0.25;
xlim_vec = [1e3,1e12];
ylim_vec = [1e3,1e12];
side_offsets = [-1,1];
sub_lim_vec = [0+side_offsets(1),1+side_offsets(2)];
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
        shared_mat{i}(shared_virulence{i}>=0.5,2),dataset_markers{i},...
        'MarkerEdgeColor',VT_colors(1,:),'MarkerFaceColor',...
        VT_colors(1,:),'MarkerEdgeAlpha',alpha,...
        'MarkerFaceAlpha',alpha,'DisplayName',VT_names{1});
    g2(i) = scatter(shared_mat{i}(shared_virulence{i}<0.5,1),...
        shared_mat{i}(shared_virulence{i}<0.5,2),dataset_markers{i},...
        'MarkerEdgeColor',VT_colors(2,:),'MarkerFaceColor',...
        VT_colors(2,:),'MarkerEdgeAlpha',alpha,...
        'MarkerFaceAlpha',alpha,'DisplayName',VT_names{2});

    axes(bulk_only_ax)
    hold on
    f1(i) = scatter(log10(bulk_only{i}(bulk_only_virulence{i}>=0.5)),...
        side_offsets(1)+rand(sum(bulk_only_virulence{i}>=0.5),1),dataset_markers{i},...
        'MarkerEdgeColor',VT_colors(1,:),'MarkerFaceColor',...
        VT_colors(1,:),'MarkerEdgeAlpha',alpha,'MarkerFaceAlpha',alpha);
    f2(i) = scatter(log10(bulk_only{i}(bulk_only_virulence{i}<0.5)),...
        side_offsets(2)+rand(sum(bulk_only_virulence{i}<0.5),1),dataset_markers{i},...
        'MarkerEdgeColor',VT_colors(2,:),'MarkerFaceColor',...
        VT_colors(2,:),'MarkerEdgeAlpha',alpha,'MarkerFaceAlpha',alpha);

    axes(VLP_only_ax)
    hold on
    f3(i) = scatter(side_offsets(1)+rand(sum(VLP_only_virulence{i}>=0.5),1),...
        log10(VLP_only{i}(VLP_only_virulence{i}>=0.5)),dataset_markers{i},...
        'MarkerEdgeColor',VT_colors(1,:),'MarkerFaceColor',...
        VT_colors(1,:),'MarkerEdgeAlpha',alpha,'MarkerFaceAlpha',alpha);
    f4(i) = scatter(side_offsets(2)+rand(sum(VLP_only_virulence{i}<0.5),1),...
        log10(VLP_only{i}(VLP_only_virulence{i}<0.5)),dataset_markers{i},...
        'MarkerEdgeColor',VT_colors(2,:),'MarkerFaceColor',...
        VT_colors(2,:),'MarkerEdgeAlpha',alpha,'MarkerFaceAlpha',alpha);

end

%Alter bulk axes
axes(bulk_only_ax)
xticks([])
yticks([])
ylim(sub_lim_vec)
xlim(log10(xlim_vec));
box on
text(4,sub_lim_vec(1)+0.5*range(sub_lim_vec),'Stool only','FontSize',GenFontSize,...
    'VerticalAlignment','bottom','Interpreter','none','FontName',FontName)

%Alter VLP axes
axes(VLP_only_ax)
xticks([])
yticks([])
xlim(sub_lim_vec)
ylim(log10(ylim_vec));
box on
text(sub_lim_vec(1)+0.25*range(sub_lim_vec),10.5,{'VLP', 'only'},'FontSize',GenFontSize,...
    'VerticalAlignment','bottom','Interpreter','none','FontName',FontName)


%Make the reference lines for the shared scatter plot
axes(shared_ax)
box on
ref_vec = [1e1,1e11];

if ~strcmp(fig_type,'no_lines')
    v(1) = plot(10000*ref_vec,ref_vec,'-.','Color',ref_color,...
        'LineWidth',LineWidth);
    v(2) = plot(1000*ref_vec,ref_vec,'-.','Color',ref_color,...
        'LineWidth',LineWidth);
    v(3) = plot(100*ref_vec,ref_vec,'-.','Color',ref_color,...
        'LineWidth',LineWidth);
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
    for k = 1:length(v)
        set( get( get( v(k), 'Annotation'), 'LegendInformation' ),...
            'IconDisplayStyle', 'off' );
    end

end


h(1) = plot([n_microbe, n_microbe],ylim_vec,'--','Color',ref_color,...
    'LineWidth',LineWidth);
h(2) = plot(xlim_vec,[n_microbe, n_microbe],'--','Color',ref_color,...
    'LineWidth',LineWidth);

%Make reference text
%text(1e4,1e9,{'Shared between', 'VLP and stool'},'FontSize',GenFontSize,...
%    'VerticalAlignment','bottom','Interpreter','tex','FontName',FontName)
text(1e7,n_microbe,'Microbial density','FontSize',GenFontSize,...
    'VerticalAlignment','bottom','Interpreter','tex','Color',ref_color,...
    'FontName',FontName)

%Disable legend representation of ref lines
for k = 1:length(h)
    set( get( get( h(k), 'Annotation'), 'LegendInformation' ),...
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

if contains(fig_type,'color')
    legend('Location','west','FontSize',GenFontSize,'Interpreter','none',...
        'FontName',FontName)
end

name_vector = ['plots/talk_plots/talk_fig_VLP_bulk_absolute_abundance_',...
    fig_type,'.pdf'];
exportgraphics(gcf,name_vector,'ContentType','vector')

