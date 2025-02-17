%% This script generates the non-MDA bulk vs. VLP abundance comparisons 
% found in Fig. S3. For pVMR analyses, this imputes a mean VLP density of
% 2 x 10^9 VLP/g for all VLP samples. 

clear;clc

%Load data
dataset_paths = {'../cohort_analyses/process_garmaeva_2024/garmaeva_2024_adult_absolute_data.mat',...
    '../cohort_analyses/process_garmaeva_2024/garmaeva_2024_infant_absolute_data.mat'};
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


%For convenience, print out the key overlap results
df = @(x) [num2str(round(x(1),2)),'/',num2str(round(x(2),2))];
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


%% Make the relative abundance scatter plot SI Figure

load('../common_util/font_config.mat');
data_colors = [100 143 255; 254 97 0]/255;
ref_color = [0.5 0.5 0.5];
alpha = 0.1;
xlim_vec = [1e-7,1];
ylim_vec = [1e-7,1];
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
ylabel('VLP phage genome rel. abundace','Interpreter','none',...
    'FontSize',LabelFontSize,'FontName',FontName)
set(gca,'FontSize',GenFontSize,'FontName',FontName)
legend('Location','west','FontSize',GenFontSize,'Interpreter','none',...
    'FontName',FontName)

name_vector = ['plots/SI_fig_no_MDA_VLP_bulk_relative_abundance.pdf'];
exportgraphics(gcf,name_vector,'ContentType','vector')

