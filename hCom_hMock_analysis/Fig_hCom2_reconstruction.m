%% This script generates Figure 3 and associated SI figures,
% which deal with the in silico hCom reconstruction

clear;clc
close all

%Import hMock and hCom
hCom_manifest = wrapped_hCom_phanta_import('phanta_output/hCom_UHGV_final_merged_outputs/');
hCom_manifest = sortrows(hCom_manifest,'mouse');
read_cutoff = 1e5;
hCom_manifest = hCom_manifest(hCom_manifest.total_reads > read_cutoff,:);

hMock_manifest = wrapped_hMock_phanta_import('phanta_output/hMock_UHGV_final_merged_outputs/');
hMock_manifest = hMock_manifest(hMock_manifest.total_reads > read_cutoff,:);
hMock_manifest = sortrows(hMock_manifest,'mouse');

%Remove cecum and SI samples
removal_ind = strcmp(hMock_manifest.site,'cecum') | strcmp(hMock_manifest.site,'small_intestine');
hMock_manifest = hMock_manifest(~removal_ind,:);
n_samples = size(hMock_manifest,1);

%Find matching samples for hMock in hCom
for i = 1:size(hMock_manifest,1)
    if i == 1
        matching_ind = strcmp(hCom_manifest.sample,hMock_manifest.sample{i});
    else
        matching_ind  = matching_ind | strcmp(hCom_manifest.sample,hMock_manifest.sample{i});
    end
end
matching_hCom = hCom_manifest(matching_ind,:);

%% Loop through and compute comparative statistics between matching hCom
% and hMock samples

%Define relevant metrics
log_mse_fun = @(v1,v2) mse(log10(v1),log10(v2));
KLdiv = @(P,Q) nansum(P.*log2(P./Q));
JSdiv = @(P,Q) 0.5*KLdiv(P,(Q+P)/2) + 0.5*KLdiv(Q,(Q+P)/2);

levels = {'species_only_read_bacteria','species_only_read_phage'};
hMock_names = hMock_manifest.species_phage{1}.Properties.RowNames;
explained_name = {'percent_explained_bacteria','percent_explained_phage'};
pseudocount = 1e-7;

%Get master labels to set up analyses
[~,bact_master_labels] = ...
    merge_phanta_tables(hCom_manifest{1,levels{1}}{1},...
    hMock_manifest{1,levels{1}}{1});

[~,phage_master_labels] = ...
    merge_phanta_tables(hCom_manifest{1,levels{2}}{1},...
    hMock_manifest{1,levels{2}}{1});

%Set up storage vectors for analysis of individual taxa
always_missing_hCom{1} = true(size(bact_master_labels));
always_missing_hCom{2} = true(size(phage_master_labels));
always_missing_hMock = always_missing_hCom;

times_missing_hCom{1} = zeros(size(bact_master_labels));
times_missing_hCom{2} = zeros(size(phage_master_labels));
times_missing_hMock = times_missing_hCom;

hCom_abun{1} = zeros(size(bact_master_labels));
hCom_abun{2} = zeros(size(phage_master_labels));
hMock_abun{1} = zeros(size(bact_master_labels));
hMock_abun{2} = zeros(size(phage_master_labels));

hCom_pres{1} = zeros(size(bact_master_labels));
hCom_pres{2} = zeros(size(phage_master_labels));
hMock_pres{1} = zeros(size(bact_master_labels));
hMock_pres{2} = zeros(size(phage_master_labels));

%Loop through hMock samples
for i = 1:n_samples

    %Loop through bacterial and viral levels
    for j = 1:length(levels)

        %Get the matching hCom and hMock composition tables
        sample = hMock_manifest.sample{i};
        hMock_tbl = hMock_manifest{i,levels{j}}{1};
        hCom_ind = strcmp(hCom_manifest.sample,sample);
        hCom_tbl = hCom_manifest{hCom_ind,levels{j}}{1};


        %Merge tables
        [out_vec,master_labels{i,j}] = merge_phanta_tables(hCom_tbl,hMock_tbl);
        out_vec_cell{i,j} = out_vec;

        %Find elements only found in hCom or hMock
        hCom_exclusive_ind{i,j} = (out_vec(:,1)' > 0) & (out_vec(:,2)' == 0);
        hMock_exclusive_ind{i,j} = (out_vec(:,1)' == 0) & (out_vec(:,2)' > 0);
        hCom_exclusive{i,j} = master_labels{i,j}(hCom_exclusive_ind{i,j});
        hMock_exclusive{i,j} = master_labels{i,j}(hMock_exclusive_ind{i,j});

        %Get mean of read abundances
        hCom_abun{j} = hCom_abun{j} + out_vec(:,1);
        hMock_abun{j} = hMock_abun{j} + out_vec(:,2);

        %Compute presence
        hCom_pres{j} = hCom_pres{j} + (out_vec(:,1)>0);
        hMock_pres{j} = hMock_pres{j} + (out_vec(:,2)>0);

        %Analyze missing taxa, identifying which are always missing and how
        %many samples a given taxa is missing from
        hCom_missing_ind{i,j} = (out_vec(:,1) == 0) & (out_vec(:,2) > 0);
        hMock_missing_ind{i,j} = (out_vec(:,1) > 0) & (out_vec(:,2) == 0);
        always_missing_hCom{j} = always_missing_hCom{j} & hCom_missing_ind{i,j};
        always_missing_hMock{j} = always_missing_hMock{j} & hMock_missing_ind{i,j};
        times_missing_hCom{j} = times_missing_hCom{j} + hCom_missing_ind{i,j};
        times_missing_hMock{j} = times_missing_hMock{j} + hMock_missing_ind{i,j};

        %Generate pseudocounted version of data
        comp_out_vec = out_vec(sum(out_vec,2) > 0,:);
        comp_out_vec_p = comp_out_vec;
        comp_out_vec_p(comp_out_vec == 0) = pseudocount;
        min_bact_abun = min(comp_out_vec_p(:));
        max_bact_abun = max(comp_out_vec_p(:));

        %Compute various comparative statistics

        %Fraction of hCom abundance nonzero in hMock
        frac_explained = sum(out_vec(out_vec(:,2) > 0,1))/sum(out_vec(:,1));
        hMock_manifest{i,explained_name{j}} = frac_explained;

        %Normal spearman correlation
        spearman_corr(i,j) = corr(comp_out_vec(:,2),comp_out_vec(:,1),'type','Spearman');

        %Jaccard similarity
        shared_ind = (comp_out_vec(:,2)>0) & (comp_out_vec(:,1)>0);
        union_ind = (comp_out_vec(:,2)>0) | (comp_out_vec(:,1)>0);
        jaccard_recon(i,j) = sum(shared_ind)/sum(union_ind);

        %Jensen-Shannon divergence
        rel_vec_1 = out_vec(:,1)/sum(out_vec(:,1));
        rel_vec_2 = out_vec(:,2)/sum(out_vec(:,2));
        JSS(i,j) = 1-JSdiv(rel_vec_1,rel_vec_2);

        %Shared spearman correlation
        if sum(shared_ind) > 1
            shared_spearman_corr(i,j) = corr(comp_out_vec(shared_ind,2),...
                comp_out_vec(shared_ind,1),'type','Spearman');
            log_shared_mse(i,j) = log_mse_fun(comp_out_vec(shared_ind,2),...
                comp_out_vec(shared_ind,1));
            log_pseudocount_mse(i,j) = log_mse_fun(comp_out_vec_p(:,2),...
                comp_out_vec_p(:,1));
        else
            shared_spearman_corr(i,j) = NaN;
            log_shared_mse(i,j) = NaN;
            log_pseudocount_mse(i,j) = NaN;
        end
    end
end

%Assign the comparative statistics to the manifest
hCom_mean{1} = hCom_abun{1}./hCom_pres{1};
hCom_mean{2} = hCom_abun{2}./hCom_pres{2};

hMock_mean{1} = hMock_abun{1}./hMock_pres{1};
hMock_mean{2} = hMock_abun{2}./hMock_pres{2};

hMock_manifest.spearman_bacteria = spearman_corr(:,1);
hMock_manifest.spearman_phage = spearman_corr(:,2);

hMock_manifest.shared_spearman_bacteria = shared_spearman_corr(:,1);
hMock_manifest.shared_spearman_phage = shared_spearman_corr(:,2);

hMock_manifest.jaccard_bacteria = jaccard_recon(:,1);
hMock_manifest.jaccard_phage = jaccard_recon(:,2);

hMock_manifest.JSS_bacteria = (JSS(:,1));
hMock_manifest.JSS_phage = (JSS(:,2));

hMock_manifest.log_shared_mse_bacteria = log_shared_mse(:,1);
hMock_manifest.log_shared_mse_phage = log_shared_mse(:,2);

hMock_manifest.log_pseudocount_mse_bacteria = log_pseudocount_mse(:,1);
hMock_manifest.log_pseudocount_mse_phage = log_pseudocount_mse(:,2);


%% Main figure: plot the hCom-hMock comparisons as abundance scatter plots. 
%this part of the script generates 180 plots, only one (mouse 3, week 1) 
%is shown as Fig. 3AB. These plots are in terms of relative read abundances

pseudocount = 1e-7;
color = [211 95 183]/255;
alpha = 0.5;
load('../common_util/font_config.mat');

%Loop through hMock samples
for i = 1:size(hMock_manifest,1)

    figure
    pos = get(gcf,'Position');
    pos(3) = 1.5*pos(3);
    pos(4) = 0.85*pos(4);
    set(gcf,'Position',pos)

    for j = 1:length(levels)
        sample = hMock_manifest.sample{i};

        out_vec = out_vec_cell{i,j};

        comp_out_vec = out_vec(sum(out_vec,2) > 0,:);

        comp_out_vec_p = comp_out_vec;
        comp_out_vec_p(comp_out_vec == 0) = pseudocount;

        subplot(1,2,j)
        scatter(comp_out_vec_p(:,2),comp_out_vec_p(:,1),'o',...
            'MarkerEdgeColor',color,'MarkerFaceColor',color,...
            'MarkerEdgeAlpha',alpha,'MarkerFaceAlpha',alpha)
        hold on

        bound_vec = [pseudocount,1];
        plot(bound_vec,bound_vec,'k-')

        xlim(bound_vec)
        ylim(bound_vec)
        tick_vec = [1e-7 1e-5 1e-3 1e-1];
        xticks(tick_vec);
        yticks(tick_vec);

        if j == 1
            prefix_str = 'Bacteria';
        else
            prefix_str = 'Phage';
        end
        mult = 0.3;
        text(10^(-6.6),10^(-1),prefix_str,...
            'Interpreter','tex','FontSize',GenFontSize,'FontName',FontName)
        text(10^(-6.6),mult*10^(-1),['JSS = ',...
            num2str(round(JSS(i,j),2))],...
            'Interpreter','tex','FontSize',GenFontSize,'FontName',FontName)
        text(10^(-6.6),mult^2*10^(-1),['Shared rho = ',...
            num2str(round(shared_spearman_corr(i,j),2))],...
            'Interpreter','tex','FontSize',GenFontSize,'FontName',FontName)
        xlabel(['hMock read fraction'],'Interpreter','tex',...
            'FontSize',LabelFontSize,'FontName',FontName)

        if j==1
            ylabel('in vivo hCom2 read fraction',...
                'Interpreter','tex','FontSize',LabelFontSize,'FontName',FontName)
        end
        set(gca,'XScale','log')
        set(gca,'YScale','log')
        set(gca,'FontSize',GenFontSize,'FontName',FontName)
        set(gca,'TickLabelInterpreter','tex')

        if j == 1
            letter = 'A';
        else
            letter = 'B';
        end

        text(-0.2,1.03,letter,'Interpreter','tex','Units','normalized',...
            'FontSize',PanelFontSize,'FontName',FontName)

    end

    name_vector = ['plots/hMock_hCom_scatter/Fig_',sample,'_hMock_hCom_scatter.pdf'];
    exportgraphics(gcf,name_vector,'ContentType','vector')

    close all


end

%% Figure 3CD: look at temporal dynamics of JSS in challenged and non-
%challenged mice.

load('../common_util/font_config.mat');

%Define various metrics to be used
var_prefixes = {'JSS'};
var_labels = {'Jensen-Shannon similarity'};
yaxis_max = 1;
unit_yaxis = 1;
main_fig = 1;
var_suffixes = {'_bacteria','_phage'};
subplot_labels = {'Bacteria','Phage'};

%Split data into challenged and non-challenged
feces_hmock = hMock_manifest(strcmp(hMock_manifest.site,'feces'),:);
ctrl_mice_ind = (feces_hmock.mouse > 5) & (feces_hmock.mouse < 11);
chal_hmock = feces_hmock(~ctrl_mice_ind,:);
ctrl_hmock = feces_hmock(ctrl_mice_ind,:);

%Plotting parameters
colors = [0.5 0.5 0.5; 1 0 0];
fill_between_lines = @(X1,X2,Y1,Y2,C)fill( [transpose(X1) fliplr(transpose(X2))],...
    [transpose(Y1) fliplr(transpose(Y2))], C,'EdgeColor','none','FaceAlpha',0.2);
LineWidth = 1.5;
ytick_vec = [0 0.25 0.5 0.75 1];
xtick_vec = 1:8;

%Make plots for each metric
for i = 1:length(var_prefixes)

    figure
    pos = get(gcf,'Position');
    pos(3) = 1.5*pos(3);
    pos(4) = 0.85*pos(4);
    set(gcf,'Position',pos)
    ylim_vec = [0,yaxis_max(i)];

    if main_fig(i)
        panel_labels = {'C','D'};
    else
        panel_labels = {'A','B'};
    end

    %Loop through phage and bacteria
    for j = 1:length(var_suffixes)
        var_name = [var_prefixes{i},var_suffixes{j}];

        [chal_g,chal_t] = findgroups(chal_hmock.week);
        chal_mean = splitapply(@mean,chal_hmock{:,var_name},chal_g);
        chal_std = splitapply(@std,chal_hmock{:,var_name},chal_g);

        [ctrl_g,ctrl_t] = findgroups(ctrl_hmock.week);
        ctrl_mean = splitapply(@mean,ctrl_hmock{:,var_name},ctrl_g);
        ctrl_std = splitapply(@std,ctrl_hmock{:,var_name},ctrl_g);

        %Plot the metric timeseries
        subplot(1,2,j)
        hold on
        ctrl_high = ctrl_mean + ctrl_std;
        ctrl_low = ctrl_mean - ctrl_std;
        chal_high = chal_mean + chal_std;
        chal_low = chal_mean - chal_std;
        h(1) = fill_between_lines(ctrl_t,ctrl_t, ctrl_high,ctrl_low,colors(1,:));
        h(2) = plot(ctrl_t,ctrl_mean,'-','LineWidth',LineWidth,'Color',colors(1,:));
        plot(ctrl_t,ctrl_mean,'o','LineWidth',LineWidth,...
            'MarkerEdgeColor',colors(1,:),'MarkerFaceColor',colors(1,:),...
            'DisplayName','Unchallenged');
        h(3) = fill_between_lines(chal_t,chal_t, chal_high,chal_low,colors(2,:));
        h(4) = plot(chal_t,chal_mean,'-','LineWidth',LineWidth,'Color',colors(2,:));
        plot(chal_t,chal_mean,'o','LineWidth',LineWidth,...
            'MarkerEdgeColor',colors(2,:),'MarkerFaceColor',colors(2,:),...
            'DisplayName','Challenged');

        h(5) = plot([4.2,4.2],ylim_vec,'r--');
        ylim(ylim_vec)

        if unit_yaxis(i)
            yticks(ytick_vec)
        end

        xticks(xtick_vec)
        for k = 1:length(h)
            set( get( get( h(k), 'Annotation'), 'LegendInformation' ),...
                'IconDisplayStyle', 'off' );
        end

        xlabel('Time (weeks)','Interpreter','tex','FontSize',LabelFontSize,...
            'FontName',FontName)

        if j == 1
            ylabel(var_labels{i},'Interpreter','tex','FontSize',LabelFontSize,...
            'FontName',FontName)
        end

        set(gca,'FontSize',GenFontSize,'FontName',FontName)
        set(gca,'TickLabelInterpreter','tex')

        text(-0.2,1.03,panel_labels{j},'Interpreter','tex','Units',...
            'normalized','FontSize',PanelFontSize,'FontName',FontName)
        text(1-6.6/7,0.6,subplot_labels{j},'Interpreter','tex','Units',...
            'normalized','FontSize',GenFontSize,'FontName',FontName)

        if j == 2
            legend('interpreter','tex','FontSize',GenFontSize,'Location','southeast')
            xlabel('Time (weeks)','Interpreter','tex',...
                'FontSize',LabelFontSize,'FontName',FontName)
        end

    end

    name_vector = ['plots/Fig_reconstruction_',var_prefixes{i},'_timeseries.pdf'];
    exportgraphics(gcf,name_vector,'ContentType','vector')

end



%% Supp figure: temporal dynamics of different metrics in challenged and non-
%challenged mice.

%Define various metrics to be used
var_prefixes = {'spearman','shared_spearman','jaccard',...
    'percent_explained','log_shared_mse','log_pseudocount_mse'};
var_labels = {{'Spearman', 'corr.'},{'Shared', 'Spearman', 'corr.'},{'Jaccard', 'index'},...
    {'Frac. of', 'hCom >0 in', 'hMock'},...
    {'MSE of log10','shared', 'entries'},{'MSE of log10', 'entries (with', 'pseudo-','count)'}};
yaxis_max = [1 1 1 1 0.2,5];
yaxis_type = [1 1 1 1 2 3];
var_suffixes = {'_bacteria','_phage'};
subplot_labels = {'Bacteria','Phage'};

%Split data into challenged and non-challenged
feces_hmock = hMock_manifest(strcmp(hMock_manifest.site,'feces'),:);
ctrl_mice_ind = (feces_hmock.mouse > 5) & (feces_hmock.mouse < 11);
chal_hmock = feces_hmock(~ctrl_mice_ind,:);
ctrl_hmock = feces_hmock(ctrl_mice_ind,:);

%Plotting parameters
colors = [0.5 0.5 0.5; 1 0 0];
fill_between_lines = @(X1,X2,Y1,Y2,C)fill( [transpose(X1) fliplr(transpose(X2))],...
    [transpose(Y1) fliplr(transpose(Y2))], C,'EdgeColor','none','FaceAlpha',0.2);
LineWidth = 1.5;
FontSize = 18;
ylim_vec = [0 1];
ytick_vec = {[0 0.25 0.5 0.75 1], [0 0.05 0.1 0.15 0.2],[0 1 2 3 4 5]};
xtick_vec = 1:8;
panel_labels = {'A','B','C','D','E','F','G','H','I','J','K','L'};

figure
pos = get(gcf,'Position');
pos(3) = 1.8*pos(3);
pos(4) = 3*pos(4);
set(gcf,'Position',pos)

%Make plots for each metric
for i = 1:length(var_prefixes)

    ylim_vec = [0,yaxis_max(i)];

    %Loop through phage and bacteria
    for j = 1:length(var_suffixes)
        var_name = [var_prefixes{i},var_suffixes{j}];

        [chal_g,chal_t] = findgroups(chal_hmock.week);
        chal_mean = splitapply(@mean,chal_hmock{:,var_name},chal_g);
        chal_std = splitapply(@std,chal_hmock{:,var_name},chal_g);

        [ctrl_g,ctrl_t] = findgroups(ctrl_hmock.week);
        ctrl_mean = splitapply(@mean,ctrl_hmock{:,var_name},ctrl_g);
        ctrl_std = splitapply(@std,ctrl_hmock{:,var_name},ctrl_g);

        %Plot the metric timeseries
        subplot(length(yaxis_max),2,2*(i-1)+j)
        hold on
        ctrl_high = ctrl_mean + ctrl_std;
        ctrl_low = ctrl_mean - ctrl_std;
        chal_high = chal_mean + chal_std;
        chal_low = chal_mean - chal_std;
        h(1) = fill_between_lines(ctrl_t,ctrl_t, ctrl_high,ctrl_low,colors(1,:));
        h(2) = plot(ctrl_t,ctrl_mean,'-','LineWidth',LineWidth,'Color',colors(1,:));
        plot(ctrl_t,ctrl_mean,'o','LineWidth',LineWidth,...
            'MarkerEdgeColor',colors(1,:),'MarkerFaceColor',colors(1,:),...
            'DisplayName','Unchallenged');
        h(3) = fill_between_lines(chal_t,chal_t, chal_high,chal_low,colors(2,:));
        h(4) = plot(chal_t,chal_mean,'-','LineWidth',LineWidth,'Color',colors(2,:));
        plot(chal_t,chal_mean,'o','LineWidth',LineWidth,...
            'MarkerEdgeColor',colors(2,:),'MarkerFaceColor',colors(2,:),...
            'DisplayName','Challenged');

        h(5) = plot([4.2,4.2],ylim_vec,'r--');
        ylim(ylim_vec)

        yticks(ytick_vec{yaxis_type(i)})

        for k = 1:length(h)
            set( get( get( h(k), 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
        end

        if j == 1
            ylabel(var_labels{i},'Interpreter','tex',...
                'FontSize',GenFontSize,'FontName',FontName)
        end

        set(gca,'FontSize',GenFontSize,'FontName',FontName)
        set(gca,'TickLabelInterpreter','tex')

        text(-0.15,1.05,panel_labels{2*(i-1)+j},'Interpreter','tex','Units',...
            'normalized','FontSize',PanelFontSize,'FontName',FontName)

        if i == 1
        text(1-6.6/7,0.3,subplot_labels{j},'Interpreter','tex','Units',...
            'normalized','FontSize',GenFontSize,'FontName',FontName)
        end

                    xticks(xtick_vec)
        if i == length(var_prefixes)
            xlabel('Time (weeks)','Interpreter','tex',...
                'FontSize',GenFontSize,'FontName',FontName)
        else
            xticklabels([])
        end

        if i == 1 && j == 1
            legend('interpreter','tex','FontSize',GenFontSize,...
                'Location','southeast','FontName',FontName)
        end

    end


end

name_vector = ['plots/SI_fig_reconstruction_timeseries.pdf'];
exportgraphics(gcf,name_vector,'ContentType','vector')
