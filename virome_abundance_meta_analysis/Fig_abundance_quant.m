%% Load in and merge study quant tables

clear;clc

table_folder = 'final_aggregate_tables';
aggregate_files = dir(table_folder);

for i = 1:length(aggregate_files)

    if contains(aggregate_files(i).name,'.xlsx')

        tmp_data = readtable([table_folder,'/',aggregate_files(i).name]);

        if ~exist('data','var')
            data = tmp_data;
        else
            data = [data;tmp_data];
        end
    end
end

%Add in the study year
study_dates = cellfun(@(x) strsplit(strrep(x,'_phanta',''),'_'),data.study,'UniformOutput',false);
study_dates = cellfun(@(x) str2double(x{end}),study_dates,'UniformOutput',false);
data.study_date = cell2mat(study_dates); 

%Standard microbial density to be used, from Sender et al.
rho = 0.92e11; %cells/g
rho_SEM = 0.19; %Relative standard error of the mean for microbes


%% Compute summary statistics across different methods, studies and
% subject states

%To get per-study averages, we take the average of all per-subject values.
%If one subject has multiple samples, we just compute the average. 

%Summary of state-method-study combinations
[G1, study_summary] = findgroups(data(:,{'state','study','method'}));
study_summary.mean = splitapply(@mean,data.mean_phage_load,G1);
study_summary.std = splitapply(@std,data.mean_phage_load,G1);
study_summary.pop_CV = study_summary.std./study_summary.mean;
study_summary.n_subjects = splitapply(@length,data.mean_phage_load,G1);
study_summary.SEM = study_summary.std./sqrt(study_summary.n_subjects);
study_summary.rel_SEM = study_summary.SEM./study_summary.mean;
study_summary.temporal_CV = splitapply(@mean,data.sttdev_phage_load...
    ./data.mean_phage_load,G1);

%Compute post-infancy summary stats for EFM
post_infancy_EFM_ind = strcmp(study_summary.method,'EFM') & ...
    (strcmp(study_summary.state,'adult') | strcmp(study_summary.state,'child_7_9_yr'));
post_infancy_EFM_mean = mean(study_summary.mean(post_infancy_EFM_ind));
post_infancy_EFM_SEM = std(study_summary.mean(post_infancy_EFM_ind))...
    /sqrt(sum(post_infancy_EFM_ind));
post_infancy_EFM_SEM_pct = post_infancy_EFM_SEM/post_infancy_EFM_mean*100;
post_infancy_EFM_CV = mean(study_summary.pop_CV(post_infancy_EFM_ind));

%Summary of state-method combinations across studies
[G2,measure_summary] = findgroups([study_summary(:,{'state','method'})]);
measure_summary.mean = splitapply(@mean,study_summary.mean,G2);
measure_summary.n_studies = splitapply(@length,study_summary.mean,G2);
measure_summary.SEM = splitapply(@std,study_summary.mean,G2)...
    ./measure_summary.n_studies;
measure_summary.rel_SEM = measure_summary.SEM./measure_summary.mean;
measure_summary.mean_pop_CV = splitapply(@mean,study_summary.pop_CV,G2);

%Get mean V/P
adult_ind = strcmp(measure_summary.state,'adult');
EFM_ind = strcmp(measure_summary.method,'EFM');
meta_VLP_ind = strcmp(measure_summary.method,'meta_VLP');
meta_bulk_ind = strcmp(measure_summary.method,'meta_bulk');
VLP_EFM_over_bulk = measure_summary.mean(adult_ind&EFM_ind)...
    ./measure_summary.mean(adult_ind&meta_bulk_ind);
VLP_seq_over_bulk = measure_summary.mean(adult_ind&meta_VLP_ind)...
    ./measure_summary.mean(adult_ind&meta_bulk_ind);


%% Plot Fig 1B

state_order = {'infant_month_0','infant_month_1','infant_month_4','child_7_9_yr','adult'};
state_labels = {'Month 0','Month 1','Month 4','Year 7-9', 'Adult'};
method_order = {'EFM','meta_VLP','meta_bulk'};
method_colors = [165 149 71; 198 74 94; 69 128 173]/255;

[G1,state_ID,study_ID,method_ID,date_ID] = ...
    findgroups(data.state,data.study,data.method,data.study_date);
g = 1:length(study_ID);

figure
hold on
x_pos = 1;
load('../common_util/font_config.mat');
study_space = 1.2;
extra_state_space = 3;
tick_record =[];

%Loop through different subject states
for i = 1:length(state_order)

    state_ind = strcmp(state_ID,state_order{i});
    matching_studies = study_ID(state_ind);
    matching_methods = method_ID(state_ind);
    g_i = g(state_ind);

    %Loop all measurement methods for this state
    x_pos_record = [];
    for j = 1:length(method_order)
        method = method_order{j};
        g_ij = g_i(find(strcmp(matching_methods,method)));

        %Proceed if this state has the matching method
        if ~isempty(g_ij)
            loop_years = date_ID(g_ij);
            [~,chrono_ind] = sort(loop_years,'ascend');
            chrono_g_ij = g_ij(chrono_ind);

            %Loop through studies using this method/state pair, ordered by
            %year
            for k = 1:length(chrono_g_ij)

                g_ijk = chrono_g_ij(k);
                phage_load_k = data.mean_phage_load(G1 == g_ijk);

                %For the Liang infant study, set zero phage count to the
                %limit of detection
                if strcmp(study_ID(g == g_ijk),'liang_stepwise_2020')
                    phage_load_k(phage_load_k == 0) = 6.6e6;
                end

                x_pos_record = [x_pos_record; x_pos];
                x_pos = x_pos + study_space; 

                violinplot_modified(phage_load_k,x_pos,ones(size(phage_load_k)),...
                    'ViolinColor',method_colors(j,:),'ViolinAlpha',0.2,'Width',0.5);
            end

        end

    end

    tick_record(i) = mean(x_pos_record);
    x_pos = x_pos + extra_state_space;

end


DisplayNames = {'Particle count',...
    'Particle genomes','Stool phage genomes',...
    'Ref. microbe density'};

h(1) = plot(800,800,'o','Color',method_colors(1,:),'MarkerFaceColor',...
    method_colors(1,:),'DisplayName',DisplayNames{1});
h(2) = plot(800,800,'o','Color',method_colors(2,:),'MarkerFaceColor',...
    method_colors(2,:),'DisplayName',DisplayNames{2});
h(3) = plot(800,800,'o','Color',method_colors(3,:),'MarkerFaceColor',...
    method_colors(3,:),'DisplayName',DisplayNames{3});

h(4) = plot([4.4,25],[0.92e11 0.92e11],'--','Color',[0.5 0.5 0.5],...
    'LineWidth',1.5,'DisplayName',DisplayNames{4});

leg = legend(h,'Interpreter','none','FontSize',GenFontSize,...
    'FontName',FontName,'Position',[0.5134    0.1357    0.3866    0.2083]);
box off
xticks(tick_record)
xticklabels(state_labels)
set(gca,'TickLabelInterpreter','none')
set(gca,'YScale','log')

xlim([0,x_pos - 3])
ylim([1e6,1e13])
yticks(10.^(6:13))
xlabel('Age group','Interpreter','none','FontSize',LabelFontSize,'FontName',FontName)
ylabel({'Estimated phage per gram feces','(particle or genome)'},...
    'Interpreter','none','FontSize',LabelFontSize,'FontName',FontName)

set(gca,'FontSize',GenFontSize,'FontName',FontName)

text(-0.17,1.03,'B','Interpreter','none','Units','normalized',...
    'FontSize',PanelFontSize,'FontName',FontName,'FontWeight','normal')

name_vector = ['plots/Fig_virome_quant_violin.pdf'];
exportgraphics(gcf,name_vector,'ContentType','vector')

