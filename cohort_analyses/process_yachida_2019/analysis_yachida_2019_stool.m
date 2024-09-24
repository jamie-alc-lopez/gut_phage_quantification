%% This script analyzes the Phanta output of Yachida et al; 2019's stool
% data. It generates a data file that is used in the
% quantification meta-analysis (Fig. 1B)

clear;clc

[manifest,virscore] = ...
    wrapped_yachida_2019_phanta_import('yachida_2019_UHGV_final_merged_outputs/');

%% Generate the Yachida et al. 2019 bulk quantification for the meta-analysis

study_name = 'yachida_2019_phanta';
variable_names = {'subject','state','n_sample','study','method','mean_phage_load','sttdev_phage_load'};

%Generate the standardized study table from this data
subject_names = cellfun(@(x) ['yachida_2019_subject_',x], ...
    manifest.sample,'UniformOutput',false);

subject_load_std = nan(size(subject_names));
subject_n_sample = ones(size(subject_names));

n_microbe = 0.92e11; %microbe/g feces
subject_mean_load = manifest.phage_to_microbe_ratio*n_microbe;

methods = cell(size(subject_names));
studies = methods;
states = methods;
methods(:) = {'meta_bulk'};
studies(:) = {study_name};
states(:) = {'adult'};

final_table = table(subject_names,states,subject_n_sample,studies,methods,...
    subject_mean_load,subject_load_std,...
    'VariableNames',variable_names);

writetable(final_table,['../../virome_abundance_meta_analysis/final_aggregate_tables/',study_name,'_aggregated.xlsx'],'WriteVariableNames',true);
