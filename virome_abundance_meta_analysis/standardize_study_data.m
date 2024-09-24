
%This script standardizes study data for abundance measurements originating
%from SI files or digitization of figures. For the Phanta quantifications
%of Liang et al. 2020, Shkoporov et al. 2019, and Yachida et al. 2019, see
%the scripts in the relevant "cohort_analysis" folders

%Units for meta_VLP is VLP genomes/g
%units for bulk is viral genomes/g
%Units for EFM is VLP/g

%% Shkoporov et al 2019 meta VLP

clear;clc

study_name = 'shkoporov_2019';
variable_names = {'subject','state','n_sample','study','method',...
    'mean_phage_load','sttdev_phage_load'};

%Import Shkoporov SI table
opts = detectImportOptions(['study_data/',study_name,...
    '/1-s2.0-S1931312819304767-mmc2.xls'],'Sheet','import_form');
opts.VariableTypes{34} = 'double';
data = readtable(['study_data/',study_name,...
    '/1-s2.0-S1931312819304767-mmc2.xls'],opts,'Sheet','import_form');

%Get only samples with VLP quant
load_data = data(~isnan(data.Total_vir_load_genomes_g),:);

%Compute subject statistics
[G,ID] = findgroups(load_data.Subject);
subject_load_averages = splitapply(@mean,load_data.Total_vir_load_genomes_g,G);
subject_load_std = splitapply(@std,load_data.Total_vir_load_genomes_g,G);
subject_n_sample = splitapply(@length,load_data.Total_vir_load_genomes_g,G);
subject_names = strcat([study_name,'_subject_'],arrayfun(@num2str,ID,'UniformOutput',false));

%Make and save table
methods = cell(size(subject_names));
studies = methods;
states = methods;
methods(:) = {'meta_VLP'};
studies(:) = {study_name};
states(:) = {'adult'};

final_table = table(subject_names,states,subject_n_sample,studies,methods,...
    subject_load_averages,subject_load_std,...
    'VariableNames',variable_names);
writetable(final_table,['final_aggregate_tables/',study_name,'_aggregated.xlsx'],'WriteVariableNames',true);

%% Liang et al 2020 (stepwise) EFM VLP
%The data for the 2-5 year old children were not provided on a per-subject
%basis are thus not included in this analysis

clear;clc
study_name = 'liang_stepwise_2020';
variable_names = {'subject','state','n_sample','study','method','mean_phage_load','sttdev_phage_load'};

%Load data and metadata
metadata = readtable(['study_data/',study_name,'/41586_2020_2192_MOESM2_ESM.xlsx'],'Sheet',1,'ReadRowNames',true);
data = readtable(['study_data/',study_name,'/41586_2020_2192_MOESM2_ESM.xlsx'],'Sheet',2,'ReadRowNames',true);
data.subject = metadata{data.Properties.RowNames,'Subject_id'};

%Reformat labels
data.state = data.Study_group;
data.state(strcmp(data.Study_group,'Month 0')) = {'infant_month_0'};
data.state(strcmp(data.Study_group,'Month 1')) = {'infant_month_1'};
data.state(strcmp(data.Study_group,'Month 4')) = {'infant_month_4'};
subject_names = strcat([study_name,'_subject_'],arrayfun(@num2str,data.subject,'UniformOutput',false));


%Make and save table 
subject_load_std = nan(size(subject_names));
subject_n_sample = ones(size(subject_names));
methods = cell(size(subject_names));
studies = methods;
methods(:) = {'EFM'};
studies(:) = {study_name};

final_table = table(subject_names,data.state,subject_n_sample,studies,methods,...
    data.VLP_count_per_gram,subject_load_std,...
    'VariableNames',variable_names);
writetable(final_table,['final_aggregate_tables/',study_name,'_aggregated.xlsx'],'WriteVariableNames',true);


%% Bikel et al. 2021

%This spreadsheet was manually constructed using data from Table S3 (see
%relevant folder). Entries where multiple subject samples were pooled were
%ignored


%% Kim et al. 2011

%This spreadsheet was manually constructed using data digitally extracted from
%Figure 3


%% Hoyles et al 2014

%VLP data digitally extracted from figure 1


