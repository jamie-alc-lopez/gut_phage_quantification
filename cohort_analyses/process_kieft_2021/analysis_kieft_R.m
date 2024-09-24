%% This script analyzes the individual prophage-host data from Kieft et al.
% to extract mean values of R. This generates the values in Table 2.

clear;clc

% Load in the data from Kieft et al. that has individual 
% prophage-host pairs, Table S3B
kieft_data = readtable('data/msystems.00084-22-st003.xlsx','Sheet','TableS3B');

%Load in infant sample metadata, Table S4A
infant_metadata = readtable('data/msystems.00084-22-st004.xlsx',...
    'Sheet','TableS4A','ReadRowNames',true);

%Get subsets of data that are human gut related and present
hg_data = kieft_data(contains(kieft_data.dataset,{'infant gut','CRC',...
    'HeQ','IjazUZ'}),:);
hg_data = hg_data(~strcmp(hg_data.active,'not present'),:);

%Perform additional QC
QC_ind = (hg_data.host_median_cov>1) & ...
    (hg_data.prophage_cov_breadth > 0.5) & ...
    (hg_data.prophage_median_cov > 1);
hg_data = hg_data(QC_ind,:);

%Relabel the infant dataset as antibiotic or control
infant_ind = strcmp(hg_data.dataset,{'infant gut'});
infant_groups = ...
    infant_metadata(hg_data.reads_set(infant_ind),'treatmentGroup').Variables;
hg_data.treatment_group(infant_ind) = infant_groups;
hg_data.dataset(strcmp(hg_data.treatment_group,'no treatment')) = {'infant_ctrl'};
hg_data.dataset(strcmp(hg_data.treatment_group,'antibiotic treatment')) = {'infant_abx'};

%Compute R statistics across different dataset grouping
dataset_groups = {{'CRC','HeQ','IjazUZ','infant_ctrl','infant_abx'},...
    {'CRC','HeQ','IjazUZ'},'CRC','HeQ','IjazUZ',...
    {'infant_ctrl','infant_abx'},'infant_ctrl','infant_abx'};
dataset_group_names = {'all','adult','CRC','HeQ','IjazUZ','all_infant',...
    'infant_ctrl','infant_abx'};

%Loop through the dataset groups
for i = 1:length(dataset_group_names)
    subdata = hg_data(contains(hg_data.dataset,dataset_groups{i}),:);

    n_samples(i) = length(unique(subdata.reads_set));
    n_pairs(i) = size(subdata,1);
    base_mean_R(i) = mean(subdata.prophage_host_ratio);
    median_R(i) = median(subdata.prophage_host_ratio);
    cond_mean_R(i) = mean(subdata.prophage_host_ratio(subdata.prophage_host_ratio >= 1));
    clipped_R = subdata.prophage_host_ratio;
    clipped_R(clipped_R < 1) = 1;
    clipped_mean_R(i) = mean(clipped_R);

end

%Make Table 2
mean_R_table = table(n_samples',n_pairs',base_mean_R',median_R',...
    cond_mean_R',clipped_mean_R','RowNames',dataset_group_names,...
    'VariableNames',{'n_samples','n_pairs','base_mean_R','median_R',...
    'cond_mean_R','clipped_mean_R'});

%% Side analysis: compute median number of induced phage per sample 
% across cohorts

kieft_sample_data = readtable('data/msystems.00084-22-st003.xlsx','Sheet','TableS3A');
hg_sample_data = kieft_sample_data(contains(kieft_sample_data.dataset,{'infant gut','CRC','HeQ','IjazUZ'}),:);

[G,ID] = findgroups(hg_sample_data.dataset);
median_active = splitapply(@median,hg_sample_data.active,G);
