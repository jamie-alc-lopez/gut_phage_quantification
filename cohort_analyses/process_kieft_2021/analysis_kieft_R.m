%% This script analyzes the individual prophage-host data from Kieft et al.
% to extract mean values of R. This generates the values in Table 2.

clear;clc

% Load in the data from Kieft et al. that has individual 
% prophage-host pairs, Table S3B
kieft_data = readtable('data/msystems.00084-22-st003.xlsx','Sheet','TableS3B');

%Get subsets of data that are human gut related and present
hg_data = kieft_data(contains(kieft_data.dataset,{'infant gut','CRC',...
    'HeQ','IjazUZ'}),:);
hg_data = hg_data(~strcmp(hg_data.active,'not present'),:);

%Perform additional QC
QC_ind = (hg_data.host_median_cov>1) & ...
    (hg_data.prophage_cov_breadth > 0.5) & ...
    (hg_data.prophage_median_cov > 1);
hg_data = hg_data(QC_ind,:);

%Relabel the infant dataset as antibiotic or control using Table S4A of
%Kieft et al.
infant_metadata = readtable('data/msystems.00084-22-st004.xlsx',...
    'Sheet','TableS4A','ReadRowNames',true);
infant_ind = strcmp(hg_data.dataset,{'infant gut'});
infant_groups = ...
    infant_metadata(hg_data.reads_set(infant_ind),'treatmentGroup').Variables;
hg_data.treatment_group(infant_ind) = infant_groups;
hg_data.dataset(strcmp(hg_data.treatment_group,'no treatment')) = {'infant_ctrl'};
hg_data.dataset(strcmp(hg_data.treatment_group,'antibiotic treatment')) = {'infant_abx'};

%Relabel the CRC dataset as cancer or control using sequence archive
%metadata from Feng et al. 2015
CRC_metadata = readtable('data/CRC_metadata/filereport_read_run_PRJEB7774_tsv.txt',...
    FileType='text',ReadRowNames=true);
CRC_ind = strcmp(hg_data.dataset,'CRC');
hg_data.treatment_group(CRC_ind) = CRC_metadata.sample_title(hg_data.reads_set(CRC_ind));
healthy_CRC_ind = CRC_ind & strcmp(hg_data.treatment_group,'Stool sample from controls');
cancer_CRC_ind = CRC_ind & ~strcmp(hg_data.treatment_group,'Stool sample from controls');
hg_data.dataset(healthy_CRC_ind) = {'CRC_healthy'};
hg_data.dataset(cancer_CRC_ind) = {'CRC_cancer'};

%Relabel the HeQ dataset as CD or control using supplementary table 1 and 
%sequence archive metadata from He et al. 2017
HeQ_metadata = readtable('data/HeQ_metadata/processed_HeQ_metadata.xlsx',...
    'ReadRowNames',true);
HeQ_ind = strcmp(hg_data.dataset,'HeQ');
hg_data.treatment_group(HeQ_ind) = HeQ_metadata.disease_status(hg_data.reads_set(HeQ_ind));
hg_data.age(HeQ_ind) = HeQ_metadata.age(hg_data.reads_set(HeQ_ind));
%Separate out and exclude patients under 18 and those without recorded ages
HeQ_adult_ind = hg_data.age >= 18; 
healthy_HeQ_ind = HeQ_ind & strcmp(hg_data.treatment_group,'Disease status:Control');
CD_HeQ_ind = HeQ_ind & strcmp(hg_data.treatment_group,'Disease status:CD');
hg_data.dataset(healthy_HeQ_ind & HeQ_adult_ind) = {'HeQ_healthy'};
hg_data.dataset(CD_HeQ_ind & HeQ_adult_ind) = {'HeQ_CD'};
hg_data.dataset(HeQ_ind & ~HeQ_adult_ind) = {'HeQ_other'};

%Relabel the Ijaz dataset as children with CD
IjazUZ_metadata = readtable('data/IjazUZ_metadata/processed_IjazUZ_metadata.xlsx',...
    'ReadRowNames',true);
IjazUZ_ind = strcmp(hg_data.dataset,'IjazUZ');
hg_data.treatment_group(IjazUZ_ind) = ...
    IjazUZ_metadata.disease_status(hg_data.reads_set(IjazUZ_ind));
hg_data.age(IjazUZ_ind) = IjazUZ_metadata.age(hg_data.reads_set(IjazUZ_ind));
minor_ind = hg_data.age < 18;
IjazUZ_child_CD_ind = IjazUZ_ind & minor_ind & ...
    contains(hg_data.treatment_group,'Crohn');
other_IjazUZ_ind = IjazUZ_ind & ~IjazUZ_child_CD_ind;
hg_data.dataset(IjazUZ_child_CD_ind) = {'IjazUZ_child_CD'};
hg_data.dataset(other_IjazUZ_ind) = {'IjazUZ_other'};

%% Compute R statistics across different dataset groupings, save table
dataset_groups = {{'CRC_healthy','HeQ_healthy'},'CRC_healthy',...
    'CRC_cancer','HeQ_healthy','HeQ_CD','IjazUZ_child_CD',...
    {'infant_ctrl','infant_abx'},'infant_ctrl','infant_abx'};
dataset_group_names = {'adult_healthy','CRC_healthy','CRC_cancer',...
    'HeQ_healthy','HeQ_CD','IjazUZ_child_CD','all_infant','infant_ctrl','infant_abx'};

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

non_sample_vars = ...
    ~contains(mean_R_table.Properties.VariableNames,{'n_samples','n_pairs'});

mean_R_table(:,non_sample_vars).Variables = ...
    round(mean_R_table(:,non_sample_vars).Variables,3,'significant');

writetable(mean_R_table,'kieft_induction_table.xlsx',...
    'WriteRowNames',true,'WriteVariableNames',true)

%% Side analysis: compute median number of induced phage per sample 
% across cohorts

kieft_sample_data = readtable('data/msystems.00084-22-st003.xlsx','Sheet','TableS3A');
hg_sample_data = kieft_sample_data(contains(kieft_sample_data.dataset,{'infant gut','CRC','HeQ','IjazUZ'}),:);

[G,ID] = findgroups(hg_sample_data.dataset);
median_active = splitapply(@median,hg_sample_data.active,G);
