%% This script analyzes the Phanta output of Liang et al 2020's stool and 
% VLP data. It generates two datafiles: one that is used in the
% quantification meta-analysis (Fig. 1B) and another that is used in the
% VLP vs. stool absolute abundance analysis (Fig. 1CD). 

clear;clc

% Import the phanta output, the VLP microscopy measurements are included in
% the import manifest
[manifest,virscore] = ...
    wrapped_liang_2020_phanta_import('liang_2020_UHGV_final_merged_outputs/');

%% Generate the Liang et al. 2020 bulk quantification data for the
% meta-analysis plot in figure 1B

%Get only the 1 month and 4 month stool samples
meta_manifest = manifest(manifest.time == 1 | manifest.time == 4,:);
meta_manifest = meta_manifest(strcmp(meta_manifest.sample_type,'stool'),:);

%Generate the standardized study table from this data
study_name = 'liang_stepwise_2020_phanta';
variable_names = {'subject','state','n_sample','study','method','mean_phage_load','sttdev_phage_load'};

subject_names = cellfun(@(x) ['liang_stepwise_2020_subject_',x], ...
    meta_manifest.subject_id,'UniformOutput',false);

states = arrayfun(@(x) ['infant_month_',num2str(x)] ,...
    meta_manifest.time,'UniformOutput',false);

subject_load_std = nan(size(subject_names));
subject_n_sample = ones(size(subject_names));

n_microbe = 0.92e11; %microbe/g feces
subject_mean_load = meta_manifest.phage_to_microbe_ratio*n_microbe;

methods = cell(size(subject_names));
studies = methods;
methods(:) = {'meta_bulk'};
studies(:) = {study_name};

final_table = table(subject_names,states,subject_n_sample,studies,methods,...
    subject_mean_load,subject_load_std,...
    'VariableNames',variable_names);

writetable(final_table,['../../virome_abundance_meta_analysis/final_aggregate_tables/',study_name,'_aggregated.xlsx'],'WriteVariableNames',true);


%% Generate the Liang et al. 2020 data for the absolute abundance 
% scatter plot in Fig 1CD

%Look only at month 4
paired_manifest = manifest(manifest.time == 4,:);
n_microbe = 0.92e11; %microbe/g feces

unique_subjects = unique(paired_manifest.subject_id);

tax_id_vec = (1:size(paired_manifest.species_phage{1},1))';
tax_cell = paired_manifest.species_phage{1}.Properties.RowNames;

%Initiate storage vectors
VLP_rel_col = [];
VLP_rel_phage_col = [];
total_VLP_dens_col = [];
VLP_abs_col = [];
stool_rel_col = [];
stool_rel_phage_col = [];
total_stool_dens_col = [];
stool_abs_col = [];
tax_id_col = [];
virscore_col = [];
subject_col = [];
study_col = [];

paired_manifest.VLP_only_frac = nan(size(paired_manifest.sample));

%Loop through subjects
for i = 1:length(unique_subjects)

    subject_i = unique_subjects{i};

    %Identify the VLP and stool samples for that subject
    VLP_ind = strcmp(paired_manifest.subject_id,subject_i) &...
        strcmp(paired_manifest.sample_type,'VLP');
    stool_ind = strcmp(paired_manifest.subject_id,subject_i) &...
        strcmp(paired_manifest.sample_type,'stool');

    %VLP-based abundances, "_rel" is taxonomic abundance including the rest
    %of the community, "_rel_phage" is taxonomic abundance among only
    %phage, "_abs" is absolute particle count
    VLP_rel = paired_manifest.species_phage{VLP_ind}.Variables;
    VLP_rel_phage = VLP_rel/sum(VLP_rel);
    total_VLP_dens = zeros(size(VLP_rel_phage));
    total_VLP_dens(:) = paired_manifest.viral_density(VLP_ind);
    VLP_abs = total_VLP_dens(1)*VLP_rel_phage;

    %Stool-based abundances, similar to above
    stool_rel = paired_manifest.species_phage{stool_ind}.Variables;
    stool_rel_phage = stool_rel/sum(stool_rel);
    total_stool_dens = zeros(size(stool_rel_phage));
    total_stool_dens(:) = paired_manifest.phage_to_microbe_ratio(stool_ind)*n_microbe;
    stool_abs = total_stool_dens(1)*stool_rel_phage;

    subject_vec = zeros(size(stool_rel));
    subject_vec(:) = str2num(subject_i);

    %Study index of 0 denotes Liang et al. 2020
    study_vec = zeros(size(stool_rel));

    %Store the individual's data in the storage vectors
    VLP_rel_col = [VLP_rel_col; VLP_rel];
    VLP_rel_phage_col = [VLP_rel_phage_col; VLP_rel_phage];
    total_VLP_dens_col = [total_VLP_dens_col; total_VLP_dens];
    VLP_abs_col = [VLP_abs_col; VLP_abs];

    stool_rel_col = [stool_rel_col; stool_rel];
    stool_rel_phage_col = [stool_rel_phage_col; stool_rel_phage];
    total_stool_dens_col = [total_stool_dens_col; total_stool_dens];
    stool_abs_col = [stool_abs_col; stool_abs];

    tax_id_col = [tax_id_col; tax_id_vec];
    virscore_col = [virscore_col; virscore];
    subject_col = [subject_col; subject_vec];
    study_col = [study_col; study_vec];

    %For side analysis, store VLP only fraction
    paired_manifest.VLP_only_frac(VLP_ind) = sum(VLP_rel_phage(stool_rel_phage == 0));

end

%Generate the table for export
VariableNames = {'VLP_rel','VLP_rel_phage','total_VLP_dens','VLP_abs',...
    'stool_rel','stool_rel_phage','total_stool_dens','stool_abs','tax_id',...
    'virscore','subject','study'};
data = table(VLP_rel_col,VLP_rel_phage_col,total_VLP_dens_col,VLP_abs_col,...
    stool_rel_col,stool_rel_phage_col,total_stool_dens_col, stool_abs_col,...
    tax_id_col, virscore_col,subject_col, study_col,...
    'VariableNames',VariableNames);

save('liang_2020_absolute_data.mat','data','tax_cell');

%% Side analysis: look at VMRs between stool and VLPs to assess possible 
% role of VLP contamination as confounder

%Get ICTV taxonomies of UHGV samples
ictv_taxonomies = get_ictv_taxonomies(paired_manifest.species_phage{1}.Properties.RowNames);

%Find microviridae 
microvir_str = 'Microviridae';
microvir_ind = contains(ictv_taxonomies,microvir_str);

%Loop through samples and compute VMR even after removing microviridae
for i = 1:size(paired_manifest,1)

    %Get relative abundance of microviridae (normalized to the total phage
    %community)
    paired_manifest.microvir_abundance(i) = ...
        sum(paired_manifest.species_phage{i}(microvir_ind,1).Variables)...
        ./sum(paired_manifest.species_phage{i}.Variables);

    %Compute VMR without microviridae
    paired_manifest.non_microvir_VMR(i) = paired_manifest.phage_to_microbe_ratio(i)...
        *(1-paired_manifest.microvir_abundance(i));

    %If VLP sample, compute VMR excluding VLP only taxa
    if ~isnan(paired_manifest.VLP_only_frac(i))
        paired_manifest.non_VLP_only_VMR(i) = paired_manifest.phage_to_microbe_ratio(i)...
            *(1-paired_manifest.VLP_only_frac(i));
    end
end

%Identify stool and VLP samples
stool_ind = strcmp(paired_manifest.sample_type,'stool');
VLP_ind = strcmp(paired_manifest.sample_type,'VLP');

%Mean microviridae abundance in VLP samples
mean_VLP_microvir_abun = mean(paired_manifest.microvir_abundance(VLP_ind));

%Compute median VMRs
median_stool_VMR = median(paired_manifest.phage_to_microbe_ratio(stool_ind));
median_VLP_VMR = median(paired_manifest.phage_to_microbe_ratio(VLP_ind));

median_stool_VMR_non_microvir = median(paired_manifest.non_microvir_VMR(stool_ind));
median_VLP_VMR_non_microvir = median(paired_manifest.non_microvir_VMR(VLP_ind));

median_VLP_VMR_non_VLP_only = median(paired_manifest.non_VLP_only_VMR(VLP_ind));
