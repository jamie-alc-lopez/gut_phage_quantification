%% This script analyzes the Phanta output of Liang et al 2019's stool and 
% VLP data. It generates two datafiles: one that is used in the
% quantification meta-analysis (Fig. 1B) and another that is used in the
% VLP vs. stool absolute abundance analysis (Fig. 1CD). 

clear;clc

%Load in the phanta quantification
[manifest,virscore] = ...
    wrapped_shkoporov_2019_phanta_import('shkoporov_2019_UHGV_final_merged_outputs/');

%Only look at month 8
manifest = manifest(strcmp(manifest.experimental_factor,'Month_8'),:);
subjects = unique(manifest.host_subject_id);

%Load in the VLP abundance data 
%Note that they did not measure VLP genomes at timepoint 8, so instead we
%using the average from the remaining timepoints. 
abun_file = '../../virome_abundance_meta_analysis/final_aggregate_tables/shkoporov_2019_aggregated.xlsx';
abun_data = readtable(abun_file);

%% Generate the Shkoporov et al. 2019 bulk quantification data for the
% meta-analysis plot in figure 1B

study_name = 'shkoporov_2019_phanta';
variable_names = {'subject','state','n_sample','study','method','mean_phage_load','sttdev_phage_load'};

data = manifest(strcmp(manifest.type,'bulk'),:);

subject_names = cellfun(@(x) ['shkoporov_2019_subject_',x(end-2:end)], data.host_subject_id,'UniformOutput',false);

subject_load_std = nan(size(subject_names));
subject_n_sample = ones(size(subject_names));

n_microbe = 0.92e11; %microbe/g feces
subject_mean_load = data.phage_to_microbe_ratio*n_microbe;

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


%% Generate the Shkoporov et al. 2019 data for the absolute abundance
% scatter plot in Fig. 1CD

unique_subjects = unique(manifest.host_subject_id);
n_microbe = 0.92e11; %microbe/g feces
tax_id_vec = (1:size(manifest.species_phage{1},1))';
tax_cell = manifest.species_phage{1}.Properties.RowNames;

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

manifest.VLP_only_frac = nan(size(manifest.sample));

%Loop through subjects
for i = 1:length(unique_subjects)

    subject_i = unique_subjects{i};

    %Identify the VLP and stool samples for that subject
    phanta_ind = strcmp(manifest.host_subject_id,subject_i);
    VLP_ind = strcmp(manifest.host_subject_id,subject_i) & strcmp(manifest.type,'MDA_VLP');
    stool_ind = strcmp(manifest.host_subject_id,subject_i) & strcmp(manifest.type,'bulk');

    %Get the corresponding viral load from the subject
    subject_num = subject_i(end-2:end);
    abun_ind  = contains(abun_data.subject,subject_num);
    V = abun_data{abun_ind,"mean_phage_load"};

    %VLP-based abundances, "_rel" is taxonomic abundance including the rest
    %of the community, "_rel_phage" is taxonomic abundance among only
    %phage, "_abs" is absolute genome count 
    VLP_rel = manifest.species_phage{VLP_ind}.Variables;
    VLP_rel_phage = VLP_rel/sum(VLP_rel);
    total_VLP_dens = zeros(size(VLP_rel_phage));
    total_VLP_dens(:) = V;
    VLP_abs = total_VLP_dens(1)*VLP_rel_phage;

    %Stool-based abundances, similar to above
    stool_rel = manifest.species_phage{stool_ind}.Variables;
    stool_rel_phage = stool_rel/sum(stool_rel);
    total_stool_dens = zeros(size(stool_rel_phage));
    total_stool_dens(:) = manifest.phage_to_microbe_ratio(stool_ind)*n_microbe;
    stool_abs = total_stool_dens(1)*stool_rel_phage;

    subject_vec = zeros(size(stool_rel));
    subject_vec(:) = str2num(subject_num);

    %Study index of 1 denotes Shkoporov et al. 2019
    study_vec = ones(size(stool_rel));

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
    manifest.VLP_only_frac(VLP_ind) = sum(VLP_rel_phage(stool_rel_phage == 0));

end

%Generate the table for export
VariableNames = {'VLP_rel','VLP_rel_phage','total_VLP_dens','VLP_abs',...
    'stool_rel','stool_rel_phage','total_stool_dens','stool_abs','tax_id',...
    'virscore','subject','study'};
data = table(VLP_rel_col,VLP_rel_phage_col,total_VLP_dens_col,VLP_abs_col,...
    stool_rel_col,stool_rel_phage_col,total_stool_dens_col, stool_abs_col,...
    tax_id_col, virscore_col,subject_col, study_col,...
    'VariableNames',VariableNames);

save('shkoporov_2019_absolute_data.mat','data','tax_cell');


%% Side analysis: look at VMRs between stool and VLPs to assess possible 
% role of VLP contamination as confounder

%Get ICTV taxonomies of UHGV samples
ictv_taxonomies = get_ictv_taxonomies(manifest.species_phage{1}.Properties.RowNames);

%Find microviridae 
microvir_str = 'Microviridae';
microvir_ind = contains(ictv_taxonomies,microvir_str);

%Loop through samples and compute VMR even after removing microviridae
for i = 1:size(manifest,1)

    %Get relative abundance of microviridae (normalized to the total phage
    %community)
    manifest.microvir_abundance(i) = ...
        sum(manifest.species_phage{i}(microvir_ind,1).Variables)...
        ./sum(manifest.species_phage{i}.Variables);

    %Compute VMR without microviridae
    manifest.non_microvir_VMR(i) = manifest.phage_to_microbe_ratio(i)...
        *(1-manifest.microvir_abundance(i));

    %If VLP sample, compute VMR excluding VLP only taxa
    if ~isnan(manifest.VLP_only_frac(i))
        manifest.non_VLP_only_VMR(i) = manifest.phage_to_microbe_ratio(i)...
            *(1-manifest.VLP_only_frac(i));
    end
end

%Identify stool and VLP samples
stool_ind = strcmp(manifest.type,'bulk');
VLP_ind = strcmp(manifest.type,'MDA_VLP');

%Mean microviridae abundance in VLP samples
mean_VLP_microvir_abun = mean(manifest.microvir_abundance(VLP_ind));

%Compute median VMRs
median_stool_VMR = median(manifest.phage_to_microbe_ratio(stool_ind));
median_VLP_VMR = median(manifest.phage_to_microbe_ratio(VLP_ind));

median_stool_VMR_non_microvir = median(manifest.non_microvir_VMR(stool_ind));
median_VLP_VMR_non_microvir = median(manifest.non_microvir_VMR(VLP_ind));

median_VLP_VMR_non_VLP_only = median(manifest.non_VLP_only_VMR(VLP_ind));
