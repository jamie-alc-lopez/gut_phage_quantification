%% This script analyzes the Phanta output of Garmaeva et al 2024's stool and
% VLP data. It generates a datafile used in the non-MDA scatterplot

clear;clc

% Import the phanta output
[manifest,virscore] = ...
    wrapped_garmaeva_2024_phanta_import('garmaeva_2024_UHGV_final_merged_outputs/');


%% Generate the Garmaeva et al. 2024 data for the relative abundance plot
% scatter plot in Fig. S3. This processing is basically identical to that
% done with Fig. 1C, but we impute n_phage as the healthy human mean 

n_microbe = 0.92e11; %microbe/g feces
n_phage = 2e9; %phage particles/g feces

%Minimum VMR for enrichment to be considered successful
VLP_VMR_min = 30;

[group_index,group_subjects,group_timepoints] = findgroups(manifest.subject_id,manifest.timepoint);
[group_counts,sample_groups] = groupcounts(group_index);

total_VLP_ind = strcmp(manifest.sample_type,'vlp');
total_bulk_ind = strcmp(manifest.sample_type,'bulk');

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
phenotype_col = [];

manifest.VLP_only_frac = nan(size(manifest.sample));

%Loop through subjects
for i = 1:length(sample_groups)

    %Check whether there is both a VLP and paired sample
    if group_counts(i) > 1

        %Get group info
        group_i = sample_groups(i);
        group_i_ind = group_index == group_i;
        subject_i = group_subjects{group_i};
        timepoint_i = group_timepoints{group_i};

        %Identify the VLP and stool samples for that subject-timepoint
        %combo
        VLP_ind = group_i_ind & total_VLP_ind;
        stool_ind = group_i_ind & total_bulk_ind;

        VLP_VMR_i = manifest.phage_to_microbe_ratio(VLP_ind);
        sample_type_i = manifest.sample_type(VLP_ind);
        phenotype_i = manifest.phenotype(VLP_ind);

        if VLP_VMR_i > VLP_VMR_min

            %VLP-based abundances, "_rel" is taxonomic abundance including the rest
            %of the community, "_rel_phage" is taxonomic abundance among only
            %phage, "_abs" is absolute particle count
            VLP_rel = manifest.species_phage{VLP_ind}.Variables;
            VLP_rel_phage = VLP_rel/sum(VLP_rel);
            total_VLP_dens = zeros(size(VLP_rel_phage));
            total_VLP_dens(:) = n_phage;
            VLP_abs = total_VLP_dens(1)*VLP_rel_phage;

            %Stool-based abundances, similar to above
            stool_rel = manifest.species_phage{stool_ind}.Variables;
            stool_rel_phage = stool_rel/sum(stool_rel);
            total_stool_dens = zeros(size(stool_rel_phage));
            total_stool_dens(:) = manifest.phage_to_microbe_ratio(stool_ind)*n_microbe;
            stool_abs = total_stool_dens(1)*stool_rel_phage;

            subject_time_cell = cell(size(stool_rel));
            subject_time_cell(:) = {[subject_i,'-',timepoint_i]};

            phenotype_cell = cell(size(stool_rel));
            phenotype_cell(:) = phenotype_i;

            %Study index of 2 denotes Garmaeva et al. 2024
            study_vec = zeros(size(stool_rel)) + 2;

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
            subject_col = [subject_col; subject_time_cell];
            study_col = [study_col; study_vec];

            phenotype_col = [phenotype_col; phenotype_cell];

            %For side analysis, store VLP only fraction
            manifest.VLP_only_frac(VLP_ind) = sum(VLP_rel_phage(stool_rel_phage == 0));
        end

    end

end

%Generate the table for export
VariableNames = {'VLP_rel','VLP_rel_phage','total_VLP_dens','VLP_abs',...
    'stool_rel','stool_rel_phage','total_stool_dens','stool_abs','tax_id',...
    'virscore','subject','study'};
all_data = table(VLP_rel_col,VLP_rel_phage_col,total_VLP_dens_col,VLP_abs_col,...
    stool_rel_col,stool_rel_phage_col,total_stool_dens_col, stool_abs_col,...
    tax_id_col, virscore_col,subject_col, study_col,...
    'VariableNames',VariableNames);

%Subset and save data
data_types = {'adult','infant'};
data_phenotypes = {'Mother','Infant'};
for i = 1:length(data_types)

    data = all_data(strcmp(phenotype_col,data_phenotypes{i}),:);

    save(['garmaeva_2024_',data_types{i},'_absolute_data.mat'],'data','tax_cell');

end

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
stool_ind = strcmp(manifest.sample_type,'bulk');
VLP_ind = strcmp(manifest.sample_type,'vlp');

%Mean microviridae abundance in VLP samples
mean_VLP_microvir_abun = mean(manifest.microvir_abundance(VLP_ind));

%Compute median VMRs
median_stool_VMR = median(manifest.phage_to_microbe_ratio(stool_ind));
median_VLP_VMR = median(manifest.phage_to_microbe_ratio(VLP_ind));

median_stool_VMR_non_microvir = median(manifest.non_microvir_VMR(stool_ind));
median_VLP_VMR_non_microvir = median(manifest.non_microvir_VMR(VLP_ind));

median_VLP_VMR_non_VLP_only = median(manifest.non_VLP_only_VMR(VLP_ind));
