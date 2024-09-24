%% This script compares hCom2 and Yachida et al. abundances at higher
%taxonomic levels than species

clear;clc
close all

%Read hCom2 data
hCom_manifest = wrapped_hCom_phanta_import('phanta_output/hCom_UHGV_final_merged_outputs/');
hCom_manifest = sortrows(hCom_manifest,'mouse');
read_cutoff = 1e5;
hCom_manifest = hCom_manifest(hCom_manifest.total_reads > read_cutoff,:);
in_vivo_ind = strcmp(hCom_manifest.site,'feces') & strcmp(hCom_manifest.challenged,'n');
hCom_manifest = hCom_manifest(in_vivo_ind,:);

%Get the hCom2 data matrices, normalize to within-type totals (i.e. all
%phage abundances sum to 1)
hCom2_mat_cell = {};
taxa_types = {'species_bacteria','species_phage'};
for j = 1:length(taxa_types)
    for i = 1:size(hCom_manifest,1)
        hCom2_mat_cell{j}(:,i) = hCom_manifest{i,taxa_types{j}}{1}.Variables...
            /sum(hCom_manifest{i,taxa_types{j}}{1}.Variables);
    end
end
hCom2_taxa = {hCom_manifest.species_bacteria{1}.Properties.RowNames,...
    hCom_manifest.species_phage{1}.Properties.RowNames};


%% Read Yachida et al. data

cd('../cohort_analyses/process_yachida_2019/')
yachida_manifest= ...
    wrapped_yachida_2019_phanta_import('yachida_2019_UHGV_final_merged_outputs/');
cd('../../hCom_hMock_analysis/')
yachida_manifest = yachida_manifest(yachida_manifest.total_reads > read_cutoff,:);

%Get the Yachida data matrices, normalize to within-type totals (i.e. all
%phage abundances sum to 1)
yachida_mat_cell = {};
taxa_types = {'species_bacteria','species_phage'};
for j = 1:length(taxa_types)
    for i = 1:size(yachida_manifest,1)
        yachida_mat_cell{j}(:,i) = yachida_manifest{i,taxa_types{j}}{1}.Variables...
            /sum(yachida_manifest{i,taxa_types{j}}{1}.Variables);
    end
end
yachida_taxa = {yachida_manifest.species_bacteria{1}.Properties.RowNames,...
    yachida_manifest.species_phage{1}.Properties.RowNames};


%% Look at genus- and family-level abundance/prevalence in the two datasets

taxa_levels = {'family','genus'};

%Loop through taxonomic levels
for i = 1:length(taxa_levels)

    taxa_str = [taxa_levels{i},'_'];

    %Loop through the taxa types (i.e. bacteria and phage)
    for j = 1:length(taxa_types)

        %Get the set of unique taxa
        yachida_species = yachida_taxa{j};
        hCom2_species = hCom2_taxa{j};
        all_taxa = [yachida_species; hCom2_species];
        unique_taxa{i,j} = get_unique_taxa(all_taxa,taxa_str);

        %Get total abundance classified at a given taxonomic level
        total_yachida_ind{i,j} = contains(yachida_species,taxa_str);
        total_hCom2_ind{i,j} = contains(hCom2_species,taxa_str);
        yachida_total(i,j) = mean(sum(yachida_mat_cell{j}(total_yachida_ind{i,j},:),1));
        hCom2_total(i,j) = mean(sum(hCom2_mat_cell{j}(total_hCom2_ind{i,j},:),1));

        %Loop through each of the taxa
        for k = 1:length(unique_taxa{i,j})
            taxon_k = unique_taxa{i,j}{k};

            %Find out which species are in taxon k
            yachida_ind{i,j}(:,k) = contains(yachida_species,[taxon_k,'|']);
            hCom2_ind{i,j}(:,k) = contains(hCom2_species,[taxon_k,'|']);

            %Get the vector of abundances of that taxon
            yachida_vec = sum(yachida_mat_cell{j}(yachida_ind{i,j}(:,k),:),1);
            hCom2_vec = sum(hCom2_mat_cell{j}(hCom2_ind{i,j}(:,k),:),1);

            %Compute prevalence and abundance (mean abundance computed only
            %in cases where it is detected)
            yachida_prevalence{i,j}(k) = sum(yachida_vec > 0)/length(yachida_vec);
            yachida_abundance{i,j}(k) = mean(yachida_vec(yachida_vec > 0));
            hCom2_prevalence{i,j}(k) = sum(hCom2_vec > 0)/length(hCom2_vec);
            hCom2_abundance{i,j}(k) = mean(hCom2_vec(hCom2_vec > 0));

        end

        %Assemble the table and sort by prevalence in the human samples
        taxa_table{i,j} = table(unique_taxa{i,j},yachida_prevalence{i,j}',...
            yachida_abundance{i,j}',hCom2_prevalence{i,j}',...
            hCom2_abundance{i,j}','VariableNames',{'taxon','yachida_prevalence',...
            'yachida_mean_total_abundance','hCom2_prevalence','hCom2_mean_total_abundance'});
        taxa_table{i,j} = sortrows(taxa_table{i,j},"yachida_prevalence",'descend');

    end

end

%Compute how many of the top 20 most prevalent phage genera in the Yachida
%cohort are present in at least one hCom sample above 0.1%
top_20_prevalent_phage_genera_in_hCom = ...
    sum(taxa_table{2,2}.hCom2_mean_total_abundance(1:20)>1e-3);
