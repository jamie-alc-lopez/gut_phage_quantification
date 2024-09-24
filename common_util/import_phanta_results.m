function [mod_manifest,virscore] = import_phanta_results(tax_file_cell,...
    manifest_file,virscore_file,total_reads_file)

%This function takes a manifest and phanta file locations and returns a
%manifest with the community compositions and stats as table entries

%Turn off warning for valid name modification, the script deals with this
%internally by using matlab.lang.makeValidName
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames')

%Load in the total reads table
total_reads_table = readtable(total_reads_file,'ReadRowNames',true,...
    'ReadVariableNames',true,'FileType','text','VariableNamingRule','modify');

%Define the superkingdom strings, this covers default MGV and
%uhggv2_uhgv_mq
sk_vir = {'superkingdom_Viruses'};
sk_bac = {'superkingdom_Bacteria','superkingdom_d__Bacteria'};
sk_euk = {'superkingdom_Eukaryota'};
sk_arc = {'superkingdom_Archaea','superkingdom_d__Archaea'};
sp_hum = {'species_Homo sapiens'};

%Define root strings needed to determine unclassified read frac
root_str = {'root_root|root_root','no rank_root|no rank_root'};
cell_root_str = {'root_root|no rank_cellular organisms','no rank_root|no rank_cellular organisms'};
non_cell_root_str = {'root_root|superkingdom_Viruses', 'no rank_root|superkingdom_Viruses'};

%Loop through taxa levels
for i = 1:length(tax_file_cell)
    tax_table{i} = readtable(tax_file_cell{i},'ReadRowNames',true,...
        'ReadVariableNames',true,'Delimiter','\t','VariableNamingRule','modify');

    phage_ind = contains(tax_table{i}.Properties.RowNames,sk_vir);
    bacteria_ind = contains(tax_table{i}.Properties.RowNames,sk_bac);
    eukarya_ind = contains(tax_table{i}.Properties.RowNames,sk_euk);
    archaea_ind = contains(tax_table{i}.Properties.RowNames,sk_arc);
    phage_table{i} = tax_table{i}(phage_ind,:);
    bacteria_table{i} = tax_table{i}(bacteria_ind,:);
    eukarya_table{i} = tax_table{i}(eukarya_ind,:);
    archaea_table{i} = tax_table{i}(archaea_ind ,:);
end

%Import manifest and initiate storage variables
mod_manifest = readtable(manifest_file,'ReadRowNames',false);
sample_names = mod_manifest.sample;
num_sample = length(sample_names); 
tax_cell = cell(num_sample,length(tax_file_cell));
total_reads_vec = nan(num_sample,1);

%Loop through all samples and add table data
for i = 1:num_sample
    
    %Get matlab valid name to control for matlab's import behavior
    orig_name = sample_names{i};
    valid_name = matlab.lang.makeValidName(orig_name);
    
    total_reads_vec(i) = total_reads_table{orig_name,'Tot_Samp_Reads'};

    %Add entries to cells from the raw tables
    for j = 1:length(tax_file_cell)
        tax_cell{i,j} = tax_table{j}(:,valid_name);
        phage_cell{i,j} = phage_table{j}(:,valid_name);
        bacteria_cell{i,j} = bacteria_table{j}(:,valid_name);
        eukarya_cell{i,j} = eukarya_table{j}(:,valid_name);
        archaea_cell{i,j} = archaea_table{j}(:,valid_name);
    end

end

%Add total reads to manifest
mod_manifest.total_reads = total_reads_vec;

%Add taxonomic data to the manifest
%All taxa
mod_manifest.read = tax_cell(:,1);
mod_manifest.species = tax_cell(:,2);
%Only phage
mod_manifest.read_phage = phage_cell(:,1);
mod_manifest.species_phage = phage_cell(:,2);
%Only bacteria
mod_manifest.read_bacteria = bacteria_cell(:,1);
mod_manifest.species_bacteria = bacteria_cell(:,2);

%Get the read abundance tables with only species-level entries
levels = {'read_phage','read_bacteria'};
for i = 1:size(mod_manifest,1)
    for j = 1:length(levels)
        species_str_ind = contains(mod_manifest{i,levels{j}}{1}.Properties.RowNames,'species_');

        species_only_read{i,j} = mod_manifest{i,levels{j}}{1}(species_str_ind,:);
    end
end
mod_manifest.species_only_read_phage = species_only_read(:,1);
mod_manifest.species_only_read_bacteria = species_only_read(:,2);

%Add in the tables for other superkingdoms
mod_manifest.species_eukarya = eukarya_cell(:,2);
mod_manifest.species_archaea = archaea_cell(:,2);

%Sort rows by sample name
mod_manifest = sortrows(mod_manifest,'sample');

%Get virulence score predictions and match to table
virscore_table = readtable(virscore_file);
virus_cell = mod_manifest.species_phage{1}.Properties.RowNames;
num_virus = length(virus_cell);
virscore = nan(size(mod_manifest.species_phage{1},1),1);
for i = 1:num_virus
    species_name = strsplit(mod_manifest.species_phage{1}.Properties.RowNames{i},'species_');
    species_name = species_name{2};
    matching_virscore = virscore_table.Var2(strcmp(virscore_table.Var1,species_name));
    if ~isempty(matching_virscore)
        virscore(i) = matching_virscore;
    end
end

%Identify which phage species are unclassified, virulent, and temperate
null_ind = isnan(virscore);
virulent_ind = (virscore > 0.5) & (~null_ind); 
temperate_ind = (virscore <= 0.5) & (~null_ind);

%Get phage and microbe index
species_names = mod_manifest.species{1}.Properties.RowNames;
phage_ind = contains(species_names,sk_vir);
bacteria_ind = contains(species_names,sk_bac);
archaea_ind = contains(species_names,sk_arc);
human_ind = contains(species_names,sp_hum);
all_eukarya_ind = contains(species_names,sk_euk);
eukarya_ind = all_eukarya_ind & (~human_ind);
microbe_ind = bacteria_ind | archaea_ind | eukarya_ind;

%Get the unclassified read index
read_names = mod_manifest.read{1}.Properties.RowNames;
all_root_ind = false(size(read_names));
cell_root_ind = false(size(read_names));
non_cell_root_ind = false(size(read_names));
for i = 1:length(root_str)
   all_root_ind = all_root_ind | strcmp(read_names,root_str{i});
   cell_root_ind = cell_root_ind | strcmp(read_names,cell_root_str{i});
   non_cell_root_ind = non_cell_root_ind | strcmp(read_names,non_cell_root_str{i});
end

%Define diversity functions
H = @(p) nansum(-(p/sum(p)).*log2(p/sum(p)));
richness = @(p) sum(p > 0);
p0 = 1e-3; 
corr_richness = @(p) sum(1 - exp(-p./p0));

%Loop through samples and compute various metrics
for i = 1:size(mod_manifest,1)

    %Fraction of unclassified reads
    read_vec = mod_manifest.read{i}.Variables;
    mod_manifest.unclassified_read_frac(i) = sum(read_vec(all_root_ind)) - ...
        sum(read_vec(cell_root_ind | non_cell_root_ind));
    mod_manifest.cell_read_frac(i) = sum(read_vec(cell_root_ind));
    mod_manifest.viral_read_frac(i) = sum(read_vec(non_cell_root_ind));

    %Ratio of virulent to non-virulent
    phage_vec = mod_manifest.species_phage{i}.Variables;
    mod_manifest.vir_temp_ratio(i) = ...
        sum(phage_vec(virulent_ind))/sum(phage_vec(temperate_ind));
    mod_manifest.null_frac(i) = sum(phage_vec(null_ind))/sum(phage_vec);
    
    %Phage to microbe ratio 
    all_vec = mod_manifest.species{i}.Variables;
    mod_manifest.phage_to_microbe_ratio(i) = ...
        sum(all_vec(phage_ind))/sum(all_vec(microbe_ind));
    mod_manifest.phage_to_bacteria_ratio(i) = ...
        sum(all_vec(phage_ind))/sum(all_vec(bacteria_ind));

    %Total phage, bacteria and microbe frac
    mod_manifest.phage_frac(i) = sum(all_vec(phage_ind));
    mod_manifest.bacteria_frac(i) = sum(all_vec(bacteria_ind));
    mod_manifest.microbe_frac(i) = sum(all_vec(microbe_ind));

    %Diversity metrics
    mod_manifest.phage_H(i) = H(all_vec(phage_ind));
    mod_manifest.bacteria_H(i) = H(all_vec(bacteria_ind));
    mod_manifest.archaea_H(i) = H(all_vec(archaea_ind));
    mod_manifest.eukarya_H(i) = H(all_vec(eukarya_ind));
    mod_manifest.H_diff(i) = mod_manifest.phage_H(i) - mod_manifest.bacteria_H(i);

    %Richness metrics
    mod_manifest.phage_richness(i) = richness(all_vec(phage_ind));
    mod_manifest.bacteria_richness(i) = richness(all_vec(bacteria_ind));
    mod_manifest.archaea_richness(i) = richness(all_vec(archaea_ind));
    mod_manifest.eukarya_richness(i) = richness(all_vec(eukarya_ind));

    %Corrected richness metrics
    mod_manifest.phage_corr_richness(i) = corr_richness(all_vec(phage_ind));
    mod_manifest.bacteria_corr_richness(i) = corr_richness(all_vec(bacteria_ind));
    mod_manifest.archaea_corr_richness(i) = corr_richness(all_vec(archaea_ind));
    mod_manifest.eukarya_corr_richness(i) = corr_richness(all_vec(eukarya_ind));

    %Evenness metrics
    mod_manifest.phage_evenness(i) = ...
        mod_manifest.phage_H(i)/log2(mod_manifest.phage_richness(i));
    mod_manifest.bacteria_evenness(i) = ...
        mod_manifest.bacteria_H(i)/log2(mod_manifest.bacteria_richness(i));
    mod_manifest.archaea_evenness(i) = ...
        mod_manifest.archaea_H(i)/log2(mod_manifest.archaea_richness(i));
    mod_manifest.eukarya_evenness(i) = ...
        mod_manifest.eukarya_H(i)/log2(mod_manifest.eukarya_richness(i));

end


end

