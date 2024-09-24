%% This portion of the script takes in the NinjaMap outputs and the genome
% lengths and produces the abundance profiles for Grinder

clear;clc

%Load in hCom2 genome lengths
genome_frag_lengths = readtable('hCom2_metadata/all_hCom2_genome_lengths.txt');
genome_frag_lengths.Properties.VariableNames = {'genome_frag','length'};

%Load in NinjaMap outputs
ninjamap_output = readtable('hCom2_metadata/8a_backfill2_20200722_MBFv2.long_Prelim_Ninjamap.csv');
sample_names = unique(ninjamap_output.sample_id);
strain_names = unique(ninjamap_output.Strain_Name);

%Manually curate strain names to agree between final genome list and NinjaMap
orig_strain_names = strain_names;
strain_names = cellfun(@(x) strrep(x,'-MAF-1',''),strain_names,'UniformOutput',false);
strain_names{39} = strain_names{39}(1:end-1);
strain_names = cellfun(@(x) strrep(x,'-MAF-','-'),strain_names,'UniformOutput',false);
strain_names = cellfun(@(x) strrep(x,'-NJ35',''),strain_names,'UniformOutput',false);
strain_names = cellfun(@(x) strrep(x,'-NJ26',''),strain_names,'UniformOutput',false);

%Generate strain information, link strain fragments to their strain genome
strain_frag_inds = false(length(strain_names),length(genome_frag_lengths.length));
total_lengths = zeros(length(strain_names),1);
for i = 1:length(strain_names)
    strain_i = strain_names{i};
    strain_frag_inds(i,:) = contains(genome_frag_lengths.genome_frag,strain_i);
    total_lengths(i) = sum(genome_frag_lengths.length(strain_frag_inds(i,:)));
end

%Note: Grinder wants relative taxonomic abundances (e.g.
%cov_i/sum_j(cov_j)) as its input!

%Loop through samples and create corresponding Grinder relative abundance
%files
for i = 1:length(sample_names)

    %Get sample name
    sample_i = sample_names{i};

    %Get subtable of ninjamap measurements for that sample
    ninjamap_subtable = ninjamap_output(strcmp(ninjamap_output.sample_id,sample_i),:);

    %Initiate storage cells/vectors for this samples
    genome_frag_cell = {};
    frag_length_vec = [];
    read_abun_vec = [];
    cov_vec = [];

    %Loop through each strain detected in the sample
    for j = 1:size(ninjamap_subtable,1)

        %Get information on strain j
        strain_j = ninjamap_subtable.Strain_Name{j};
        strain_j_cov = ninjamap_subtable.Coverage_Depth(j);
        strain_j_ind = find(strcmp(orig_strain_names,strain_j));
        strain_j_genome_name = strain_names{strain_j_ind};
        strain_j_genome_size = total_lengths(strain_j_ind);

        %Get strain j genome fragment names and read abundance of 
        % each genome fragment
        strain_j_frags = genome_frag_lengths.genome_frag(strain_frag_inds(strain_j_ind,:));

        %Get the lengths of the fragments
        strain_j_frag_lengths = genome_frag_lengths.length(strain_frag_inds(strain_j_ind,:));

        %Get the relative read abundance of the fragments. If there are
        %multiple fragments, we assume the individual relative read abundance 
        %is proportional to the fragment length
        strain_j_frag_read_abun = strain_j_frag_lengths.*ninjamap_subtable.Read_Fraction(j)./strain_j_genome_size;

        %Store strain j info
        genome_frag_cell = [genome_frag_cell; strain_j_frags];
        frag_length_vec = [frag_length_vec; strain_j_frag_lengths];
        read_abun_vec = [read_abun_vec; strain_j_frag_read_abun];
        cov_vec = [cov_vec; repmat(strain_j_cov,size(strain_j_frags))];

    end

    %Normalize read abundance to fraction
    read_abun_vec = read_abun_vec/sum(read_abun_vec);

    %Compute taxonomic abundance starting from relative read abundances
    taxon_abun_vec_from_reads = (read_abun_vec./frag_length_vec)...
        ./sum((read_abun_vec./frag_length_vec));

    %Compute taxonomic abundance directly from NinjaMap coverages 
    taxon_abun_vec_from_cov = cov_vec/sum(cov_vec);

    %Save to abundance file using the taxon abundances computed from
    %coverages, as these are the most direct (the read and coverage-derived
    %values are very similar)
    sample_i_table = table(genome_frag_cell,taxon_abun_vec_from_cov);
    filename = ['grinder_abundance_files/',sample_i,'_grinder_abundances.txt'];
    writetable(sample_i_table,filename,'Delimiter',' ','WriteVariableNames',false);

end
