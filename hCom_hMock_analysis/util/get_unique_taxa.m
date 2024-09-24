function unique_taxa = get_unique_taxa(species_cell,taxa_str)

%This function takes in a cell of species names and produces a cell of
%unique taxa names at the specified high taxonomic level (e.g. genus)

%Split taxa by specified taxa string
unique_taxa = cellfun(@(x) strsplit(x,taxa_str),species_cell,'UniformOutput',false);

%Include only cases where higher taxon is found
taxa_len = cellfun(@length,unique_taxa);
unique_taxa = unique_taxa(taxa_len > 1);

%Get only the higher order taxa and find unique set
unique_taxa = cellfun(@(x) strsplit(x{2},'|'),unique_taxa,'UniformOutput',false);
unique_taxa = cellfun(@(x) x{1},unique_taxa,'UniformOutput',false);
unique_taxa = unique(unique_taxa);

end