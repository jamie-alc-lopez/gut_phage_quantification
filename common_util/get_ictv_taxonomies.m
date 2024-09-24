function ictv_taxonomies = get_ictv_taxonomies(species_cell)

%This script takes in a list of UHGV phage taxa and produces to corresponding set
%of ictv taxonomies. This is used for identifying certain taxa (e.g.
%microviridae)

%Get vOTUs from the uhgv taxonomy
votus = cellfun(@(x) strsplit(x,'species_'),species_cell,'UniformOutput',false);
votus = cellfun(@(x) x{2},votus,'UniformOutput',false);

%Load UHGV vOTU metadata
fun_path = mfilename('fullpath');
uhgv_metadata_loc = strrep(fun_path,'get_ictv_taxonomies','metadata/votus_metadata.tsv');
metadata = readtable(uhgv_metadata_loc,'FileType','text');

%Make the table indexed by the votu names
metadata.Properties.RowNames = metadata.uhgv_votu;

%Get corresponding ictv taxonomies
ictv_taxonomies = metadata{votus,'ictv_taxonomy'};

end