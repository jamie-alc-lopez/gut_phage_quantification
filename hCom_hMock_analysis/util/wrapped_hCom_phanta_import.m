function manifest = wrapped_hCom_phanta_import(folder,vir_type)

%This script is a wrapper of import_phanta_results for the hCom dataset

if ~exist('vir_type')
    vir_type = 'bacphlip';
end

manifest_file = 'hCom_hMock_manifests/hCom_manifest.xlsx';

tax_file_cell = {'relative_read_abundance.txt','relative_taxonomic_abundance.txt'};

tax_file_cell = strcat(folder,tax_file_cell);

if strcmp(vir_type,'bacphlip')
    virscore_file = [folder,'/species_name_to_vir_score.txt'];
elseif strcmp(vir_type,'phatyp')
    virscore_file = [folder,'/species_name_to_vir_score_phatyp.txt'];
end

total_reads_file = [folder,'/total_reads.tsv'];

manifest = import_phanta_results(tax_file_cell,manifest_file,...
    virscore_file,total_reads_file);


end

