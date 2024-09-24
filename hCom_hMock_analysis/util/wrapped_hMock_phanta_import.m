function manifest = wrapped_hMock_phanta_import(folder)

%This script is a wrapper of import_phanta_results for the hMock dataset

manifest_file = 'hCom_hMock_manifests/hMock_manifest.xlsx';

tax_file_cell = {'relative_read_abundance.txt','relative_taxonomic_abundance.txt'};
tax_file_cell = strcat(folder,tax_file_cell);

virscore_file = [folder,'/species_name_to_vir_score.txt'];

total_reads_file = [folder,'/total_reads.tsv'];

manifest = import_phanta_results(tax_file_cell,manifest_file,...
    virscore_file,total_reads_file);


end

