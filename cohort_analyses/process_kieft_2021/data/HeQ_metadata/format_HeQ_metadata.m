clear;clc

%This script formats the HeQ metadata. The only available metadata was in a
%rough form on the Gigascience DB, so it must be manually reformatted

%Load in raw metadata from Gigascience repo
raw_HeQ_meta = readtable("raw_HeQ_gigascience_repo_metadata.xlsx");

sample_ind = contains(raw_HeQ_meta.SampleID,'ERS');
sample_positions = find(sample_ind);

%Loop through and find the disease status information corresponding to each
%sample
for i = 1:length(sample_positions)

    pos_i = sample_positions(i);
    
    if i < length(sample_positions)
    final_i = sample_positions(i+1) - 1;
    else
        final_i = length(sample_ind);
    end
       
    attributes_i = raw_HeQ_meta{pos_i:final_i,'SampleAttributes'};

    disease_status_i = attributes_i(contains(attributes_i,'Disease status'));
    age_i = attributes_i(contains(attributes_i,'Age:'));
    age_i = str2double(age_i{1}(5:end));

    raw_HeQ_meta(pos_i,'disease_status') = disease_status_i; 
    raw_HeQ_meta{pos_i,'age'} = age_i;

end

HeQ_metadata1 = raw_HeQ_meta(sample_ind,:);
HeQ_metadata1.Properties.RowNames = HeQ_metadata1.SampleID;

%Now load in SRA data to link disease status to read run ID
HeQ_metadata2 = readtable('filereport_read_run_PRJEB15371_tsv.txt',...
    'FileType','text','ReadRowNames',true,'Delimiter','\t');
HeQ_metadata2.disease_status = ...
    HeQ_metadata1{HeQ_metadata2.secondary_sample_accession,'disease_status'};
HeQ_metadata2.age = ...
    HeQ_metadata1{HeQ_metadata2.secondary_sample_accession,'age'};

%Save final table
writetable(HeQ_metadata2,'processed_HeQ_metadata.xlsx',...
    'WriteRowNames',true,'WriteVariableNames',true)
