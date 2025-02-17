clear;clc

%This script processes the IjazUZ metadata into a single table

%Load in read archive metadata
IjazUZ_metadata1 = readtable('filereport_read_run_PRJEB18780_tsv.txt',...
    'FileType','text','Delimiter','\t','ReadRowNames',true);
temp_names = cellfun(@(x) strsplit(x,'_'), IjazUZ_metadata1.sample_alias,'UniformOutput',false);
IjazUZ_metadata1.sample = cellfun(@(x) x{1},temp_names,'UniformOutput',false);

%Load in sample data from Ijaz et al. 2017
IjazUZ_metadata2 = readtable('pone.0172605.s005.csv',...
    'FileType','text','Delimiter',',','ReadRowNames',true);

%Add disease status to read archive metadata
IjazUZ_metadata1.disease_status = ...
    IjazUZ_metadata2.SampleTypeExplained(IjazUZ_metadata1.sample);
IjazUZ_metadata1.age = ...
    IjazUZ_metadata2.Age(IjazUZ_metadata1.sample);

%Save final table
writetable(IjazUZ_metadata1,'processed_IjazUZ_metadata.xlsx',...
    'WriteRowNames',true,'WriteVariableNames',true)