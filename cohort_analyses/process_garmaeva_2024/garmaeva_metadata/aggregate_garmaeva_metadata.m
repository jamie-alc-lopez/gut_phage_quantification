%% Make the Garmaeva et al. Phanta samp table using a list of sample files

clear;clc

%Load file list
file_list = readtable('garmaeva_sample_list.txt','Delimiter','\t',...
    'ReadVariableNames',false);

%Construct table columns
read1_ind = contains(file_list.Var1,'pair_1.fastq.gz');
read1_cell = file_list.Var1(read1_ind);
read2_cell = cell(size(read1_cell));
EGA_id_1_cell = read2_cell;
EGA_id_2_cell = read2_cell;
sample_type_cell = read2_cell;
sample_name_cell = read2_cell;


%Loop through and make the sample table entries
loc_prefix = '/oak/stanford/groups/kchuang/jglopez/garmaeva_2024_phanta_analysis/';
for i = 1:length(read1_cell)

    sample_i_1 = read1_cell{i};

    split_str1 = strsplit(sample_i_1,'/');
    EGA_id_1_cell{i} = split_str1{2};
    sample_type_cell{i} = strrep(split_str1{1},'_fastq','');
    sample_name_cell{i} = strrep(split_str1{3},'_kneaddata_cleaned_pair_1.fastq.gz','');

    sample_2_str = ['/',sample_name_cell{i},'_kneaddata_cleaned_pair_2.fastq.gz'];
    sample_i_2 = file_list.Var1{contains(file_list.Var1,sample_2_str)};
    read2_cell{i} = [loc_prefix,sample_i_2];

    split_str2 = strsplit(sample_i_2,'/');
    EGA_id_2_cell{i} = split_str2{2};

    read1_cell{i} = [loc_prefix,sample_i_1];

end


%Generate and save table
phanta_samp_table = table(sample_name_cell,read1_cell,read2_cell,...
    'VariableNames',{'sample','read1','read2'});
writetable(phanta_samp_table,'garmaeva_2024_samp.txt','Delimiter','\t',...
    'WriteVariableNames',false,'WriteRowNames',false);


%% Make the import manifest for processing Phanta results

%Load in metadata pulled from the EGA
bulk_metadata = readtable('garmaeva_bulk_metadata.csv','Delimiter',',');
bulk_metadata.sample_type(:) = {'bulk'};
vlp_metadata = readtable('garmaeva_VLP_metadata.csv','Delimiter',',');
vlp_metadata.sample_type(:) = {'vlp'};

%Select which variables to keep from the metadata
saved_variables = {'accession_id','alias','subject_id','phenotype',...
    'submission_accession_id','sample_type','extra_attributes'};
all_metadata = [bulk_metadata(:,saved_variables);vlp_metadata(:,saved_variables)];

%The timepoint metadata is in a strangely formatted "extra attributes"
%field, this code manually extracts it row-by-row
for i = 1:size(all_metadata,1)

extra_att = strsplit(all_metadata.extra_attributes{i},'},');
timepoint_data = extra_att{contains(extra_att,'Timepoint')};
timepoint_value = strsplit(timepoint_data,'"value":');
timepoint_value = strrep(timepoint_value{2},'"','');
timepoint_value = strrep(timepoint_value,' ','');

all_metadata.timepoint(i) = {timepoint_value};

end

%Fix labeling issues in the sample alias field
fixed_alias = strrep(all_metadata.alias,'VL_','_VL_');
fixed_alias = strrep(fixed_alias,'LN','LN_');
fixed_alias = strrep(fixed_alias,'__','_');
all_metadata.Properties.RowNames = fixed_alias;

%Match the EGA metadata to the files actually downloaded from EGA
manifest = phanta_samp_table; 
manifest.Properties.RowNames = manifest.sample;
matching_manifest = manifest(intersect(all_metadata.Properties.RowNames,manifest.sample),:);
matching_manifest.phenotype = all_metadata.phenotype(matching_manifest.sample);
matching_manifest.subject_id = all_metadata.subject_id(matching_manifest.sample);
matching_manifest.sample_type = all_metadata.sample_type(matching_manifest.sample);
matching_manifest.timepoint = all_metadata.timepoint(matching_manifest.sample);

%Save import manifest
writetable(matching_manifest,'../garmaeva_2024_manifest.xlsx','WriteRowNames',false);