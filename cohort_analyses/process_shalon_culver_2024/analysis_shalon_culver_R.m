%% This script analyzes the individual prophage-host data from Shalon and 
%Culver et al. to extract mean values of R.

% Import and concatenate data from all of the individual subfolders

clear;clc

data_folder = 'shalon_culver_data/propagate_subset_output/';
subfolders = dir(data_folder);
master_table = [];
VariableTypes = ["char","char","char","double","double","double","double",...
    "double","double","double","double","double","double","double",...
    "double"];

%Loop though all of the subfolders and add them to the mastert able. 
for i = 1:length(subfolders)

    samp_name = subfolders(i).name;

    if ~strcmp(samp_name(1),'.')

        file_i = [data_folder,samp_name,'/',samp_name,'.tsv'];

        opts = detectImportOptions(file_i,'FileType','text',...
            'ReadVariableNames',true,'Delimiter','\t');
        opts.VariableTypes = VariableTypes;

        table_i = readtable(file_i,opts);

        table_i.sample(:) = {samp_name};

        master_table = [master_table; table_i];

    end
end


%Perform additional QC
QC_ind = (master_table.host_median_cov>1) & ...
    (master_table.prophage_cov_breadth > 0.5) & ...
    (master_table.prophage_median_cov > 1) & ...
    ~strcmp(master_table.active,'not present');
master_table = master_table(QC_ind,:);


%% Load in metadata and merge with propagate output

metadata_file = 'shalon_culver_data/sample_metadata_matlab_export.xlsx';
opts = detectImportOptions(metadata_file,...
    'ReadVariableNames',true,'VariableNamesRange','A1','DataRange','A2');

metadata = readtable(metadata_file,opts);
metadata.full_type = strcat(metadata.Type,metadata.Putative_location);

%Fix metadata labeling issues
metadata.full_type = strrep(metadata.full_type,'stoolStool','Stool');
metadata.full_type = strrep(metadata.full_type,'stool','Stool');

%Add sample types to master table
metadata.Properties.RowNames = metadata.Abrreviation;
master_table.sample_type = metadata.full_type(master_table.sample);


%% Compute summary statistics

%Get clipped R
master_table.clipped_R = master_table.prophage_host_ratio;
master_table.clipped_R(master_table.clipped_R < 1) = 1;

%Splitapply to compute stats
[G,id] = findgroups(master_table.sample_type);
mean_R = splitapply(@mean,master_table.prophage_host_ratio,G);
median_R = splitapply(@median,master_table.prophage_host_ratio,G);
mean_clipped_R = splitapply(@mean,master_table.clipped_R,G);
n_pairs = splitapply(@length,master_table.clipped_R,G);
n_samples = splitapply(@(x) length(unique(x)),master_table.sample,G);

%Compute mean induction using model function
B = 10; %burst size
delta = 1; %dilution rate
gamma = delta; %lysis time
xi_fun = @(R_minus_1) (R_minus_1).*(gamma + delta)./(B-1 - R_minus_1); %From mgx fold changes
mean_induction = xi_fun(mean_clipped_R-1);

%Build and export summary table
summary_array = round([n_samples,n_pairs,mean_R,median_R,...
    mean_clipped_R,mean_induction],3,'significant');
summary_table = array2table(summary_array,...
    'RowNames',id,'VariableNames',{'n_samples','n_pairs','base_mean_R',...
    'median_R','clipped_mean_R','mean_induction_rate'});

writetable(summary_table,'shalon_culver_induction_table.xlsx',...
    'WriteRowNames',true,'WriteVariableNames',true)
