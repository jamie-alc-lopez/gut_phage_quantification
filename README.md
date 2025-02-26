# gut_phage_quantification
Code repository for Lopez et al. 2024 "Abundance measurements reveal the balance between lysis and lysogeny in the human gut microbiome". Below, we describe scripts within the repository by folder. 

### common_util
Utility scripts and files used by multiple analyses.

-`generate_default_font_config.m`: generates `font_config.mat`, which contains default formatting settings for generating figures

-`get_ictv_taxonomies.m`: uses `metadata/votous_metadata.tsv` to get ICTV taxonomies given a list of UHGV vOTUs

-`import_phanta_results.m`: imports output from Phanta into a MATLAB table

-`Violinplot-Matlab-master`: folder containing third-party MATLAB code for making violin plots (Bastian Bechtold, 10.5281/zenodo.4559847). Slightly modified for this manuscript's plots

### virome_abundance_meta_analysis
Code for generating Fig. 1 and associated SI figures

-`Fig_VLP_bulk_absolute_abundance.m`: generates Fig. 1CD and the associated SI figures using data files `cohort_analyses/process_shkoporov_2019/shkoporov_2019_absolute_data.mat` and `cohort_analyses/process_liang_2020/liang_2020_absolute_data.mat`

-`Fig_abundance_quant.m`: generates Fig. 1B using data in the `final_aggregate_tables` folder. Also generates a table of summary statistics

-`standardize_study_data.m`: for studies where the needed quantification data can be immediately taken from SI file or manuscript figure, this script formats quantification data into a standardized format. It takes in files from `study_data` folder and outputs standardized tables in `final_aggregate_tables`

-`SI_fig_no_MDA_VLP_bulk_relative_abundance.m`: Generates Fig. S3. 

### hCom_hMock_analysis
Code for generating Fig. 3-4 and associated SI figures.

-`Fig_hCom2_human_comparison.m`: generates Fig. 4 and S4

-`Fig_hCom2_reconstruction.m`: generates Fig. 3 and S3

-`SI_fig_phatyp_virulence.m`: generates Fig. S5

-`SI_table_hCom2_human_taxa_comparison.m`: analyzes shared taxa between hCom2 mouse stool and the Yachida et al. stool

-`util/compute_relative_stats.m`: compute statistics of a value x relative to a dataset V

-`util/get_unique_taxa.m`: given a list of species produces list of unique higher order taxa

-`util/merge_phanta_table.m`: merges two Phanta tables

-`util/wrapped_*_phanta_import.m`: dataset-specific wrappers for `import_phanta_results.m`

-`generate_grinder_specs/generate_grinder_specs.m`: generates the abundance files used by Grinder to generate hMock


### induction_estimate_analysis
Code for generating Fig. 2 and associated SI figures

-`Fig_induction_fun.m`: generates Fig. 2C

### cohort_analyses
This folder contains analyses and secondary quantifications of different gut virome studies. All files of the name `wrapped_*_phanta_import.m` are dataset-specific wrappers for `import_phanta_results.m`.

-`process_kieft_2021/analysis_kieft_R.m`: generates Table 2 and computes the median number of induced phage per sample

-`process_liang_2020/analysis_liang_VLP_stool.m`: analyzes the Phanta quantification of Liang et al. 2020. Generates `virome_abundance_meta_analysis/final_aggregate_tables/liang_stepwise_2020_phanta_aggregated.xlsx`, `liang_2020_absolute_data.mat`, and computes VMRs of the samples to assess bacterial contamination. 

-`process_shkoporov_2019/analysis_shkoporov_VLP_stool.m`: analyzes the Phanta quantification of Shkoporov et al. 2019. Generates `virome_abundance_meta_analysis/final_aggregate_tables/shkoporov_2019_phanta_aggregated.xlsx`, `shkoporov_2019_absolute_data.mat`, and computes VMRs of the samples to assess bacterial contamination. 

-`process_yachida_2019/analysis_yachida_VLP_stool.m`: analyzes the Phanta quantification of Yachida et al. 2019. Generates `virome_abundance_meta_analysis/final_aggregate_tables/yachida_2019_phanta_aggregated.xlsx`.

-`process_garmaeva_2024/analysis_garmaeva_VLP_stool.m`: analyzes the Phanta quantification of Garmaeva et al. 2024. Generates `garmaeva_2024_*_absolute_data.mat`, and computes VMRs of the samples to assess bacterial contamination.

-`process_shalon_culver_2024/analysis_shalon_culver_R.m`: analyzes phage induction in the Shalon and Culver et al. 2024 data. Generates Table 3. 