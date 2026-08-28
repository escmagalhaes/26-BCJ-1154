## TITLE ## 
High MECOM Protein Levels Confer Adverse Prognosis in Acute Myeloid Leukemia Independent of Chromosome 3 Abnormalities

## RUNNING TITLE ##
Wildtype MECOM Protein is adverse in AML

## AUTHORS ##
Steven M. Kornblau [1] *, Eduardo Sabino de Camargo Magalhães [2] *, Brandon D. Brown [3], Samanta Soledad Catueno [3], Yihua Qiu [1], Eitan Kugler [1], Abhishek Maiti [1], Warren C. Fiskus [1], Tapan Kadia [1] 

## AFFILIATIONS ##
[1] Department of Leukemia, The University of Texas M.D. Anderson Cancer Center, Houston, TX, USA 

[2] Department of Ageing Biology/ERIBA, University Medical Center Groningen, Groningen, Netherlands; Glial Cell Biology Laboratory, Biomedical Sciences Institute, Universidade Federal do Rio de Janeiro, Rio de Janeiro, Brazil

[3] Division of Pediatrics, The University of Texas M.D. Anderson Cancer Center, Houston, TX, USA 

*These authors contributed equally to this work

## CORRESPONDING AUTHOR ##
Steven Kornblau, Address 1515 Holcombe Blvd, Box 448, Houston, TX 77030-4009 (skornblau@mdanderson.org)

## ABSTRACT ##
MDS1 and EVI1 Complex locus (MECOM) gene rearrangements (MECOM-R) occur in ~2% of AML, with EVI1 protein overexpression causing aggressive disease. EVI1 overexpression in non-MECOM-R cases also occurs and is prognostically adverse, suggesting that EVI1 protein overexpression could induce similar biological programs without MECOM-R. We measured MECOM protein in 810 untreated AML cases, (MECOM-R=2%), observing low protein in 86%, and high expression in 14%. High protein was equally prognostically adverse with MECOM-R, other Chromosome 3 (Chr3) abnormalities (Chr3-abn) or Chr3-wild-type (WT). Despite high protein expression co-associating with known adverse features, high MECOM was independently prognostic in CoxPH models. High MECOM protein overrides the benefit of mutant CEBPα protein. Protein-protein correlations and networks were similar between MECOM-R, Chr3-WT-high and Chr3-abn-high, with common activated pathways, but also distinct ones. Chr3-WT-high cases also overexpressed MDS1/EVI1, and expressed protein favoring epigenetic deregulation, explaining some observed differences. Collectively, we showed that elevated MECOM protein levels in non-MECOM-R cases are equally adverse to MECOM-R cases, and share similar biological programs, suggesting that therapies targeting MECOM-related proteins may be efficacious beyond those with MECOM-R. Clinically measuring MECOM protein levels could be a useful strategy to stratify patients and guide therapy recommendation.

## Licenses ##

**Code** (files under `R/`) is licensed under [MIT License](LICENSE.md).

**Data** (files under `inputs/`) is licensed under [CC-BY-4.0](DATA_LICENSE.md).

Specific code snippets from third parties are licensed under [MIT License](THIRD_PARTY_LICENSES.md)

## INSTRUCTIONS ##
# Requirements #
R version 4.5.2

RStudio (recommended)

Python version 3.10.11 (if using Synapse data)

Cytoscape software version 3.10.4

Synapse account (recommended)

Synapse Personal Access Token (PAT) (recommended)

Internet connection 

# Getting Started #
1) Clone or download this repository
2) Open the project file (.Rproj) located inside main project folder
3) Open and run the script "load_all.R" located inside folder "R/setup/"

Note: Synapse data is not necessary to reproduce Figures/Tables in this manuscript. 
This data was only used for exploratory purposes. 

If assessing Synapse data is necessary, visit this website to create an account and a Personal Access Tokens (PAT):
https://docs.synapse.org/synapse-docs/managing-your-account#ManagingYourAccount-PersonalAccessTokens

Add the PAT to `config.R` script as and object called `.synapse_token` and restart R to run scripts.

# To Reproduce All Figures and Tables in sequence #
Run "run_all.R" located inside main project folder 
(REMEMBER TO OPEN CYTOSCAPE SOFTWARE BEFORE STARTING!)

# To Reproduce specific Figures or Tables #
1) Open and run the script "load_all.R" located inside folder "R/setup/"
2) Open and run the script "00_preprocessing_inputs.R" located inside folder "R/analyses/" 
3) Run the required analysis script(s) located inside folder "R/analyses/" according to the table below. 
4) Run the desired workflow script from located inside folder "R/workflows/" according to the table below
Note: if multiple analysis scripts are needed, run in numerical order (e.g. "01_" before "02_", etc.).

Figure/Table		 | R/workflows script		 | R/analyses script	
-------------------------|-------------------------------|------------------------------------------
Figure1			 | figure_metanalysis.R 	   | 01_online_datasets	
Figure2			 | figure_dens_km_plots_os.R | 02_KM_plots	
Figure3			 | figure_mv_models.R		     | 04_UV_MV_models	
Figure4			 | figure_prog_mods.R		     | 06_prognostic_models	
Figure5			 | figure_scatter_plots.R	   | 07_correlation_analysis, 08_scatter_plots	
Figure6			 | figure_paths_networks.R	 | 07_correlation_analysis, 09_processing_correlations_output, 10_pathway_analysis	
Table1			 | table_descr_sig.R		     | 03_descriptive_tables	
Supplemental_Figure_S1	 | figure_km_plots_cohorts.R	 | 01_online_datasets	
Supplemental_Figure_S2	 | figure_uv_models_cohorts.R	 | 01_online_datasets	
Supplemental_Figure_S3	 | figure_km_plots_rd.R 	 | 02_KM_plots	
Supplemental_Figure_S4	 | figure_km_plots_qt.R		 | 02_KM_plots	
Supplemental_Figure_S5	 | figure_boxplots.R		 | 05_boxplots	
Supplemental_Figure_S6	 | figure_scatter_plots_extra.R	 | 07_correlation_analysis, 08_scatter_plots	
Supplemental_Table_S1	 | N/A				 | N/A	
Supplemental_Table_S2	 | dataset_multi_cohort.R	 | 01_online_datasets	
Supplemental_Table_S3	 | table_descr_ext.R		 | 03_descriptive_tables	
Supplemental_Table_S4	 | table_descr_alternate.R 	 | 03_descriptive_tables	
Supplemental_Table_S5	 | dataset_boxplot_tests.R	 | 05_boxplots	
Supplemental_Table_S6	 | table_uv_models.R		 | 04_UV_MV_models	
Supplemental_Table_S7	 | dataset_prog_mods_bxp_tests.R | 06_prognostic_models	
Supplemental_Table_S8	 | dataset_corr_data.R		 | 07_correlation_analysis	
Supplemental_Table_S9	 | dataset_combined_paths.R	 | 07_correlation_analysis, 09_processing_correlations_output, 10_pathway_analysis
	
Note: Supplemental_Figure_S1 and Figure 6E were generated externally  

## MAPPING OF INPUTS ##
Dataset    | Meaning
-----------|----------------------------------------------------------------------------------------
data1.csv  | table with clinical variables and protein expression data (RPPA) of AML patients
data2.csv  | table with metadata from clinical variables
data3.csv  | table with protein expression values of CD34+ normal bone marrow-derived cells
data4.csv  | table with Protein Functional Groups (PFG) assignment for each protein
data5.csv  | table with MECOM mRNA expression data, survival data and cytogenetic characteristics from UTMDACC samples (N=20)
data6.csv  | table with MECOM mRNA expression data, survival data and cytogenetic characteristics from Bottomly, D. et al. (BeatAML) cohort
data7.csv  | table with MECOM mRNA expression data, survival data and cytogenetic characteristics from TCGA et al. (LAML) cohort
data8.csv  | table with MECOM mRNA expression data, protein expression data (Mass Spectrometry), survival data and cytogenetic characteristics from Kramer, M. et al. cohort
data9.csv  | table with MECOM protein expression data (Mass Spectrometry), survival data and cytogenetic characteristics from Jayavelu, A. et al. cohort
data10.csv | table with MECOM protein expression data (Mass Spectrometry), survival data and cytogenetic characteristics from UTMDACC samples (N=47)

## MAPPING OF DATASET VARIABLE TO CORRESPONDING LABEL ##
Variable Name        | Label
---------------------|-------------------------------------------------
group                | Group 
upi                  | Unique Patient Identifier Number 
status               | Deceased 
relapse              | Relapse 
event                | Event 
surv_time            | Overall Survival (years) 
rem_time             | Remission Duration (years) 
event_time           | Event-free Survival (years) 
age                  | Age at diagnosis (years) 
age_gp               | Age group (years) 
age_40_minus         | Age < 40 (years) 
age_41_55            | Age 41 - 55 (years) 
age_56_70            | Age 56 - 70 (years) 
age_70_plus          | Age > 70 (years) 
gender               | Gender 
race                 | Race 
race_white           | White (race) 
race_black           | Black (race) 
race_hisp            | Hispanic (race) 
race_asian           | Asian (race) 
race_other           | Other (race) 
wbc                  | White Blood Cell count (K/uL) 
blasts               | Blasts (%) 
hgb                  | Hemoglobin (g/dl) 
plt                  | Platelets (K/uL) 
fab2                 | FAB classification 
aml_gp               | Primary or Secondary AML 
tx_aml               | Therapy-related AML 
pre_mlgn             | Prior Malignancy 
ps                   | Performance status (PS) 
eln_risk	           | ELN risk (2022)
chem2                | Therapeutic regimen 
chem3                | Therapeutic regimen expanded 
tx1_resp2            | Initial Therapy Response 
arac                 | AraC-based therapy 
hma                  | HMA-based therapy 
vtx                  | Venetoclax therapy 
flt3i                | FLT3 inhibitor therapy 
idhi                 | IDH inhibitor therapy 
hsct                 | Stem-cell Transplant 
karyotype            | Cytogenetics
cyto_risk            | Cytogenetic Risk 
cyto_unfav           | Unfavorable Cytogenetics 
cyto_cat             | Cytogenetic Category 
complex_kar          | Complex Karyotype 
diploid              | Diploid Karyotype 
del5_5q              | -5/5q- 
del7_7q              | -7/7q- 
inv16                | Inv 16 
del12                | Del 12 
MECOM_R              | MECOM-R
chr3_abn             | Other Chromosome 3 abnormalities
chr3_status_cat      | Chromosome 3 alteration type
chr3_status          | Chromosome 3 status 
chr3_status_simple   | Chromosome 3 status simplified 
t_8_21               | t(8;21) 
t_15_17              | t(15;17) 
t_11q23              | t(11q23) 
t_9_11               | t(9;11) 
t_6_9                | t(6;9) 
trisomy8             | Trisomy 8 
trisomy6             | Trisomy 6 
t_9_22_Ph_pos        | t(9;22)Ph+ 
mut_asxl1            | ASXL1 Mutation 
mut_bcor             | BCOR Mutation 
mut_cebpa            | CEBPA Mutation 
mut_cebpa_db         | CEBPA double Mutation 
mut_dnmt3a           | DNMT3a Mutation 
mut_ezh2             | EZH2 Mutation 
mut_flt3             | FLT3 Mutation 
mut_flt3_itd         | FLT3-ITD Mutation 
mut_flt3_itd_high_ar | FLT3-ITD Mutation with high allelic ratio 
mut_flt3_d835        | FLT3-TKD Mutation 
mut_gata2            | GATA2 Mutation 
mut_idh              | IDH Mutation 
mut_idh1             | IDH1 Mutation 
mut_idh2             | IDH2 Mutation 
mut_jak2             | JAK2 Mutation 
mut_kit              | KIT Mutation 
mut_kmt2a            | KMT2A Mutation 
mut_mll              | MLL Mutation 
mut_npm1             | NPM1 Mutation 
mut_phf6             | PHF6 Mutation 
mut_ptpn11           | PTPN11 Mutation 
mut_rad21            | RAD21 Mutation 
mut_ras              | RAS Mutation 
mut_nras             | KRAS Mutation 
mut_kras             | NRAS Mutation 
mut_runx1            | RUNX1 Mutation 
mut_sf3b1            | SF3B1 Mutation 
mut_smc1a            | SMC1A Mutation 
mut_smc3             | SMC3 Mutation 
mut_srsf2            | SRSF2 Mutation 
mut_stag2            | STAG2 Mutation 
mut_tet2             | TET2 Mutation 
mut_tp53             | TP53 Mutation 
mut_u2af1            | U2AF1 Mutation 
mut_wt1              | WT1 Mutation 
mut_zrsr2            | ZRSR2 Mutation

## PROJECT STRUCTURE ##
```
Main project folder/
│
├── README.md│
├── .gitignore
├── renv.lock
├── .Rprofile
│
├── run_all.R
│
├── bio_cache/
│   └── .gitkeep
│
├── renv/
│   ├── activate.R
│   ├── settings.json
│   └── .gitignore
│
├── inputs/
│   ├── data1.csv
│   ├── data2.csv
│   ├── data3.csv
│   ├── data4.csv
│   ├── data5.csv
│   ├── data6.csv
│   ├── data7.csv
│   ├── data8.csv
│   ├── data9.csv
│   └── data10.csv
│
├── R/
│   ├── setup/
│   │     ├── config.R
│   │     ├── load_all.R
│   │     ├── packages.R
│   │     └── renv.R
│   │  
│   ├── helper_functions/
│   │   └── helper_functions.R
│   │
│   ├── analyses/
│   │    ├── 00_preprocessing_inputs.R
│   │    ├── 01_online_datasets.R
│   │    ├── 02_KM_plots.R
│   │    ├── 03_descriptive_tables.R
│   │    ├── 04_UV_MV_models.R
│   │    ├── 05_boxplots.R
│   │    ├── 06_prognostic_models.R
│   │    ├── 07_correlation_analysis.R
│   │    ├── 08_scatter_plots.R
│   │    ├── 09_processing_correlations_output.R
│   │    └── 10_pathway_analysis.R
│   │
│   └── workflows/
│         ├── figure_metanalysis.R 
│         ├── figure_dens_km_plots_os.R
│         ├── figure_mv_models.R
│         ├── figure_prog_mods.R
│         ├── figure_scatter_plots.R
│         ├── figure_paths_networks.R
│         ├── table_descr_sig.R
│         ├── figure_km_plots_cohorts.R
│         ├── figure_uv_models_cohorts.R
│         ├── figure_km_plots_rd.R 
│         ├── figure_km_plots_qt.R
│         ├── figure_boxplots.R
│         ├── figure_scatter_plots_extra.R
│         ├── dataset_multi_cohort.R
│         ├── table_descr_ext.R
│         ├── table_descr_alternate.R 
│         ├── dataset_boxplot_tests.R
│         ├── table_uv_models.R
│         ├── dataset_prog_mods_bxp_tests.R
│         ├── dataset_corr_data.R
│         └── dataset_combined_paths.R
│
└── output/
     ├── Analyses/
     │     └── .gitkeep
     │
     └── Exported/
           └── .gitkeep
```
  
