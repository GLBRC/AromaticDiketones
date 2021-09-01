# Aromatic dimer dehydrogenases from Novosphingobium aromaticivorans reduce monoaromatic diketones

Alexandra M. Linz, Yanjun Ma, Jose M. Perez, Kevin S. Myers, Wayne S. Kontur, Daniel R. Noguera, Timothy J. Donohue

DOE Great Lakes Bioenergy Research Center, Madison, WI, USA
Wisconsin Energy Institute, University of Wisconsin-Madison, Madison, WI, USA
Department of Civil and Environmental Engineering, University of Wisconsin-Madison, Madison, WI, USA.
Department of Bacteriology, University of Wisconsin-Madison, Madison, WI, USA

Direct correspondence to: Timothy J. Donohue, email: tdonohue@bact.wisc.edu, phone: 6082624663

## About this page

In this study, we describe the ability of Novosphingobium aromaticivorans, a sphingomonad with a versatile aromatic metabolism, to process aromatic diketones, a potential substrate resulting from lignin depolymerization. We used RNA-Seq to generate a hypothesis about which genes encode the enzymes responsible for G-diketone degradation, then confirmed this hypothesis via in vitro enzyme assays.

Here, you can find data, scripts, and manuscript materials related to this study. If you have any questions, please do not hesitate to contact us.

## Repo Structure

- Enzyme_kinetics (Data and calculations for Michaelis-Menten kinetics)
  - compiled_kinetics_2021-01-12.xlsx (Measured slopes of enzyme oxidation/reduction and results)
  - ligLND_enzymes_2020-07-09.csv (HPLC measurements of diketone degradation)
  - GP-HPLC-samplekey_2021-07-18.csv (HPLC measurements of GP-1 degradation)
  - plot_kinetics_curves.R (Script for generating Figure S5)
  - mm-kinetics_curve_fitting.R (Script for fitting Michaelis-Menten co-efficients based on slopes)
- Growth_curves (Growth of -N. aromaticivorans- on G-diketone
  - Diketone_growth_and_HPLC.csv (Detection of substrate and intermediates during growth)
  - ligLNDO-deletion_growth_curves_cleaned.csv (Growth of deletion mutants on glucose and G-diketone)
- Manuscript_2021-06-03.zip (Materials initially submitted for publication)
- Manuscript_revised (Materials for resubmission after initial review)
- Plots (Intermediate plots generated for figures)
- RNA-Seq
  - Code (Scripts used to process RNA-Seq data)
    - tidy_data.R (Import RNA-Seq data into R datasets)
    - figure_list.R (Script for taking all these datasets and intermediate plots and generating manuscript figures) (Note there was some further processing in Adobe Illustrator in some cases)
  - Data (.csv or .xlsx to open in Excel, .rds to import into R)
    - AromaticDiketones_compiled_RPKM_data.xlsx (compiled RPKM, differential expression, and gene info. This is the file available through GEO.)
    - AromaticDiketones_GEO_submission.xlsx (metadata file for GEO submission)
    - Differential_expression_testing_results (log fold change and FDR-corrected p-values for gene expression comparisons between substrates)
    - Gene_information (locus tags, gene location, and annotation information for -N. aromaticivorans-)
    - gene_names, gene_names2locus_tags (List of gene names and locus tags of published aromatic genes)
    - RNA-Seq_HPLC_results (HPLC data of aromatics detected in supernatant of cultures harvested for RNA)
    - RPKM_normalized_read_counts (Read counts normalized by sample size

- LICENSE
- README