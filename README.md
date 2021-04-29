# Aromatic dimer dehydrogenases from Novosphingobium aromaticivorans reduce aromatic diketones

Alexandra M. Linz, Yanjun Ma, Jose M. Perez, Kevin Myers, Wayne S. Kontur, Daniel R. Noguera, Timothy J. Donohue

DOE Great Lakes Bioenergy Research Center, Madison, WI, USA
Wisconsin Energy Institute, University of Wisconsin-Madison, Madison, WI, USA
Department of Civil and Environmental Engineering, University of Wisconsin-Madison, Madison, WI, USA.
Department of Bacteriology, University of Wisconsin-Madison, Madison, WI, USA

Direct correspondence to: Timothy J. Donohue, email: tdonohue@bact.wisc.edu, phone: 6082624663

## About this page

In this study, we describe the ability of Novosphingobium aromaticivorans, a sphingomonad with a versatile aromatic metabolism, to process aromatic diketones, a potential substrate resulting from lignin depolymerization. We used RNA-Seq to generate a hypothesis about which genes encode the enzymes responsible for G-diketone degradation, then confirmed this hypothesis via -in vitro- enzyme assays.

Here, you can find data, scripts, and manuscript materials related to this study. If you have any questions, please do not hesitate to contact us.

## Repo Structure

- Enzyme_kinetics (Data and calculations for Michaelis-Menten kinetics)
  - compiled_kinetics_2021-01-12.xlsx (Measured slopes of enzyme oxidation/reduction and results)
  - enzyme_assays_2021-01-12.xlsx
  - ligLND_enzymes_2020-07-09.csv (HPLC measurements of diketone degradation)
  - mm-kinetics_curve_fitting.R (Script for fitting Michaelis-Menten co-efficients based on slopes)
- Growth_curves (Growth of -N. aromaticivorans- on G-diketone
  - Diketone_degradation_HPLC_2020-07-01.csv (Detection of substrate and intermediates during growth)
  - HPLC_diketone_growth_curves.xlsx
  - ligLNDO-deletion_growth_curves_cleaned.csv (Growth of deletion mutants on glucose and G-diketone)
- Manuscript (Materials submitted for publication)
- Plots (Intermediate plots generated for figures)
- RNA-Seq
  - Code (Scripts used to process RNA-Seq data)
    - tidy_data.R (Import RNA-Seq data into R datasets)
  - Data (.csv to open in Excel, .rds to import into R)
    - Differential_expression_testing_results (log fold change and FDR-corrected p-values for gene expression comparisons between substrates)
    - GDK_downreg_genes.csv/GDK_upreg_genes.csv (list of genes significantly up or downregulated on the G-diketone compared to glucose, but not other substrates)
    - Gene_information (locus tags, gene location, and annotation information for -N. aromaticivorans-)
    - gene_names (List of gene names and locus tags of published aromatic genes)
    - RNA-Seq_HPLC_results (HPLC data of aromatics detected in supernatant of cultures harvested for RNA)
    - RPKM_normalized_read_counts (Read counts normalized by sample size)
- figure_list.R
  - Script for taking all these datasets and generating manuscript figures
- LICENSE
- README