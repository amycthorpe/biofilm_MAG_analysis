# Biofilm MAG analysis
 
R scripts to analyse bacterial metagenome assembled genomes (MAGs) from river biofilms.  

### Initial processing and analysis 
1. Raw metagenomic reads are first processed according to [this workflow](https://github.com/amycthorpe/metag_analysis_EA) to assemble MAGs
2. Assembled MAGs are then analysed according to [this workflow](https://github.com/amycthorpe/EA_metag_post_analysis) to identify metabolic and functional traits

- Metagenomic reads and assembled MAGs have been deposited in the European Nucleotide Archive (ENA) at EMBL-EBI under accession number PRJEB85861

### Downstream analysis and visualisation
Outputs (available in the [Zenodo archive](https://doi.org/10.5281/zenodo.14762144)) can then be analysed with the provided R scripts to:
   * Assess the biogeography and taxonomic composition of biofilm bacterial communities
   * Determine the metabolic and functional potential of biofilm bacterial communities
   * Identify environmental drivers shaping biofilm bacterial communities with variance partitioning

### Citation
https://www.biorxiv.org/content/10.1101/2025.03.05.641783v1
