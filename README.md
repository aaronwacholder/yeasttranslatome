# YeastTranslatome
## Analyses
YeastTranslatome is a command line program consisting multiple subprocesses. 

In order to run it you first need to compile from source file (`YeastTranslatome.cpp`). Note that several processes require external data or programs to run successfully. See dependencies below.

```bash
g++ -o YeastTranslatome YeastTranslatome.cpp 
```

The general scheme for running the program is to execute YeastTranslatome followed by a single argument corresponding to a subprocess to run.

```bash
YeastTranslatome -ProcessNameFromListBelow
```

The following processes will map and combine Riboseq data from several studies and is required to be completed before other processes.

- MapRiboseqReads
- CombineRiboseqReads

The output of these processes is a pair of files mapping processed ribo-seq reads to both strands of the S. cerevisiae genome:

- `mapped_ribseqs_combo_all_r` // there are 24 subprocesses and no dependencies are known for those. 
- `mapped_ribseqs_combo_all_f`


As downloading and processing all required fastq files takes approximately 3 weeks on a cluster, These processed files could be accessed from ***figshare_link***. 

Other output files required to generate plots shown in our manuscript are already included in the `input_files` directory.


Processes must be run in the given order as they often depend on successful completion of previous steps. Run the following processes in order to generate all required files to perform the analysis conducted in this study:

- GetAllORFs
- IdentifyTranslatedORFs 
- GetSubsMatrix
- SyntenyAlign
- SyntenyMatch
- SyntenyCheck
- SensuStrictoBlast
- BlastAlign
- BlastMatch
- ResolveAlignments
- ValidateAlignments
- MultAlignOrthologs
- PhyloORFs
- PhyloBlast
- PhyloAnalyze
- SelfBlast
- CodingScores
- MakeCleanAlignments
- MakeJoinedCleanAlignments
- IndividualDNDS
- StrainSetup
- StrainAnalyze
- NucDiverse
- CollectiveDNDS
- AntiDNDS

##Input Files
The following files are included in the input_files directory and are necessary for successfully running the program:

Genomes: 
GCF_001298625.1_SEUB3.0_genomic.fna
Sbay.ultrascaf
GCF_000292725.1_SacArb1.0_genomic.fna
Skud.ultrascaf
GCA_900290405.1_SacJureiUoM1_genomic.fna
S288C_reference_sequence_R64-2-1_20150113.fsa
Smik.ultrascaf
Spar.ultrascaf

Saccharomyces Genome Database (https://www.yeastgenome.org/) genome annotations:
saccharomyces_cerevisiae.gff
other_features_genomic.fasta
rna_coding.fasta
orf_coding_all.fasta

Tif-seq data from  Pelechano et al. 2013 (PMID 23615609), downloaded from SGD:
tsedAnno_V2.txt

Information prepared for this study:
budding_yeast_genomes_filenames.txt
riboseq_experiments_dataset.txt

Figures:
flowchart_evolve14.eps
Figure1_schema4_nolabel.eps
Schema_RiboInt6.eps

### Dependencies 
The last two processes (`MapRiboseqReads` and `CombineRiboseqReads`) require a collection of fastq files with ribosome profiling data. A full list of the fastq files analyzed in this study is included in: `input_files/riboseq_experiments_dataset.txt`

The file 1011Matrix.gvcf is required for population analyses. Download from: http://1002genomes.u-strasbg.fr/files/ 

The budding yeast genomes collected by Shen et al. 2018 is also required. Download 0_332yeast_genomes.zip from the figshare repository: https://figshare.com/articles/dataset/Tempo_and_mode_of_genome_evolution_in_the_budding_yeast_subphylum/5854692

Then place all the genome assemblies (the .fas files) in the folder budding_yeast_genomes. The file budding_yeast_genomes_filenames.txt should list all genomes to process.

The BLAST, PRANK and MUSCLE programs are called for sequence alignments and must be available.

## Figure generation
`translatome.R` script uses files given in `input_files` and others generated when above processes run.

Following libraries are required to run `translateome.R`
```R
require("ggpubr")
require("grImport") 
require("grid")
require("gridExtra")
require("ggplot2")
require("cowplot")
require("scales")
require("heatmap3")
require("ggcorrplot")
require("tidyr")
require("plotrix")
library("tidyverse")
library("glue")
library("igraph")
library("ggsignif")
library("cowplot")
library("readxl")
library("coin")
```

## Function Input and Outputs

-GetAllORFs
Ouput: orf_comp, orf_overlap, shared_transcripts, orf_neighbors
-IdentifyTranslatedORFs 
Input: orfs_comp
Output: many files with names starting with "riboseq_orfs" corresponding to different sets of studies 
-GetSubsMatrix
Output: sub_rates_nongenic
-SyntenyAlign
Output: orfs_aa_Suva.fas, orfs_aa_Spar.fas, orfs_aa_Smik.fas, orfs_aa_Skud.fas, orfs_aa_Sjur.fas, orfs_aa_Seub.fas, orfs_aa_Scer.fas, orfs_aa_Sarb.fas
Output: blocks7_0, blocks6_0, blocks5_0, blocks4_0, blocks3_0, blocks2_0, blocks1_0
Output: blocks_Scer_Suva_0, blocks_Scer_Spar_0, blocks_Scer_Smik_0, blocks_Scer_Skud_0, blocks_Scer_Sjur_0, blocks_Scer_Seub_0, blocks_Scer_Sarb_0
Output: aligned_blocks1_0, aligned_blocks2_0, aligned_blocks3_0, aligned_blocks4_0, aligned_blocks5_0, aligned_blocks6_0, aligned_blocks7_0
-SyntenyMatch
Input: sub_rates_nongenic, aligned_blocks1_0, aligned_blocks2_0, aligned_blocks3_0, aligned_blocks4_0, aligned_blocks5_0, aligned_blocks6_0, aligned_blocks7_0
Output: orf_bounds_synteny, orfs_Scer_synteny_matched
-SyntenyCheck
Input: orfs_Scer_synteny_matched
Output: orfs_scer_synteny_matched_checked
-SensuStrictoBlast
Output: blasts_genome_nuc_Sarb_3.out, blasts_genome_nuc_Scer_3.out, blasts_genome_nuc_Seub_3.out, blasts_genome_nuc_Sjur_3.out, blasts_genome_nuc_Skud_3.out, blasts_genome_nuc_Smik_3.out, blasts_genome_nuc_Spar_3.out, blasts_genome_nuc_Suva_3.out
-BlastAlign
Input: blasts_genome_nuc_Sarb_3.out, blasts_genome_nuc_Scer_3.out, blasts_genome_nuc_Seub_3.out, blasts_genome_nuc_Sjur_3.out, blasts_genome_nuc_Skud_3.out, blasts_genome_nuc_Smik_3.out, blasts_genome_nuc_Spar_3.out, blasts_genome_nuc_Suva_3.out
Output: blocks_blast_1, blocks_blast_2, blocks_blast_3, blocks_blast_4, blocks_blast_5, blocks_blast_6, blocks_blast_7
Output: aligned_blocks_blast_1, aligned_blocks_blast_2, aligned_blocks_blast_3, aligned_blocks_blast_4, aligned_blocks_blast_5, aligned_blocks_blast_6, aligned_blocks_blast_7
-BlastMatch
Input: aligned_blocks_blast_1, aligned_blocks_blast_2, aligned_blocks_blast_3, aligned_blocks_blast_4, aligned_blocks_blast_5, aligned_blocks_blast_6, aligned_blocks_blast_7
Input: orfs_scer_synteny_matched_checked
Output: orf_bounds_blast, orfs_scer_scer_blast_matched
-ResolveAlignments
Input: orfs_scer_blast_matched, orf_bounds_blast
Output: blast_matches_to_validate.fasta, blast_matches_to_validate.out
-ValidateAlignments
Input: orfs_scer_blast_matched, blast_matches_to_validate.out
Output: orfs_scer_blast_matched_validated
-MultAlignOrthologs
Input: orfs_scer_blast_matched_validated
Output: orfs_mult_align
-PhyloORFs
Output: ORF data files in directory ORFS_332
-PhyloBlast -TBLASTN:
Output: BLAST result files in directory tblastn_align
-PhyloBlast -TBLASTN_SCRAMBLED:
Output: BLAST result files in directory tblastn_align_scrambled_and_real
-PhyloBlast -ORF_BLAST
Output: BLAST result files in directory orf_blast_align
-PhyloAnalyze -TBLASTN
Input: orfs_comp, blast result files
Output: phylo_pres
-PhyloAnalyze -TBLASTN_SCRAMBLED
Input: orfs_comp, blast result files
Output: phylo_pres_tblastn_scrambled
-SelfBlast
Input: orf_genomic_all.fasta
Output: gene_sim
-CodingScores
Input: orfs_comp
Output: coding_scores_
-MakeCleanAlignments
Input: orfs_scer_blast_matched_validated
Output: clean_alignments
-MakeJoinedCleanAlignments
Input: orfs_comp
Output: joined_tree.nwk
-IndividualDNDS
Input: orfs_comp, sub_rates_nongenic
Output: tree_dnds
-StrainSetup
Input: orfs_comp
Output: align_positions
-StrainAnalyze
Input: sub_rates_nongenic, orfs_comp, align_positions
Output: anti_pnps_pop_1, anti_pnps_pop_2, anti_pnps_pop_3
-NucDiverse
Input: orfs_scer_blast_matched_validated
Output: nuc_diverse
-CollectiveDNDS
Input: orfs_scer_blast_matched_validated
Output: spar_dnds_all_
-AntiDNDS
Input: orfs_scer_blast_matched_validated
Output: spar_antisyn_all_1, spar_antisyn_all_2, spar_antisyn_all_3, spar_antinonsyn_all_1, spar_antinonsyn_all_2, spar_antinonsyn_all_3
