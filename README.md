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
