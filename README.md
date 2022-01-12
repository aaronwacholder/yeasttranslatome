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



### Dependencies 
The last two processes (`MapRiboseqReads` and `CombineRiboseqReads`) require a collection of fastq files with ribosome profiling data. A full list of the fastq files analyzed in this study is included in: `input_files/riboseq_experiments_dataset.txt`

The file 1011Matrix.gvcf is required for population analyses. Download from: http://1002genomes.u-strasbg.fr/files/ 

The PRANK and MUSCLE programs are called for sequence alignments and must be available.

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
