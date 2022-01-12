# YeastTranslatome
YeastTranslatome is a command line program consisting multiple subprocesses. 
The general scheme for running the program is to execute YeastTranslatome followed by a single argument corresponding to a subprocess to run.

```bash
YeastTranslatome -ProcessNameFromListBelow
```


Note that several processes require external data or programs. Processes must be run in the given order as they often depend on successful completion of previous steps. Run the following processes in order to generate all required files to perform the analysis conducted in this study:

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

The following processes will map and combine Riboseq data from several studies.
- MapRiboseqReads
- CombineRiboseqReads

The output of these processes is a pair of files mapping processed ribo-seq reads to both strands of the S. cerevisiae genome:

- `mapped_ribseqs_combo_all_r`
- `mapped_ribseqs_combo_all_f`

These files could be accessed from ***figshare_link***

As downloading and processing all required fastq files takes approximately 3 weeks on a cluster, the output files for the analysis conducted and required to generate plots shown in our manuscript are already included in the `input_files` directory.  


## Dependencies 
The last two processes (`MapRiboseqReads` and `CombineRiboseqReads`) require a collection of fastq files with ribosome profiling data. A full list of the fastq files analyzed in this study is included in: `input_files/riboseq_experiments_dataset.txt`

The file 1011Matrix.gvcf is required for population analyses. Download from: http://1002genomes.u-strasbg.fr/files/ 

The PRANK and MUSCLE programs are called for sequence alignments and must be available.