# Trimviz
**Visualisation of trimming or soft-clipping of short sequence reads**

for help and dependencies:

`python ./trimviz.py -h`

<br>

### Quickstart:

**fastq-fastq comparison (trimming analysis):**

`python path/to/trimviz.py FQ -o tv_outdir -u untrimmed_R1.fastq.gz -t trimmed_R1.fastq.gz`

[Example report - FQ mode](http://MonashBioinformaticsPlatform.github.io/trimviz/example_reports/tvFQ_simple/trimvis_report.html)

<br>

**fastq-bam comparison (soft-clipping analysis):**

`python path/to/trimviz.py SC -o tv_outdir -u aligner_input_R1.fastq.gz -b aligned.bam -g path/to/reference_fasta.fa`

[Example report - SC mode](http://MonashBioinformaticsPlatform.github.io/trimviz/example_reports/tvSC/trimvis_report.html)

<br>

**fastq-fastq comparison (trimming analysis) with downstream alignment context in diff mode:**

`python path/to/trimviz.py FQ -o tv_outdir -u untrimmed_R1.fastq.gz -t trimmed_R1.fastq.gz -b aligned.bam -g path/to/reference_fasta.fa -d`

[Example report - FQ + bam mode](http://MonashBioinformaticsPlatform.github.io/trimviz/example_reports/tvFQ_withbam/trimvis_report.html)

<br>

**Problematic read-2 dataset: soft-clipping analysis with diff mode:**

`python path/to/trimviz.py SC -o tv_outdir -u aligner_input_R2.fastq.gz -b aligned.bam -g path/to/reference_fasta.fa -d -e `

[Example report - FQ + bam mode](http://MonashBioinformaticsPlatform.github.io/trimviz/example_reports/tvSC_R2/trimvis_report.html)

<br>

The above analysis suggests a problem with the original sequence data (41,250 of 50,000 reads were 3'-soft-clipped by the aligner). 
The input fastq had already been trimmed for quality and adapter sequences. However the last 4 bp still almost never aligned with the reference, 
suggesting that they originated elsewhere.

**Note on paired-end reads:**  if the fastq file is Read 2, set the `-e/--read2only` flag so that trimviz knows to extract the Read 2 alignment from the bam file.

### Command-line arguments:

Trimviz takes a random sample of untrimmed reads from a fastq file,
looks up the same reads in an trimmed fastq file and visualises the
trimmed reads with respect to surrounding base call quality values and
adapter sequence. In soft-clipping mode, Trimviz will instead
visualize the soft-clipping of reads by an aligner. To examine Read 2
from paired-end data, set the `-e/--read2only` flag.

```   
    Usage:
    ./trimviz.py FQ -o/-O output_dir -u untrimmed_R1.fq.gz -t trimmed_R1.fq.gz [ -b align.bam -g reference.fa ]
    ./trimviz.py SC -o/-O output_dir -u unaligned_R1.fq.gz -b align.bam -g reference.fa
    
    trimsviz.py FQ        Fastq-fastq comparison. Bam file and genome fasta file can be optionally given to view the mapping outcomes for trimmed reads.
    trimsviz.py SC        Treat soft clipping as the trimming of interest. Bam and genome fasta file are required, with only one fastq file.
    
    options:
    -o/--out_dir          Directory for output. If it already exists, an error will be generated.
    -O/--out_dir_fat      Directory for output (retain temporary files; choose this option to keep the sub-sampled fastq files)
    -u/--untrimmed        Untrimmed (original) fastq file name. In SC mode, this should be the fastq file that was directly input into the aligner.
    -t/--trimmed          Trimmed fastq file name
    -U/--untrimmed_r2     Read2 fastq file corresponding to -u file (not implemented yet)
    -T/--trimmed_r2       Read2 fastq file corresponding to -t file (not implemented yet)
    -b/--bam              Bam file (optional in FQ mode; required in SC mode). Only R1 alignments will be extracted (or only R2 if -e is set).
    -g/--genome_fasta     Fasta file of genome sequence (required if using .bam alignment)
    -c/--classes          [uncut,3pcut,removed,5pcut] Comma-separated trim-classes to visualise in individual read visualisation.
                          One/several of 'uncut','5pcut','3pcut','removed','generated_warning','indel'
    -a/--adapt:           comma-separated adapter sequences to highlight (multiple adapters not supported yet - defaults to first adapter in list)
    -A/--adaptfile        Text file containing adapter sequences (multiple adapters not supported yet - defaults to first adapter in list)
    -n/--sample_size      [50000] internal parameter: max reads to subsample in file (should be >> -m and -v, especially if only a small proportion are trimmed)
    -v/--nvis             [20] number of reads in each category (or in total if -R is set) to use for detailed individual plots
    -w/--heatmap_reads    [200] number of reads to plot in heatmaps
    -f/--agg_flank        [20] number of flanking nucleotides around trim point to plot in heatmaps 
    -r/--rid_file         File of read-ids to select, instead of using random sampling
    -s/--rseed            [1] random seed for sampling
    
    flags:
    -k/--skim             Speed up by skimming the reads from the tops of the fastq files (warning: these will be edge-of-flowcell reads)
    -R/--representative   Ignore read-classes and take a representative sample. (This often results in untrimmed reads dominating the visualization)
    -z/--gzipped          Assume fastq files are gzipped (default behaviour is to guess via .gz file extension)
    -e/--read2only        If set, only alignments with the 'read-2' flag will be extracted from the .bam file (default: only read-1 is extracted)
    -d/--diff             When displaying genomic alignment context from bam file, only display nucleotides that differ from the read sequence
    -h/--help:            Print help
```

<br>

### Dependencies:

**Command-line programs:**

Rscript <br>
seqtk <br>
zcat <br>
fgrep <br>

**python libraries:**

getopt, subprocess, random, re, sys, os, gzip, pipes, pysam

**R libraries:**

ggplot2, ape, reshape2, gridExtra

Tested with R3.6.0 and Python 2.7.6
