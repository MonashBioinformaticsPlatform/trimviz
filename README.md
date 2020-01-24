# trimviz
**trimviz: intuitive visualization and QC of trimming and soft-clipping of short DNA sequencing reads**

for help and dependencies:

`python ./trimviz.py -h`

### Quickstart:

**fastq-fastq comparison (trimming analysis):**

`python path/to/trimviz.py FQ -o tv_outdir -u untrimmed.fastq.gz -t trimmed.fastq.gz`

**fastq-bam comparison (soft-clipping analysis):**

`python path/to/trimviz.py SC -o tv_outdir -u aligner_input.fastq.gz -b aligned.bam -g path/to/reference_fasta.fa`
