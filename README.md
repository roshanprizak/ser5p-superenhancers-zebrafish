# Guide

## Downloading and preparing raw data
1. Download raw data as `fastq.gz` files from **SRA**. 
2. Run **FastQC** on the `fastq.gz` files to check the quality of the raw sequencing reads, and use `trimmomatic` to cut any adaptor sequences present.
3. Make a *bowtie2 index* of the zebrafish genome (danRer10) if you do not have it already, or download a premade index from here.
4. Align `fastq.gz` reads to the indexed genome using `bowtie2` to obtain `SAM` files.
5. Filter `SAM` files using `samtools` to retain reads that align with a maximum of one mismatch as `BAM` files.
6. Sort and index `BAM` files using `samtools`.
7. Estimate fragment length using `macs2 predictd`.
8. Make `bigWig` genome-wide coverage tracks from the *sorted* `BAM` files by extending reads to the predicted fragment length.
9. Call ChIP-peaks on *sorted* `BAM` files using `macs2 callpeak`. Use sorted `BAM` files from *input DNA* as control. 

## Superenhancers

## Detailed steps

### Create folder structure
After entering the directory where you want to store all the data, create a variable called `DATAPATH` to store the directory path.

`export DATAPATH=$(pwd)`

Create directories to store raw and processed data. 

`mkdir raw_data fastqc_output danRer10_index aligned_data coverage_data peakcalled_data`

### Download reads
We get **H3K27ac** ChIP-seq data at *Dome* (`SRR352249.fastq.gz`) and *80% Epiboly* (`SRR352250.fastq.gz`) stages from [GSE32483](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE32483), **Pol II Ser5P** ChIP-seq data at *Dome* (`SRR711350.fastq.gz`) stage from [GSE44269](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE44269), **Input** control ChIP-seq data at *Dome* (`SRR3932160.fastq.gz` and `SRR3932161.fastq.gz`) from [GSE84602](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84602
), and **Input** control ChIP-seq data at *80% Epiboly* (`SRR585264.fastq.gz`) from [GSE41458](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41458).

**H3K27ac** *Dome*: `wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR352/SRR352249/SRR352249.fastq.gz -P $DATAPATH/raw_data/`
**H3K27ac** *80% Epiboly*:  `wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR352/SRR352250/SRR352250.fastq.gz -P $DATAPATH/raw_data/`
**Pol II Ser5P** *Dome*:  `wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR711/SRR711350/SRR711350.fastq.gz -P $DATAPATH/raw_data/`
**Input** *Dome*: `wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR393/000/SRR3932160/SRR3932160.fastq.gz -P $DATAPATH/raw_data/`
`wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR393/001/SRR3932161/SRR3932161.fastq.gz -P $DATAPATH/raw_data/`
**Input** *80% Epiboly*: `wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR585/SRR585264/SRR585264.fastq.gz -P $DATAPATH/raw_data/`

### Running **FastQC**
Run **FastQC** on each `fastq.gz` file to check quality of sequencing reads. For example,
`fastqc -o $DATAPATH/fastqc_output/ $DATAPATH/raw_data/SRR352249.fastq.gz`


 