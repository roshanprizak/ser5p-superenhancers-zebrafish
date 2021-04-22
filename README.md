# Guide

## Downloading and preparing raw data
1. Download raw data as `fastq.gz` files from SRA. 
2. Run `fastqc` on the `fastq.gz` files to check the quality of the raw sequencing reads, and use `trimmomatic` to cut any adaptor sequences present.
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

`mkdir raw_data aligned_data coverage_data peakcalled_data danRer10_index`

### Download reads
We get **H3K27ac** ChIP-seq data (`SRR352249.fastq.gz`) at *Dome* and *80% Epiboly* stages from [GSE32483](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE32483), **Pol II Ser5P** ChIP-seq data at *Dome* stage from [GSE44269](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE44269), **Input** control ChIP-seq data at *Dome* from [GSE84602](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84602
), and **Input** control ChIP-seq data at *80% Epiboly* from [GSE41458](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41458).

`wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR352/SRR352249/SRR352249.fastq.gz -P $DATAPATH/raw_data/`
