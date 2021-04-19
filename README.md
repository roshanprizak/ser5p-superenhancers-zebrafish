# ser5p-superenhancers-zebrafish

## Downloading and preparing raw data
1. Download raw data as `fastq.gz` files from SRA. 
3. Run `fastqc` on the `fastq.gz` files to check the quality of the raw sequencing reads, and use `trimmomatic` to cut any adaptor sequences present.
4. Make a bowtie2 index of the zebrafish genome (danRer10) if you do not have it already, or download a premade index from here.