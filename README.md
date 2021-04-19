# ser5p-superenhancers-zebrafish

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
> aksjdakjsjka 