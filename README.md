# Guide

## Downloading and preparing raw data
1. Download raw data as **fastq.gz** files from **SRA**. 
2. Run **FastQC** on the **fastq.gz** files to check the quality of the raw sequencing reads.
3. Make a **bowtie2 index** of the zebrafish genome (danRer10) if you do not have it already, or download a premade index from here.
4. Align **fastq.gz** reads to the indexed genome using `bowtie2` to obtain **SAM** files.
5. Filter **SAM** files using **samtools** to retain reads that align with a maximum of one mismatch as **BAM** files.
6. Sort and index **BAM** files using **samtools**.
7. Estimate fragment length using `macs2 predictd` from **MACS2**.
8. Make **bigWig** genome-wide coverage tracks from the *sorted* **BAM** files by extending reads to the predicted fragment length.
9. Call ChIP-peaks on *sorted* **BAM** files using `macs2 callpeak` from **MACS2**. Use *sorted* **BAM** files from *input DNA* as control. 

## Superenhancers

## Detailed steps

### Create folder structure
After entering the directory where you want to store all the data, create a variable called `DATAPATH` to store the directory path.

`export DATAPATH=$(pwd)`

Create directories to store raw and processed data. 

`mkdir raw_data fastqc_output genome_files danRer10_index aligned_data fragment_length_predict coverage_data peakcalled_data`

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
Run **FastQC** on each **fastq.gz** file to check quality of sequencing reads. For example,
`fastqc -o $DATAPATH/fastqc_output/ $DATAPATH/raw_data/SRR352249.fastq.gz`

**Pol II Ser5P** *Dome* data reveals a large fraction of duplicated reads. Is this expected? We also deduplicate reads later to check if peak detection and coverage around peaks is different.

### Creating genome index
Obtain the DNA reference sequence in `fasta` format from **NCBI**, unzip it, and use `bowtie2-build` to build a genome index.

`wget https://hgdownload.soe.ucsc.edu/goldenPath/danRer10/bigZips/danRer10.fa.gz -P $DATAPATH/genome_files/`

`gunzip $DATAPATH/genome_files/danRer10.fa.gz`

`bowtie2-build $DATAPATH/genome_files/danRer10.fa $DATAPATH/danRer10_index/danRer10`

### Aligning reads
Aligning reads to the indexed genome returns **SAM** files.

`bowtie2 -q -x $DATAPATH/danRer10_index/danRer10 -U raw_data/SRR352249.fastq.gz -S $DATAPATH/aligned_data/H3K27ac_Dome.sam`

> 11861357 reads; of these:
>   11861357 (100.00%) were unpaired; of these:
>     659229 (5.56%) aligned 0 times
>     6953084 (58.62%) aligned exactly 1 time
>     4249044 (35.82%) aligned >1 times
> 94.44% overall alignment rate

`bowtie2 -q -x $DATAPATH/danRer10_index/danRer10 -U raw_data/SRR352250.fastq.gz -S $DATAPATH/aligned_data/H3K27ac_80Epi.sam`

`bowtie2 -q -x $DATAPATH/danRer10_index/danRer10 -U raw_data/SRR711350.fastq.gz -S $DATAPATH/aligned_data/Ser5P_Dome.sam`

`bowtie2 -q -x $DATAPATH/danRer10_index/danRer10 -U raw_data/SSRR3932160.fastq.gz -S $DATAPATH/aligned_data/Input1_Dome.sam`

`bowtie2 -q -x $DATAPATH/danRer10_index/danRer10 -U raw_data/SSRR3932161.fastq.gz -S $DATAPATH/aligned_data/Input2_Dome.sam`

`bowtie2 -q -x $DATAPATH/danRer10_index/danRer10 -U raw_data/SRR585264.fastq.gz -S $DATAPATH/aligned_data/Input_80Epi.sam`

### Filtering, sorting and indexing aligned reads
Next, use **samtools** to filter aligned reads to retain those with a maximum of 1 mismatch, sort and index them. This will return **BAM** files and a **BAM.BAI** file. We use the following filename endings.

1. *\*_ut1m.bam*: Filtered reads with upto 1 mismatch.
2. *\*_ut1m\_sorted.bam*: Filtered and sorted.
3. *\*_ut1m\_sorted.bam.bai*: Indexed.

**Filtering**: `samtools view -bq 40 $DATAPATH/aligned_data/Ser5P_Dome.sam > $DATAPATH/aligned_data/Ser5P_Dome_ut1m.bam`
**Sorting**: `samtools sort $DATAPATH/aligned_data/Ser5P_Dome_ut1m.bam > $DATAPATH/aligned_data/Ser5P_Dome_ut1m_sorted.bam`
**Indexing**: `samtools index $DATAPATH/aligned_data/Ser5P_Dome_ut1m.bam`
 
For **Input** *Dome* data, we merge the BAM files before sorting. 
`samtools merge $DATAPATH/aligned_data/Input_Dome_ut1m.bam $DATAPATH/aligned_data/Input1_Dome_ut1m.bam $DATAPATH/aligned_data/Input2_Dome_ut1m.bam`

### Marking duplicates
**Pol II Ser5P** *Dome* data has a lot of duplicate reads. So we use `MarkDuplicates` from **picard** to mark potential duplicates. ![](/Users/rprizak/Documents/GitHub/ser5p-superenhancers-zebrafish/additional_files/Ser5P_Dome_duplication_levels.png)

`MarkDuplicates -I $DATAPATH/aligned_data/Ser5P_Dome_ut1m_sorted.bam -O $DATAPATH/aligned_data/Ser5P_Dome_ut1m_sorted_dm.bam -M $DATAPATH/aligned_data/Ser5P_Dome_marked_dup_metrics.txt`

### Predicting fragment length
To produce **bigWig** coverage tracks, we need an estimate of the fragment length from the ChIP-seq experiments. We use `macs2 predictd` from **MACS2** to do this. We choose effective genome size based on read length using [this table](https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html) from **deepTools**.

**H3K27ac** *Dome*
`macs2 filterdup -i $DATAPATH/aligned_data/H3K27ac_Dome_ut1m_sorted.bam --keep-dup=1 -o $DATAPATH/fragment_length_predict/H3K27ac_Dome_filterdup.bed`
`macs2 predictd -g 1195445591 --outdir $DATAPATH/fragment_length_predict/ -i $DATAPATH/fragment_length_predict/H3K27ac_Dome_filterdup.bed --rfile H3K27ac_Dome -m 5 50`

> tag size = 49, total tags in alignment file: 5772149, number of paired peaks: 50433, **predicted fragment length is 287** bps, alternative fragment length(s) may be 287 bps 

**H3K27ac** *80Epi*
`macs2 filterdup -i $DATAPATH/aligned_data/H3K27ac_80Epi_ut1m_sorted.bam --keep-dup=1 -o $DATAPATH/fragment_length_predict/H3K27ac_80Epi_filterdup.bed`
`macs2 predictd -g 1195445591 --outdir $DATAPATH/fragment_length_predict/ -i $DATAPATH/fragment_length_predict/H3K27ac_80Epi_filterdup.bed --rfile H3K27ac_80Epi -m 5 50`

> tag size = 49, total tags in alignment file: 5580260, number of paired peaks: 43941, **predicted fragment length is 297** bps, alternative fragment length(s) may be 297 bps 

**Ser5P** *Dome*
`macs2 filterdup -i $DATAPATH/aligned_data/Ser5P_Dome_ut1m_sorted.bam --keep-dup=1 -o $DATAPATH/fragment_length_predict/Ser5P_Dome_filterdup.bed`
`macs2 predictd -g 1251132686 --outdir $DATAPATH/fragment_length_predict/ -i $DATAPATH/fragment_length_predict/Ser5P_Dome_filterdup.bed --rfile Ser5P_Dome -m 5 50`

> tag size = 76, total tags in alignment file: 10681375, number of paired peaks: 23456, **predicted fragment length is 224** bps, alternative fragment length(s) may be 224 bps 

**Input** *Dome*
`macs2 filterdup -i $DATAPATH/aligned_data/Input_Dome_ut1m_sorted.bam --keep-dup=1 -o $DATAPATH/fragment_length_predict/Input_Dome_filterdup.bed`
`macs2 predictd -g 1195445591 --outdir $DATAPATH/fragment_length_predict/ -i $DATAPATH/fragment_length_predict/Input_Dome_filterdup.bed --rfile Input_Dome -m 3 20`

> tag size = 49, total tags in alignment file: 21866194, number of paired peaks: 9428, **predicted fragment length is 191** bps, alternative fragment length(s) may be 191 bps 

**Input** *80Epi*
`macs2 filterdup -i $DATAPATH/aligned_data/Input_80Epi_ut1m_sorted.bam --keep-dup=1 -o $DATAPATH/fragment_length_predict/Input_80Epi_filterdup.bed`
`macs2 predictd -g 1195445591 --outdir $DATAPATH/fragment_length_predict/ -i $DATAPATH/fragment_length_predict/Input_80Epi_filterdup.bed --rfile Input_80Epi -m 3 20`

> tag size = 35, total tags in alignment file: 8355884, number of paired peaks: 34209, **predicted fragment length is 181** bps, alternative fragment length(s) may be 181 bps 

### Making bigWig coverage files
Next we use `bamCoverage` from **deepTools** to extend reads and make a **bigWig** coverage file using the RPGC normalization method. This normalizes the reads using 1x coverage, making comparison with other tracks better. We use a bin size of 30 bps and extend reads to the fragment length predicted using `macs2 predictd`.

`bamCoverage -bs 30 -p max/2 --normalizeUsing RPGC --ignoreDuplicates --effectiveGenomeSize 1195445591 --extendReads 287 -b $DATAPATH/aligned_data/H3K27ac_Dome_ut1m_sorted.bam -o $DATAPATH/coverage_data/H3K27ac_Dome_ut1m.bigWig`
`bamCoverage -bs 30 -p max/2 --normalizeUsing RPGC --ignoreDuplicates --effectiveGenomeSize 1195445591 --extendReads 297 -b $DATAPATH/aligned_data/H3K27ac_80Epi_ut1m_sorted.bam -o $DATAPATH/coverage_data/H3K27ac_80Epi_ut1m.bigWig`
`bamCoverage -bs 30 -p max/2 --normalizeUsing RPGC --ignoreDuplicates --effectiveGenomeSize 1251132686 --extendReads 224 -b $DATAPATH/aligned_data/Ser5P_Dome_ut1m_sorted.bam -o $DATAPATH/coverage_data/Ser5P_Dome_ut1m.bigWig`
`bamCoverage -bs 30 -p max/2 --normalizeUsing RPGC --ignoreDuplicates --effectiveGenomeSize 1195445591 --extendReads 191 -b $DATAPATH/aligned_data/Input_Dome_ut1m_sorted.bam -o $DATAPATH/coverage_data/Input_Dome_ut1m.bigWig`
`bamCoverage -bs 30 -p max/2 --normalizeUsing RPGC --ignoreDuplicates --effectiveGenomeSize 1195445591 --extendReads 181 -b $DATAPATH/aligned_data/Input_80Epi_ut1m_sorted.bam -o $DATAPATH/coverage_data/Input_80EPi_ut1m.bigWig`

### Calling peaks
Next we call peaks using `macs2 callpeak`. 
`macs2 callpeak -g 1369631918 -q 0.05 --broad -n test --nomodel --extsize 200 --keep-dup all -t X_upto1mismatch_sorted.bam -c I_upto1mismatch_sorted.bam`