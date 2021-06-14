# Guide

## Tools needed

1. bowtie2 [build OK, alignment OK] [2.3.5.1 -> 2.4.2]
2. FastQC [OK]
3. samtools [filtering slightly different. can 1.9 be installed? yes but filtering slightly different again!] [1.9 -> 1.12]
4. MACS2 [filterdup, predictd lsightly different but OK, same length predicted]
5. deepTools [bamCoverage OK]
6. bedtools
7. pyGenomeTracks (requires Python)
8. R
	* ggplot2
	* ggExtra

You can either install all of these in a conda environment, or if you already have them installed, use those versions. We briefly outline the steps needed to put them all in a conda environment.

```bash
conda create -n ser5p-se
conda activate ser5p-se
conda install -c bioconda samtools bowtie2 fastqc macs2 deeptools bedtools pygenometracks
conda install -c bioconda samtools=1.12 bowtie2=2.4.2 fastqc=0.11.9 macs2 deeptools bedtools pygenometracks
conda install 
```

## Downloading and preparing raw data
1. Download raw data as **fastq.gz** files from **SRA**. 
2. Run **FastQC** on the **fastq.gz** files to check the quality of the raw sequencing reads.
3. Make a **bowtie2 index** of the zebrafish genome (danRer10) if you do not have it already, or download a premade index from here.
4. Align **fastq.gz** reads to the indexed genome using `bowtie2` to obtain **SAM** files.
5. Filter **SAM** files using **samtools** to retain reads that align with a maximum of one mismatch as **BAM** files.
6. Optional duplicate removal. Sort and index **BAM** files using **samtools**.
7. Estimate fragment length using `macs2 predictd` from **MACS2**.
8. Make **bigWig** genome-wide coverage tracks from the *sorted* **BAM** files by extending reads to the predicted fragment length.
9. Call ChIP-peaks on *sorted* **BAM** files using `macs2 callpeak` from **MACS2**. Use *sorted* **BAM** files from *input DNA* as control. 

## Detailed steps

### Create folder structure
After entering the directory where you want to store all the data, create a variable called `DATAPATH` to store the directory path.

```bash
export DATAPATH=$(pwd)
```

Create directories to store raw and processed data. 

```bash
mkdir raw_data fastqc_output genome_files danRer10_index aligned_data fragment_length_predict coverage_data peakcalled_data superenhancers coverage_tracks
```
	
### Download reads

| Target | Stage | GSE Accession | Filename(s) |
| ----------- | ----------- | ----------- | ----------- |
| H3K27ac | Dome | [GSE32483](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE32483) | SRR352249.fastq.gz | 
| H3K27ac | 80% Epiboly | [GSE32483](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE32483) | SRR352250.fastq.gz |
| Pol II Ser5P | Dome | [GSE44269](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE44269) | SRR711350.fastq.gz |
| Input | Dome | [GSE84602](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84602) | SRR3932160.fastq.gz, SRR3932161.fastq.gz |
| Input | 80% Epiboly | [GSE41458](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41458) | SRR585264.fastq.gz |


```bash
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR352/SRR352249/SRR352249.fastq.gz -P $DATAPATH/raw_data/
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR352/SRR352250/SRR352250.fastq.gz -P $DATAPATH/raw_data/
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR711/SRR711350/SRR711350.fastq.gz -P $DATAPATH/raw_data/
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR393/000/SRR3932160/SRR3932160.fastq.gz -P $DATAPATH/raw_data/
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR393/001/SRR3932161/SRR3932161.fastq.gz -P $DATAPATH/raw_data/
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR585/SRR585264/SRR585264.fastq.gz -P $DATAPATH/raw_data/
```
### Running **FastQC**
Run **FastQC** on each **fastq.gz** file to check quality of sequencing reads. For example,

```bash
fastqc -o $DATAPATH/fastqc_output/ $DATAPATH/raw_data/SRR352249.fastq.gz
```
**Pol II Ser5P** *Dome* data reveals a large fraction of duplicated reads. Is this expected? We also deduplicate reads later to check if peak detection and coverage around peaks is different.

### Creating genome index
Obtain the DNA reference sequence in `fasta` format from **NCBI**, unzip it, and use `bowtie2-build` to build a genome index.

```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/danRer10/bigZips/danRer10.fa.gz -P $DATAPATH/genome_files/
gunzip $DATAPATH/genome_files/danRer10.fa.gz
bowtie2-build $DATAPATH/genome_files/danRer10.fa $DATAPATH/danRer10_index/danRer10
```
### Aligning reads
Aligning reads to the indexed genome returns **SAM** files.

```bash
bowtie2 -q -x $DATAPATH/danRer10_index/danRer10 -U raw_data/SRR352249.fastq.gz -S $DATAPATH/aligned_data/H3K27ac_Dome.sam
```
```bash
> 11861357 reads; of these:
>   11861357 (100.00%) were unpaired; of these:
>     659229 (5.56%) aligned 0 times
>     6953084 (58.62%) aligned exactly 1 time
>     4249044 (35.82%) aligned >1 times
> 94.44% overall alignment rate
```

```bash
bowtie2 -q -x $DATAPATH/danRer10_index/danRer10 -U raw_data/SRR352250.fastq.gz -S $DATAPATH/aligned_data/H3K27ac_80Epi.sam
bowtie2 -q -x $DATAPATH/danRer10_index/danRer10 -U raw_data/SRR711350.fastq.gz -S $DATAPATH/aligned_data/Ser5P_Dome.sam
bowtie2 -q -x $DATAPATH/danRer10_index/danRer10 -U raw_data/SSRR3932160.fastq.gz -S $DATAPATH/aligned_data/Input1_Dome.sam
bowtie2 -q -x $DATAPATH/danRer10_index/danRer10 -U raw_data/SSRR3932161.fastq.gz -S $DATAPATH/aligned_data/Input2_Dome.sam
bowtie2 -q -x $DATAPATH/danRer10_index/danRer10 -U raw_data/SRR585264.fastq.gz -S $DATAPATH/aligned_data/Input_80Epi.sam
```

### Filtering, sorting and indexing aligned reads
Next, use **samtools** to filter aligned reads to retain those with a maximum of 1 mismatch, sort and index them. This will return **BAM** files and a **BAM.BAI** file. We use the following filename endings.

1. *\*_ut1m.bam*: Filtered reads with upto 1 mismatch.
2. *\*_ut1m\_sorted.bam*: Filtered and sorted.
3. *\*_ut1m\_sorted.bam.bai*: Indexed.

**Filtering**: 

```bash
samtools view -bq 40 $DATAPATH/aligned_data/Ser5P_Dome.sam > $DATAPATH/aligned_data/Ser5P_Dome_ut1m.bam
```
**Sorting**: 

```bash
samtools sort $DATAPATH/aligned_data/Ser5P_Dome_ut1m.bam > $DATAPATH/aligned_data/Ser5P_Dome_ut1m_sorted.bam
```
**Indexing**: 

```bash
samtools index $DATAPATH/aligned_data/Ser5P_Dome_ut1m_sorted.bam
```
For **Input** *Dome* data, we merge the BAM files before sorting. 

```bash
samtools merge $DATAPATH/aligned_data/Input_Dome_ut1m.bam $DATAPATH/aligned_data/Input1_Dome_ut1m.bam $DATAPATH/aligned_data/Input2_Dome_ut1m.bam
```
### Removing duplicates
**Pol II Ser5P** *Dome* data has a lot of duplicate reads. So we use `markdup` from **samtools** to mark and remove potential duplicates. ![](/Users/rprizak/Documents/GitHub/ser5p-superenhancers-zebrafish/additional_files/Ser5P_Dome_duplication_levels.png)

After aligning and filtering, we remove duplicates and index.

```bash
samtools sort -n -o $DATAPATH/aligned_data/Ser5P_Dome_ut1m_ns.bam $DATAPATH/aligned_data/Ser5P_Dome_ut1m.bam
samtools fixmate -m $DATAPATH/aligned_data/Ser5P_Dome_ut1m_ns.bam $DATAPATH/aligned_data/Ser5P_Dome_ut1m_fm.bam
samtools sort -o $DATAPATH/aligned_data/Ser5P_Dome_ut1m_ps.bam $DATAPATH/aligned_data/Ser5P_Dome_ut1m_fm.bam
samtools markdup -r -s $DATAPATH/aligned_data/Ser5P_Dome_ut1m_ps.bam $DATAPATH/aligned_data/Ser5P_Dome_ut1m_dr.bam
samtools index $DATAPATH/aligned_data/Ser5P_Dome_ut1m_dr.bam
```

### Predicting fragment length
To produce **bigWig** coverage tracks, we need an estimate of the fragment length from the ChIP-seq experiments. We use `macs2 predictd` from **MACS2** to do this. We choose effective genome size based on read length using [this table](https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html) from **deepTools**.

**H3K27ac** *Dome*

```bash
macs2 filterdup -i $DATAPATH/aligned_data/H3K27ac_Dome_ut1m_sorted.bam --keep-dup=1 -o $DATAPATH/fragment_length_predict/H3K27ac_Dome_filterdup.bed
macs2 predictd -g 1195445591 --outdir $DATAPATH/fragment_length_predict/ -i $DATAPATH/fragment_length_predict/H3K27ac_Dome_filterdup.bed --rfile H3K27ac_Dome -m 5 50
```
> tag size = 49, total tags in alignment file: 5772149, number of paired peaks: 50433, **predicted fragment length is 287** bps, alternative fragment length(s) may be 287 bps 

**H3K27ac** *80Epi*

```bash
macs2 filterdup -i $DATAPATH/aligned_data/H3K27ac_80Epi_ut1m_sorted.bam --keep-dup=1 -o $DATAPATH/fragment_length_predict/H3K27ac_80Epi_filterdup.bed
macs2 predictd -g 1195445591 --outdir $DATAPATH/fragment_length_predict/ -i $DATAPATH/fragment_length_predict/H3K27ac_80Epi_filterdup.bed --rfile H3K27ac_80Epi -m 5 50
```
> tag size = 49, total tags in alignment file: 5580260, number of paired peaks: 43941, **predicted fragment length is 297** bps, alternative fragment length(s) may be 297 bps 

**Ser5P** *Dome*

```bash
macs2 filterdup -i $DATAPATH/aligned_data/Ser5P_Dome_ut1m_sorted.bam --keep-dup=1 -o $DATAPATH/fragment_length_predict/Ser5P_Dome_filterdup.bed
macs2 predictd -g 1251132686 --outdir $DATAPATH/fragment_length_predict/ -i $DATAPATH/fragment_length_predict/Ser5P_Dome_filterdup.bed --rfile Ser5P_Dome -m 5 50
```
> tag size = 76, total tags in alignment file: 10681375, number of paired peaks: 23456, **predicted fragment length is 224** bps, alternative fragment length(s) may be 224 bps 
 
**Input** *Dome*

```bash
macs2 filterdup -i $DATAPATH/aligned_data/Input_Dome_ut1m_sorted.bam --keep-dup=1 -o $DATAPATH/fragment_length_predict/Input_Dome_filterdup.bed
macs2 predictd -g 1195445591 --outdir $DATAPATH/fragment_length_predict/ -i $DATAPATH/fragment_length_predict/Input_Dome_filterdup.bed --rfile Input_Dome -m 3 20
```
> tag size = 49, total tags in alignment file: 21866194, number of paired peaks: 9428, **predicted fragment length is 191** bps, alternative fragment length(s) may be 191 bps 

**Input** *80Epi*

```bash
macs2 filterdup -i $DATAPATH/aligned_data/Input_80Epi_ut1m_sorted.bam --keep-dup=1 -o $DATAPATH/fragment_length_predict/Input_80Epi_filterdup.bed
macs2 predictd -g 1195445591 --outdir $DATAPATH/fragment_length_predict/ -i $DATAPATH/fragment_length_predict/Input_80Epi_filterdup.bed --rfile Input_80Epi -m 3 20
```

> tag size = 35, total tags in alignment file: 8355884, number of paired peaks: 34209, **predicted fragment length is 181** bps, alternative fragment length(s) may be 181 bps 

### Making bigWig coverage files
Next we use `bamCoverage` from **deepTools** to extend reads and make a **bigWig (.bw)** coverage file using the RPGC normalization method. This normalizes the reads using 1x coverage, making comparison with other tracks better. We use a bin size of 30 bps and extend reads to the fragment length predicted using `macs2 predictd`.

```bash
bamCoverage -bs 30 -p max/2 --normalizeUsing RPGC --effectiveGenomeSize 1195445591 --extendReads 287 -b $DATAPATH/aligned_data/H3K27ac_Dome_ut1m_sorted.bam -o $DATAPATH/coverage_data/H3K27ac_Dome_ut1m.bw
bamCoverage -bs 30 -p max/2 --normalizeUsing RPGC --effectiveGenomeSize 1195445591 --extendReads 297 -b $DATAPATH/aligned_data/H3K27ac_80Epi_ut1m_sorted.bam -o $DATAPATH/coverage_data/H3K27ac_80Epi_ut1m.bw
bamCoverage -bs 30 -p max/2 --normalizeUsing RPGC --effectiveGenomeSize 1251132686 --extendReads 224 -b $DATAPATH/aligned_data/Ser5P_Dome_ut1m_sorted.bam -o $DATAPATH/coverage_data/Ser5P_Dome_ut1m.bw
bamCoverage -bs 30 -p max/2 --normalizeUsing RPGC --effectiveGenomeSize 1195445591 --extendReads 191 -b $DATAPATH/aligned_data/Input_Dome_ut1m_sorted.bam -o $DATAPATH/coverage_data/Input_Dome_ut1m.bw
bamCoverage -bs 30 -p max/2 --normalizeUsing RPGC --effectiveGenomeSize 1195445591 --extendReads 181 -b $DATAPATH/aligned_data/Input_80Epi_ut1m_sorted.bam -o $DATAPATH/coverage_data/Input_80EPi_ut1m.bw
```
We also make **bigWig** files of **Pol II Ser5P** duplicate-removed-reads at *Dome* stage.

```bash
bamCoverage -bs 30 -p max/2 --normalizeUsing RPGC --effectiveGenomeSize 1251132686 --extendReads 224 -b $DATAPATH/aligned_data/Ser5P_Dome_ut1m_dr.bam -o $DATAPATH/coverage_data/Ser5P_Dome_ut1m_dr.bw
```
### Calling peaks
Next we call peaks using `macs2 callpeak`. 

**Pol II Ser5P** peaks:

```bash
macs2 callpeak --outdir $DATAPATH/peakcalled_data/ --nomodel --keep-dup all -q 0.0001 -n Ser5P_Dome_alldup --extsize 224 -g 1251132686 -t $DATAPATH/aligned_data/Ser5P_Dome_ut1m_sorted.bam -c $DATAPATH/aligned_data/Input_Dome_ut1m_sorted.bam
macs2 callpeak --outdir $DATAPATH/peakcalled_data/ --nomodel --keep-dup auto -q 0.0001 -n Ser5P_Dome_autodup --extsize 224 -g 1251132686 -t $DATAPATH/aligned_data/Ser5P_Dome_ut1m_sorted_dm.bam -c $DATAPATH/aligned_data/Input_Dome_ut1m_sorted.bam
macs2 callpeak --outdir $DATAPATH/peakcalled_data/ --nomodel --keep-dup 1 -q 0.0001 -n Ser5P_Dome_1dup --extsize 224 -g 1251132686 -t $DATAPATH/aligned_data/Ser5P_Dome_ut1m_sorted_dm.bam -c $DATAPATH/aligned_data/Input_Dome_ut1m_sorted.bam
```
When we later plot **Pol II Ser5P** coverage tracks and the peaks detected by **MACS2**, we will see that retaining only 1 duplicate per position results in the best peak detection.

**H3K27ac** peaks:

```bash
macs2 callpeak --outdir $DATAPATH/peakcalled_data/ --nomodel --keep-dup auto -q 0.05 --broad -n H3K27ac_Dome --extsize 287 -g 1195445591 -t $DATAPATH/aligned_data/H3K27ac_Dome_ut1m_sorted.bam -c $DATAPATH/aligned_data/Input_Dome_ut1m_sorted.bam
macs2 callpeak --outdir $DATAPATH/peakcalled_data/ --nomodel --keep-dup auto -q 0.05 --broad -n H3K27ac_80Epi --extsize 297 -g 1195445591 -t $DATAPATH/aligned_data/H3K27ac_80Epi_ut1m_sorted.bam -c $DATAPATH/aligned_data/Input_Dome_ut1m_sorted.bam
```

### Identifying superenhancers
To identify superenhancers, we use the [ROSE program](https://github.com/stjude/ROSE) from the Young lab. We cloned the repository on 23 April 2020, modified the code slightly (`bin/ROSE_main.py` and `bin/ROSE_geneMapper.py`), and also pass a **danRer10** annotation file: `annotation/danRer10_refseq.ucsc` as a custom genome. This file has been downloaded from the [UCSC Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables) using *UCSC RefSeq (refGene)* table. We provide this version as part of the package. 

```bash
PATHTO=/path/to/ROSE
PYTHONPATH=$PATHTO/lib
export PYTHONPATH
export PATH=$PATH:$PATHTO/bin
```

```bash
cp $DATAPATH/peakcalled_data/H3K27ac_Dome_peaks.broadPeak $DATAPATH/peakcalled_data/H3K27ac_Dome_peaks.bed
cp $DATAPATH/peakcalled_data/H3K27ac_80Epi_peaks.broadPeak $DATAPATH/peakcalled_data/H3K27ac_80Epi_peaks.bed
cd $DATAPATH/ROSE
python bin/ROSE_main.py --custom annotation/danRer10_refseq.ucsc -i $DATAPATH/peakcalled_data/H3K27ac_Dome_peaks.bed -r $DATAPATH/aligned_data/H3K27ac_Dome_ut1m_sorted.bam -o $DATAPATH/superenhancers/ -s 12500 -t 2500 -c $DATAPATH/aligned_data/Input_Dome_ut1m_sorted.bam
python bin/ROSE_main.py --custom annotation/danRer10_refseq.ucsc -i $DATAPATH/peakcalled_data/H3K27ac_80Epi_peaks.bed -r $DATAPATH/aligned_data/H3K27ac_80Epi_ut1m_sorted.bam -o $DATAPATH/superenhancers/ -s 12500 -t 2500 -c $DATAPATH/aligned_data/Input_80Epi_ut1m_sorted.bam
cd ..
sed '1,6d' $DATAPATH/superenhancers/H3K27ac_Dome_peaks_AllStitched.table.txt | awk -F'\t' -v OFS="\t" 'int($10)==0' | awk -F'\t' -v OFS="\t" '{print $2,$3,$4,$1,$7}' > $DATAPATH/superenhancers/Enhancers_Dome.bed
sed '1,6d' $DATAPATH/superenhancers/H3K27ac_Dome_peaks_AllStitched.table.txt | awk -F'\t' -v OFS="\t" 'int($10)==1' | awk -F'\t' -v OFS="\t" '{print $2,$3,$4,$1,$7}' > $DATAPATH/superenhancers/Superenhancers_Dome.bed
sed '1,6d' $DATAPATH/superenhancers/H3K27ac_80Epi_peaks_AllStitched.table.txt | awk -F'\t' -v OFS="\t" 'int($10)==0' | awk -F'\t' -v OFS="\t" '{print $2,$3,$4,$1,$7}' > $DATAPATH/superenhancers/Enhancers_80Epi.bed
sed '1,6d' $DATAPATH/superenhancers/H3K27ac_80Epi_peaks_AllStitched.table.txt | awk -F'\t' -v OFS="\t" 'int($10)==1' | awk -F'\t' -v OFS="\t" '{print $2,$3,$4,$1,$7}' > $DATAPATH/superenhancers/Superenhancers_80Epi.bed
```

From this, we have the following files:

1. `Enhancers_Dome.bed`: Locations of identified enhancers at *Dome* stage.
2. `Superenhancers_Dome.bed`: Locations of identified superenhancers at *Dome* stage.
3. `Enhancers_80Epi.bed`: Locations of identified enhancers at *80% Epiboly* stage.
4. `Superenhancers_80Epi.bed`: Locations of identified enhancers at *80% Epiboly* stage.
5. `H3K27ac_Dome_peaks_SuperStitched_GENE_TO_REGION.txt`: List of overlapping and proximal superenhancers (within 50 kb) for genes at *Dome*.
6. `H3K27ac_Dome_peaks_SuperStitched_REGION_TO_GENE.txt`: List of overlapping and proximal genes (within 50 kb) for identified superenhancers at *Dome*.
7. `H3K27ac_80Epi_peaks_SuperStitched_GENE_TO_REGION.txt`: List of overlapping and proximal superenhancers (within 50 kb) for genes at *80% Epiboly*.
6. `H3K27ac_80Epi_peaks_SuperStitched_REGION_TO_GENE.txt`: List of overlapping and proximal genes (within 50 kb) for identified superenhancers at *80% Epiboly*.


### Making coverage tracks around superenhancers
We use the package [**pyGenomeTracks**](https://pygenometracks.readthedocs.io/en/latest/) to make coverage tracks around the identified superenhancer regions. First, we make **BED** files of a 1 Mb region around the identified superenhancers at *Dome* and *80% Epiboly*.

```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/danRer10/bigZips/danRer10.chrom.sizes -P $DATAPATH/genome_files/
awk -v OFS="\t" 'FNR==NR{a[$1]=$2;next} ($1 in a) {print $1,$2,$3,$4,$5,a[$1]}' $DATAPATH/genome_files/danRer10.chrom.sizes $DATAPATH/superenhancers/Superenhancers_Dome.bed | awk -F'\t' -v OFS="\t" '{print $1,(int((int($2)+int($3))/2)-50000 > 0 ? int((int($2)+int($3))/2)-50000: 1),(int((int($2)+int($3))/2)+50000 < $6 ? int((int($2)+int($3))/2)+50000: $6),$4,$5}' > $DATAPATH/superenhancers/Superenhancers_Dome_trackregions.bed
awk -v OFS="\t" 'FNR==NR{a[$1]=$2;next} ($1 in a) {print $1,$2,$3,$4,$5,a[$1]}' $DATAPATH/genome_files/danRer10.chrom.sizes $DATAPATH/superenhancers/Superenhancers_80Epi.bed | awk -F'\t' -v OFS="\t" '{print $1,(int((int($2)+int($3))/2)-50000 > 0 ? int((int($2)+int($3))/2)-50000: 1),(int((int($2)+int($3))/2)+50000 < $6 ? int((int($2)+int($3))/2)+50000: $6),$4,$5}' > $DATAPATH/superenhancers/Superenhancers_80Epi_trackregions.bed
```
`pyGenomeTracks`, the program that makes the coverage track files requires the **BED** files made above, and it produces one file per line (or region) in the **BED** file. It also requires a **.ini** file which gives the details of all the tracks to be plotted. Put the provided **coverage_tracks.ini** file in the `coverage_tracks/` folder.

```bash
mdkir $DATAPATH/coverage_tracks/SE_Dome $DATAPATH/coverage_tracks/SE_80Epi
pyGenomeTracks --tracks $DATAPATH/coverage_tracks/coverage_tracks.ini --BED $DATAPATH/superenhancers/Superenhancers_Dome_trackregions.bed --outFileName $DATAPATH/coverage_tracks/SE_Dome/tracks.pdf
pyGenomeTracks --tracks $DATAPATH/coverage_tracks/coverage_tracks.ini --BED $DATAPATH/superenhancers/Superenhancers_80Epi_trackregions.bed --outFileName $DATAPATH/coverage_tracks/SE_80Epi/tracks.pdf
```
From these, we selected a few regions to design Oligopaint probes targetting them. 

```bash
awk -v OFS="\t" 'FNR==NR{a[$1]=$2;next} ($1 in a) {print $1,$2,$3,$4,$5,a[$1]}' $DATAPATH/genome_files/danRer10.chrom.sizes $DATAPATH/coverage_tracks/OP_coords_sorted.bed | awk -F'\t' -v OFS="\t" '{print $1,(int((int($2)+int($3))/2)-100000 > 0 ? int((int($2)+int($3))/2)-100000: 1),(int((int($2)+int($3))/2)+100000 < $6 ? int((int($2)+int($3))/2)+100000: $6),$4,$5}' > $DATAPATH/coverage_tracks/OP_trackregions.bed 
mdkir $DATAPATH/coverage_tracks/OP_regions
pyGenomeTracks --tracks $DATAPATH/coverage_tracks/coverage_tracks_OP.ini --BED $DATAPATH/coverage_tracks/OP_trackregions.bed --outFileName $DATAPATH/coverage_tracks/OP_regions/tracks.pdf
```

### Making Fig. X

Making genomic windows of 50 kb:

```bash
bedtools makewindows -g $DATAPATH/genome_files/danRer10.chrom.sizes -w 50000 -s 25000 > $DATAPATH/genome_files/danRer10_windows_w50k_s25k.bed
sort -k1,1 -k2,2n $DATAPATH/genome_files/danRer10_windows_w50k_s25k.bed > $DATAPATH/genome_files/danRer10_windows_w50k_s25k_sorted.bed
```
Keeping all duplicates: 

```bash
bedtools multicov -bams $DATAPATH/aligned_data/Ser5P_Dome_ut1m_sorted.bam $DATAPATH/aligned_data/Input_Dome_ut1m_sorted.bam $DATAPATH/aligned_data/H3K27ac_Dome_ut1m_sorted.bam -bed $DATAPATH/genome_files/danRer10_windows_w50k_s25k_sorted.bed > $DATAPATH/FigX/Genomic_50k_windows.txt
bedtools multicov -bams $DATAPATH/aligned_data/Ser5P_Dome_ut1m_sorted.bam $DATAPATH/aligned_data/Input_Dome_ut1m_sorted.bam $DATAPATH/aligned_data/H3K27ac_Dome_ut1m_sorted.bam -bed $DATAPATH/superenhancers/Superenhancers_Dome_sorted.bed > FigX/Superenhancers_Dome.txt
bedtools multicov -bams $DATAPATH/aligned_data/Ser5P_Dome_ut1m_sorted.bam $DATAPATH/aligned_data/Input_Dome_ut1m_sorted.bam $DATAPATH/aligned_data/H3K27ac_Dome_ut1m_sorted.bam -bed $DATAPATH/coverage_tracks/OP_coords_sorted.bed > FigX/OP.txt
awk -v OFS="\t" '{print $1,$2,$3,$4*1251132686/(224*35513591),$5*1195445591/(191*23824857),$6*1195445591/(287*6256116)}' $DATAPATH/FigX/Genomic_50k_windows.txt | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,($5>0 ? $4/$5: 0), 1}' > $DATAPATH/FigX/Genomic_50k_windows_normalized.txt
awk -v OFS="\t" '{print $1,$2,$3,$6*1251132686/(224*35513591),$7*1195445591/(191*23824857),$8*1195445591/(287*6256116)}' $DATAPATH/FigX/Superenhancers_Dome.txt | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,($5>0 ? $4/$5: 0), 2}' > $DATAPATH/FigX/Superenhancers_Dome_normalized.txt
awk -v OFS="\t" '{print $1,$2,$3,$5*1251132686/(224*35513591),$6*1195445591/(191*23824857),$7*1195445591/(287*6256116),$4}' $DATAPATH/FigX/OP.txt | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,($5>0 ? $4/$5: 0), $7}' > $DATAPATH/FigX/OP_normalized.txt 
```

After removing duplicates: 

```bash
bedtools multicov -bams $DATAPATH/aligned_data/Ser5P_Dome_ut1m_dr.bam $DATAPATH/aligned_data/Input_Dome_ut1m_sorted.bam $DATAPATH/aligned_data/H3K27ac_Dome_ut1m_sorted.bam -bed $DATAPATH/genome_files/danRer10_windows_w50k_s25k_sorted.bed > $DATAPATH/FigX/Genomic_50k_windows_dr.txt
bedtools multicov -bams aligned_data/Ser5P_Dome_ut1m_dr.bam aligned_data/Input_Dome_ut1m_sorted.bam aligned_data/H3K27ac_Dome_ut1m_sorted.bam -bed superenhancers/Superenhancers_Dome_sorted.bed > FigX/Superenhancers_Dome_dr.txt
bedtools multicov -bams aligned_data/Ser5P_Dome_ut1m_dr.bam aligned_data/Input_Dome_ut1m_sorted.bam aligned_data/H3K27ac_Dome_ut1m_sorted.bam -bed coverage_tracks/OP_coords_sorted.bed > FigX/OP_dr.txt
awk -v OFS="\t" '{print $1,$2,$3,$6*1251132686/(224*10681375),$7*1195445591/(191*23824857),$8*1195445591/(287*6256116)}' FigX/Superenhancers_Dome_dr.txt | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,($5>0 ? $4/$5: 0), 2}' > FigX/Superenhancers_Dome_dr_normalized.txt
awk -v OFS="\t" '{print $1,$2,$3,$4*1251132686/(224*10681375),$5*1195445591/(191*23824857),$6*1195445591/(287*6256116)}' FigX/Genomic_50k_windows_dr.txt | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,($5>0 ? $4/$5: 0), 1}' > FigX/Genomic_50k_windows_dr_normalized.txt 
awk -v OFS="\t" '{print $1,$2,$3,$5*1251132686/(224*10681375),$6*1195445591/(191*23824857),$7*1195445591/(287*6256116),$4}' $DATAPATH/FigX/OP_dr.txt | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,($5>0 ? $4/$5: 0), $7}' > $DATAPATH/FigX/OP_dr_normalized.txt 
```

```bash
pyGenomeTracks --tracks $DATAPATH/FigX/coverage_tracks.ini --BED $DATAPATH/FigX/FigX_regions.bed --outFileName $DATAPATH/FigX/FigX.pdf
```

### Making Fig. Y

```bash
awk -v OFS="\t" '{print $0,($2>$3? $2-$3: $3-$2)}' danRer10_genes.bed | sort -k4,4 -k7,7rn | sort -uk4,4 | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}' | sort -k 1,1 -k2,2n > danRer10_genes_unique_sorted.bed
```

```bash
awk -v OFS="\t" '{print $1,($6=="+" ? int($2)-2000: int($3)-200), ($6=="+" ? int($2)+200: int($3)+2000),$4,$5,$6}' $DATAPATH/genome_files/danRer10_genes_sorted.bed > $DATAPATH/genome_files/danRer10_genes_promoters.bed
awk -v OFS="\t" 'FNR==NR{a[$1]=$2;next} ($1 in a) {print $1,$2,$3,$4,$5,$6,a[$1]}' $DATAPATH/genome_files/danRer10.chrom.sizes genome_files/danRer10_genes_promoters.bed | awk -F'\t' -v OFS="\t" '{print $1,(int($2) > 0 ? int($2): 1),(int($3) < $7 ? int($3): $7),$4,$5,$6}' > $DATAPATH/genome_files/danRer10_promoters.bed
awk -v OFS="\t" '{print $1,($6=="+" ? int($2)+200: int($2)), ($6=="+" ? int($3): int($3)-200),$4,$5,$6}' $DATAPATH/genome_files/danRer10_genes_sorted.bed > $DATAPATH/genome_files/danRer10_genes_genebody.bed
awk -v OFS="\t" 'FNR==NR{a[$1]=$2;next} ($1 in a) {print $1,$2,$3,$4,$5,$6,a[$1]}' $DATAPATH/genome_files/danRer10.chrom.sizes $DATAPATH/genome_files/danRer10_genes_genebody.bed | awk -F'\t' -v OFS="\t" '{print $1,(int($2) > 0 ? int($2): 1),(int($3) < $7 ? int($3): $7),$4,$5,$6}' | awk -F'\t' -v OFS="\t" '{if(int($2)<int($3)) print $0}'> $DATAPATH/genome_files/danRer10_genebody.bed
```

```bash
bedtools multicov -bams $DATAPATH/aligned_data/Ser5P_Dome_ut1m_sorted_dr.bam $DATAPATH/aligned_data/Input_Dome_ut1m_sorted.bam $DATAPATH/aligned_data/H3K27ac_Dome_ut1m_sorted.bam -bed $DATAPATH/genome_files/danRer10_promoters.bed > $DATAPATH/FigY/Promoters.txt
bedtools multicov -bams $DATAPATH/aligned_data/Ser5P_Dome_ut1m_sorted_dr.bam $DATAPATH/aligned_data/Input_Dome_ut1m_sorted.bam $DATAPATH/aligned_data/H3K27ac_Dome_ut1m_sorted.bam -bed $DATAPATH/genome_files/danRer10_genebody.bed > $DATAPATH/FigY/Genebody.txt
awk -v OFS="\t" '{print $1,$2,$3,$7*1251132686/(224*10681375),$8*1195445591/(191*23824857),$9*1195445591/(287*6256116),$4}' $DATAPATH/FigY/Promoters.txt | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,($5>0 ? $4/$5: 0), $7}' > $DATAPATH/FigY/Promoters_normalized.txt 
awk -v OFS="\t" '{print $1,$2,$3,$7*1251132686/(224*10681375),$8*1195445591/(191*23824857),$9*1195445591/(287*6256116),$4}' $DATAPATH/FigY/Genebody.txt | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,($5>0 ? $4/$5: 0), $7}' > $DATAPATH/FigY/Genebody_normalized.txt 
```

```bash
pyGenomeTracks --tracks $DATAPATH/FigY/coverage_tracks.ini --BED $DATAPATH/FigY/FigY_regions.bed --outFileName $DATAPATH/FigY/FigY.pdf --width 13 --height 4
```

## Oligopaint design

We obtained genome region homology oligos for danRer10 from https://oligopaints.hms.harvard.edu/genome-files/danrer10. We used the `balance` oligo set because it resulted in the highest probe density. Sort these BED files and place these in the `$DATAPATH/oligopaints/hybridizing_probes` folder.

We used the following design to make oligopaint probes:

UniversalMainstreet-LocusSpecificStreet-UniversalFlourophoreMainstreet-HybridizingOligo-UniversalBackstreet

Primers against UniversalMainstreet and UniversalBackstreet are used to amplify the whole Oligo library, primers against LocusSpecificStreet can be used to amplify a specific oligopaint locus from the oligopaint library, UniversalFlourophoreMainstreet binds to the flourophore containing primer that we use to image the sample, and HybridizingOligo contains the gene region homology oligos obtained above. 

We used `OligoLego` https://github.com/gnir/OligoLego to design streets and pick sets of streets that can function well together. First we used `MakingStreets.m` from `OligoLego` to design a large list of possible streets.  

The list of genomic regions for which we designed oligopaint probes is available in the BED file `OP.bed`. Using this and the hybridizing probe set, we run `bedtools intersect` to obtain `Main_isected.txt` which contains the coordinates and sequences of hybridizing probes that overlap each region in `OP.bed`. This file should contain 9 columns as in the file provided. 

One can use `ApOPs.m` from `OligoLego` to design the oligopaint library. However, because we use a different design from what is offered in `OligoLego`, we instead use our own program `AppendGivenStreets.m` and supply it with a list of streets `GivenStreets.txt` (selected using the same criteria as `OligoLego`'s programs), and the universal streets in `Universal.txt` to construct the library from. This program results in a list of files that we can use to order the oligopaint library. 

`MSDensity.txt` contains information about the oligopaint regions -- coordinates of each region, number of probes for each region, and the per-region density of probes. `PrimersToOrder.txt` contains the list of universal and locus-specific probes to be ordered. `Oligopaints.txt` contains the final oligopaint library with each oligopaint probe and its sequence. We use `sed 's/ \+/,/g' Oligopaints.txt > Oligopaints.csv` to prepare a CSV file to order from.
