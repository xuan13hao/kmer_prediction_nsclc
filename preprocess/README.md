# DOWNLOAD
If you have installed SRA Toolkit `https://github.com/ncbi/sra-tools` in order to use `download.sh`.
you need to change FastqDUMP to /path/to/your/fastq-dump/bin
find it from the directory of your sratoolkit/bin/

Alternatively, you can use any method you prefer to download the raw data from the NCBI bio projects using accession number: 
1. PRJNA866654 (DS1) `https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA866654&o=acc_s%3Aa`
1. PRJEB22863 (DS2) `https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJEB22863&o=acc_s%3Aa`
1. PRJNA751792 (DS3) `https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA751792&o=acc_s%3Aa`

Note:
 - not all samples from DS2 are related to NSCLC. We only used the samples that were taken at timepoint T0.
 - It requires more than **2 TeraByte** hard disk to host the raw data. 

# ANNOTATION
In order to run our annotation pipeline, you should install the following prerequisites and change the path to each software:
1. FastQC (v0.11.8)
1. Trimmomatic (0.38)
1. Fiona (0.2.10)
1. BWA (0.7.16a-r1185-dirty)
1. MetaPhlAn (version 4.0.3)
1. UProC (1.2.0)
1. KMC3 (3.2.1)

Before running our annotation pipeline, you should perform quality control, removal of human contamination and error correction.

## Quality Control
1. quality score before adapter removal:
 - For DS1, using the following as it's paired-end experiment, assuming
  - > /path/to/fastqc -t 16 $SAMPLE_ID.read1.fastq $SAMPLE_ID.read2.fastq
 - For DS2 and DS3: 
  - > /path/to/fastqc -t 16 $SAMPLE_ID.fastq

2. adapter removal:
 - we knew the adapter sequences used for DS1 and we found out that the per-base sequence quality of about 20 nt from the beginning is not ideal.
    > java -jar /path/to/trimmomatic-0.38.jar PE $SAMPLE_ID.read1.fastq $SAMPLE_ID.read2.fastq $SAMPLE_ID.read1.paired.fq $SAMPLE_ID.read1.unpaired.fq $SAMPLE_ID.read2.fq $SAMPLE_ID.read2.unpaired.fq ILLUMINACLIP:SRR20881982.fa:2:30:10  LEADING:30 TRAILING:30 SLIDINGWINDOW:4:30 MINLEN:35 -threads 16;

 - we don't know the adapter sequences and we found out that the per-base sequence quality after 280 nt drops in DS2 and DS3.
    > java -jar /path/to/trimmomatic-0.38.jar SE $SAMPLE_ID.fastq $SAMPLE_ID.trimed.fastq ILLUMINACLIP:/path/to/trimmomatic/adapters/TruSeq2-SE.fa:2:30:10 CROP:280 LEADING:30 TRAILING:30 SLIDINGWINDOW:4:15 MINLEN:35 -threads 16;

3. use fastqc to verify that the per-base quality is better.

## Removal of Human Contamination
1. Index human genome:
   > /path/to/bwa index hg38.fasta

1. For DS1:
   - > /path/to/bwa -t 16 hg38.fa $SAMPLE_ID.read1.paired.fq $SAMPLE_ID.read2.paired.fq > $SAMPLE_ID.sam
   - > python3 removeReadSAM.py $SAMPLE_ID.read1.paired.fq $SAMPLE_ID.sam > $SAMPLE_ID.read1.paired.clean.fq
   - > python3 removeReadSAM.py $SAMPLE_ID.read2.paired.fq $SAMPLE_ID.sam > $SAMPLE_ID.read2.paired.clean.fq
  
1. For DS2 and DS3:
   - > /path/to/bwa -t 16 hg38.fa $SAMPLE_ID.trimed.fastq > $SAMPLE_ID.sam
   - > python3 removeReadSAM.py $SAMPLE_ID.trimed.fastq $SAMPLE_ID.sam > $SAMPLE_ID.trimed.clean.fq

## Error Correction
We recommend performing an error correction on the reads after the removal of human contamination:
> /path/to/fiona -nt 16 -g 4639675 $SAMPLE_ID.trimed.clean.fq $SAMPLE_ID.trimed.clean.ec.fq

## Annotation
After that, It's as simple as `bash annotate.sh $FILENAMEs $SAMPLE_ID` to run our annotation pipeline, where 
 - $FILENAMEs is the `"$SAMPLE_ID.read1.paired.clean.fq $SAMPLE_ID.read2.paired.clean.fq"` for DS1 and `$SAMPLE_ID.trimed.clean.fq` for DS2 and DS3.
 - $SAMPLE_ID is the ID you want to use to label the sample, you can just use the sample id when downloading the data
 - This script will create three folders:
   - `otu` for MetaPhlAn annotation of taxonomical profile,
   - `pfam` for UproC annotation of functional profile on Pfam database, and
   - `K30` for 30-mer counting file.
 - This script will automatically annotate one sample and put the corresponding profiles in its folder.
   
Note:
 - Be aware that for DS1, `$SAMPLE_ID.read1.paired.clean.fq $SAMPLE_ID.read2.paired.clean.fq` must be inside **""**!
 - you should be familiar with the supplemental materials to fully understand our annotation pipeline.
 - The 30-mer counting profiles of all 417 samples require more than **2 TeraByte** hard disk.

## Generate OTU-like Table
Now you can merge the annotation profiles into one OTU-like table for study. We included the scripts:
 1. getMTPTable.py
    - Python script that merges all profiles in one folder and makes an OTU-like table of the taxonomic profile of all samples at the species level.
    - Call it by `python3 getMTPTable.py metadata.csv otu`, the output is "otu.count.csv".
    - The "otu.count.csv" has the same format described in `data/README.md`.
    - You can use `python3 getMTPTable.py` to check its help information.
    - You need to install pandas to use getMTPTable.py.
 1. getUproCTable.py
    - Python script that merges all profiles in one folder and makes an OTU-like table of the functional profile generated by UproC.
    - Call it by `python3 getUproCTable.py metadata.csv pfam`, the output is "pfam.count.tsv".
    - The "pfam.count.tsv" has the same format described in `data/README.md`.
    - You can use `python3 getUproCTable.py` to check its help information.
    - You need to install pandas to use getUproCTable.py.
 1. mergeKMCTable.zip
    - The C++ code to merge all the 30-mer counting files. You need to install the boost library to use it.
    - We provided the makefile, you can just `make` at its directory.
    - After compiling, you should move "getKMCTable.exe" and the "KMC_list.txt" to the `K30` folder.
    - At the `K30` folder, run the getKMCTable.exe by `./getKMCTable.exe KMC_list.txt`, the final 30-mer table is called "k30.count.tsv", which is described in `data/README.md`
    - The mergeMCTable.exe implemented a naive two-way external merging algorithm. You can modify the procedure according to the resources you have to speed it up.

