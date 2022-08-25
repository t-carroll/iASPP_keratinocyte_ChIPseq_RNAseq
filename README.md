# Recreating analyses for JNK-iASPP-AP1 axis paper (Al Moussawi et al. *Cell Reports*, 2022)
## Recreate files needed for downstream analysis
First, download (untrimmed) RNA-seq and ChIP-seq fastq files from GEO repositories into your analysis directory, BASE_DIR:

    cd $BASE_DIR
    ##download raw FASTQs for each accesion in GSE188448
    ##use your preferred client, e.g. fastq-dump --split-files $ACCESSION

### RNA-seq processing
Download RNA reference files (GENCODE primary assembly recommended for STAR. I use GENCODE relase M25 here).

    cd $BASE_DIR/RNA_REF
    wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.primary_assembly.annotation.gtf.gz
    gunzip gencode.vM25.primary_assembly.annotation.gtf.gz
    wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.primary_assembly.genome.fa.gz
    gunzip GRCm38.primary_assembly.genome.fa.gz

Generate STAR (2.7.3a) index for 42bp reads 

    STAR --runMode genomeGenerate --runThreadN 16 --genomeDir STAR_mm10_gencode_42 --genomeFastaFiles GRCm38.primary_assembly.genome.fa --sjdbGTFfile gencode.vM25.primary_assembly.annotation.gtf --sjdbOverhang 41

Trim FASTQs from BASE_DIR with cutadapt (v2.10), output trimmed FASTQs into RNA_WORK

    cd $BASE_DIR/RNA_WORK
    cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -o $RNA_NAME.fastq.gz -j 6 -m 30 ../$RNA_NAME.fastq.gz

Run STAR (2.7.3a)

    STAR --genomeDir $BASE_DIR/RNA_REF/STAR_mm10_gencode_42 --outSAMtype BAM Unsorted --readFilesCommand zcat --readFilesIn $RNA_NAME.fastq.gz --outFileNamePrefix $RNA_NAME

Create counts tables with featureCounts (v2.0.0)

    featureCounts -T 6 -t exon -g gene_name -a $RNA_NAME/gencode.vM25.primary_assembly.annotation.gtf -s2 -o MKC_iASPP_rawCounts.txt *.bam

### ChIP-seq processing and MAnorm/motif calls
Download ChIP reference files (I used UCSC version of mm10 from iGenomes here, as ChIP-seq work was conducted before initiation of RNA-seq analysis; GENCODE primary assembly used for RNA-seq may also be suitable). 

    cd $BASE_DIR/CHIP_REF
    wget http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Mus_musculus/UCSC/mm10/Mus_musculus_UCSC_mm10.tar.gz
    tar -zxvf Mus_musculus_UCSC_mm10.tar.gz
    wget https://raw.githubusercontent.com/t-carroll/Blacklist/master/lists/mm10-blacklist.v2.bed
    wget --no-check-certificate jaspar2018.genereg.net/download/data/2018/CORE/JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt
    grep "MOTIF" JASPAR2018_CORE_vertebrates_non-redundant.meme | cut -d " " -f2 > names.jaspar2018.core.vertebrate.txt

Trim raw ChIP-seq files in BASE_DIR with cutadapt (v2.10), output into CHIP_WORK

    cd $BASE_DIR/CHIP_WORK
    cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o $NAME_R1.fastq.gz -p $NAME_R2.fastq.gz -j 6 -m 30 ../$NAME_R1.fastq.gz ../$NAME_R2.fastq.gz

Align with bwa (v0.7.17-r1188)

    bwa aln -R 2 -q 20 -t 2 $BASE_DIR/CHIP_REF/Mus_musculus/UCSC/mm10/Sequence/BWAIndex/genome.fa $NAME_R1.fastq.gz > $NAME_R1.sai
    bwa aln -R 2 -q 20 -t 2 $BASE_DIR/CHIP_REF/Mus_musculus/UCSC/mm10/Sequence/BWAIndex/genome.fa $NAME_R2.fastq.gz > $NAME_R2.sai
    bwa sampe $BASE_DIR/CHIP_REF/Mus_musculus/UCSC/mm10/Sequence/BWAIndex/genome.fa $NAME_R1.sai $NAME_R2.sai $NAME_R1.fastq.gz $NAME_R2.fastq.gz > $NAME.sam

Low-quality alignments, duplicates, and alignments to mitochondrial, random, and unplaced contigs are filtered out prior to MACS2

    samtools fixmate -@ 2 -O sam -m $NAME.sam - | samtools sort -O sam -@ 2 -o - - | samtools markdup -d 100 -@ 2 -s -f $NAME.stats.txt - - | samtools view -h -f 3 -F 4 -F 8 -F 256 -F 512 -F 1024 -F 2048 - | awk '$1~"@" || $7 == "="' | sed '/chrM/d;/random/d;/chrUn/d' | samtools sort -@ 2 -o $NAME_dedup.bam -
    samtools index $NAME_dedup.bam

Run macs2 on filtered BAM file (for further downstream analysis)

    macs2 callpeak -B -g mm -f BAMPE -s 43 --call-summits -q 0.05 -t $CHIP -c $CONTROL -n $OUTPUT_NAME
Run macs2 with SPMR (for normalized signal files, output bedgraphs can be converted to bigwigs in UCSC genome browser session with e.g. bedGraphToBigWig):

    macs2 callpeak -B --SPMR -g mm -f BAMPE -s 43 --call-summits -q 0.05 -t $CHIP -c $CONTROL -n $OUTPUT_NAME

Filter out blacklisted peaks (from non-SPMR macs2 call):

    touch $OUTPUT_NAME_peaks.filt.xls
    grep "#" $OUTPUT_NAME_peaks.xls >> $OUTPUT_NAME_peaks.filt.xls
    grep -v "#" $OUTPUT_NAME_peaks.xls | head -2 | tail -1 >> $OUTPUT_NAME_peaks.filt.xls
    grep -v "#" $OUTPUT_NAME_peaks.xls | grep -v "pileup" | bedtools intersect -a stdin -b $BASE_DIR/CHIP_REF/mm10-blacklist.v2.bed -v >> $OUTPUT_NAME_peaks.filt.xls

Run MAnorm (v1.3.0) with p63 ChIP BAM and filt_peaks.xls files for each sample (Sample 1: iASPP KO, sample 2: iASPP WT):
*MAnorm can be found at https://github.com/shao-lab/MAnorm. Using their default install is recommended, although I use my non-standard fork (https://github.com/t-carroll/MAnorm) which also output unnormalized M-A values.*

    manorm --p1 $PEAK_KO_FILT_PEAKS --p2 $PEAK_WT_FILT_PEAKS --r1 $BAM_KO --r2 $BAM_WT --n1 iASPP_KO --n2 iASPP_WT --pf macs2 --rf bam --pe -w 1000 --summit-dis 100 -o .

Get FASTA for all peak regions in MAnorm output:

    cat iASPP_KO_vs_iASPP_WT_all_MAvalues.xls | cut -f1-3 | grep -v end | bedtools getfasta -fi $BASE_DIR/CHIP_REF/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa -bed - -fo out.fasta

Get FIMO results for each motif entry in JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt
*Recommend to use your HPC job submission with $BASE_DIR/CHIP_REF/names.jaspar2018.core.vertebrate.txt to submit one job for each motif in this list. But can use this while loop to run sequentially.*

    mkdir fimo && cd fimo
    while read a; do fimo --text --skip-matched-sequence --thresh 1 --verbosity 1 --max-strand --motif $a $BASE_DIR/CHIP_REF/JASPAR2018_CORE_vertebrates_non-redundant.meme.txt ../out.fasta | awk 'NR>1' | sort -k7,7nr | awk -F "\t" '!seen[$3]++' > $a.txt; done < $BASE_DIR/CHIP_REF/names.jaspar2018.core.vertebrate.txt
    ##combine into fimo_out.txt file
    while read a; do cat $a.txt >> ../fimo_out.txt; done < $BASE_DIR/CHIP_REF/names.jaspar2018.core.vertebrate.txt

## Downstream analysis and figure generation
With all input files generated, we can now perform recreate downstream analysis (including RNA-seq/ChIP-seq integration) and generation of relevant figures following the code in paper.Rmd.
### STREME motif output (Fig. S2E)
To run recreate S2E, STREME needs to be run on output bed files from integrative analysis in paper.Rmd. First get sequences from selected bed entries:

    bedtools getfasta -fi $BASE_DIR/CHIP_REF/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa -bed $git_repo/output/p63_peaks_enriched_iASPP_WT.bed > $git_repo/output/p63_peaks_enriched_iASPP_WT.fa
    bedtools getfasta -fi $BASE_DIR/CHIP_REF/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa -bed $git_repo/output/p63_peaks_enriched_iASPP_KO.bed > $git_repo/output/p63_peaks_enriched_iASPP_KO.fa
    bedtools getfasta -fi $BASE_DIR/CHIP_REF/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa -bed $git_repo/output/p63_peaks_not_enriched.bed > $git_repo/output/p63_peaks_not_enriched.fa
    

Then can run STREME (v5.4.1) at the web server (https://meme-suite.org/meme/tools/streme) with each output fasta file from above. Default options are used, except max motif width of 25 and number of motifs = 1 rather than p-value cut-off (both found in "Advanced options"). Or for CLI, the following command can be used (repeat with each fasta):

    streme --verbosity 1 --oc . --dna --totallength 4000000 --time 14400 --minw 8 --maxw 25 --nmotifs 1 --align center --p p63_peaks_enriched_iASPP_WT.fa

## Citation
If using the data or analyses from this repository or the underlying publication in your work, please make sure to cite the following publication:
"Mutant Ras and inflammation-driven skin tumorigenesis is suppressed via a JNK-iASPP-AP1 axis". Al Moussawi *et al*. *Cell Reports*, 2022. 

## Need help?
Ping me a message @t-carroll, or raise an issue on the tab above!
