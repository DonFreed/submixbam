#!/usr/bin/env bash

REFERENCE=$1

# Download the fastq files from NA12878 and NA12882 using [ascp](https://www.ebi.ac.uk/ega/about/ftp-aspera). These files are large.

mkdir Fastq

ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174324/ERR174324_1.fastq.gz ./Fastq/
ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174324/ERR174324_2.fastq.gz ./Fastq/

ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174325/ERR174325_1.fastq.gz ./Fastq/
ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174325/ERR174325_2.fastq.gz ./Fastq/

ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174326/ERR174326_1.fastq.gz ./Fastq/
ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174326/ERR174326_2.fastq.gz ./Fastq/

ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174327/ERR174327_1.fastq.gz ./Fastq/
ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174327/ERR174327_2.fastq.gz ./Fastq/

ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174328/ERR174328_1.fastq.gz ./Fastq/
ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174328/ERR174328_2.fastq.gz ./Fastq/

ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174329/ERR174329_1.fastq.gz ./Fastq/
ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174329/ERR174329_2.fastq.gz ./Fastq/

ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174330/ERR174330_1.fastq.gz ./Fastq/
ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174330/ERR174330_2.fastq.gz ./Fastq/

ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174331/ERR174331_1.fastq.gz ./Fastq/
ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174331/ERR174331_2.fastq.gz ./Fastq/

ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174332/ERR174332_1.fastq.gz ./Fastq/
ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174332/ERR174332_2.fastq.gz ./Fastq/

ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174333/ERR174333_1.fastq.gz ./Fastq/
ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174333/ERR174333_2.fastq.gz ./Fastq/

ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174334/ERR174334_1.fastq.gz ./Fastq/
ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174334/ERR174334_2.fastq.gz ./Fastq/

ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174335/ERR174335_1.fastq.gz ./Fastq/
ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174335/ERR174335_2.fastq.gz ./Fastq/

ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174336/ERR174336_1.fastq.gz ./Fastq/
ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174336/ERR174336_2.fastq.gz ./Fastq/

ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174337/ERR174337_1.fastq.gz ./Fastq/
ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174337/ERR174337_2.fastq.gz ./Fastq/

ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174338/ERR174338_1.fastq.gz ./Fastq/
ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174338/ERR174338_2.fastq.gz ./Fastq/

ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174347/ERR174347_1.fastq.gz ./Fastq/
ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174347/ERR174347_2.fastq.gz ./Fastq/

ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174348/ERR174348_1.fastq.gz ./Fastq/
ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174348/ERR174348_2.fastq.gz ./Fastq/

ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174349/ERR174349_1.fastq.gz ./Fastq/
ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174349/ERR174349_2.fastq.gz ./Fastq/

ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174350/ERR174350_1.fastq.gz ./Fastq/
ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174350/ERR174350_2.fastq.gz ./Fastq/

ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174351/ERR174351_1.fastq.gz ./Fastq/
ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174351/ERR174351_2.fastq.gz ./Fastq/

ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174368/ERR174368_1.fastq.gz ./Fastq/
ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174368/ERR174368_2.fastq.gz ./Fastq/

ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174369/ERR174369_1.fastq.gz ./Fastq/
ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174369/ERR174369_2.fastq.gz ./Fastq/

ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174370/ERR174370_1.fastq.gz ./Fastq/
ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174370/ERR174370_2.fastq.gz ./Fastq/

ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174371/ERR174371_1.fastq.gz ./Fastq/
ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174371/ERR174371_2.fastq.gz ./Fastq/

ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174372/ERR174372_1.fastq.gz ./Fastq/
ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174372/ERR174372_2.fastq.gz ./Fastq/

ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174373/ERR174373_1.fastq.gz ./Fastq/
ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174373/ERR174373_2.fastq.gz ./Fastq/

ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174374/ERR174374_1.fastq.gz ./Fastq/
ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174374/ERR174374_2.fastq.gz ./Fastq/

ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174375/ERR174375_1.fastq.gz ./Fastq/
ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174375/ERR174375_2.fastq.gz ./Fastq/

ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174376/ERR174376_1.fastq.gz ./Fastq/
ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174376/ERR174376_2.fastq.gz ./Fastq/

ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174377/ERR174377_1.fastq.gz ./Fastq/
ascp -T -l 400M -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR174/ERR174377/ERR174377_2.fastq.gz ./Fastq/


# Align these files using BWA #
mkdir BAM
for fq1 in `ls ./Fastq/*[23][0-9]_1.fastq.gz`
do
    run=$(basename ${fq1%%_1.fastq.gz})
    fq2=${fq1%%_1.fastq.gz}_2.fastq.gz
    bwa mem $REFERENCE -R \'@RG\tID:${run}\tSM:NA12878\' -t 8 $fq1 $fq2 | samblaster | samtools sort -o BAM/${run}_sorted.bam -T BAM/tmp_sort -O BAM
done

for fq1 in `ls ./Fastq/*[4567][0-9]_1.fastq.gz`
do
    run=$(basename ${fq1%%_1.fastq.gz})
    fq2=${fq1%%_1.fastq.gz}_2.fastq.gz
    bwa mem $REFERENCE -R \'@RG\tID:${run}\tSM:NA12878\' -t 8 $fq1 $fq2 | samblaster | samtools sort -o BAM/${run}_sorted.bam -T BAM/tmp_sort -O BAM
done

# Merge the aligned BAM files #
samtools merge BAM/NA12878_merged_sorted.bam BAM/*[23][0-9]_sorted.bam
samtools merge BAM/NA12882_merged_sorted.bam BAM/*[4567][0-9]_sorted.bam

# Index #
samtools index BAM/NA12878_merged_sorted.bam
samtools index BAM/NA12882_merged_sorted.bam

# Grab the BAM headers #
# Change NA12882 to NA12878 for post-mixing variant calling #
mkdir Tmp
samtools view -H BAM/NA12878_merged_sorted.bam > Tmp/header.sam
samtools view -H BAM/NA12882_merged_sorted.bam | grep "^RG\|^PG" | sed 's/NA12882/NA12878/g' >> Tmp/header.sam

samtools idxstats BAM/NA12878_merged_sorted.bam | awk 'BEGIN{SUM=0} {SUM += $3} END {print SUM}'
# 6218931678
samtools idxstats BAM/NA12882_merged_sorted.bam | awk 'BEGIN{SUM=0} {SUM += $3} END {print SUM}'
# 6615728881

# 6218931678 * 101 / 3100000000 = ~203x
# 6615728881 * 101 / 3100000000 = ~216x

# 30 / 203 = .1478
# 30 / 216 = .1389

# 50 / 203 = .2463
# 50 / 216 = .2315

# Download GIAB regions to mix over #
mkdir BED
cd BED
wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.2.2/NA12878_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2.2_highconf.bed
cd ..

# Merge nearby GIAB regions #
bedtools merge -d 150 BED/NA12878_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2.2_highconf.bed > BED/giab_150bp_merge.bed

# Subsample to 30x with random mixing #
submixbam -d 0.1478 -e 0.1389 -h Tmp/header.sam -o BAM/plat_30x_mix.bam -b BED plat_30x_midx.bed BAM/NA12878_merged_sorted.bam BAM/NA12882_merged_sorted.bam BED/giab_150bp_merge.bed

# Subsample to 50x with random mixing #
submixbam -d 0.2463 -e 0.2315 -h Tmp/header.sam -o BAM/plat_50x_mix.bam -b BED plat_50x_midx.bed BAM/NA12878_merged_sorted.bam BAM/NA12882_merged_sorted.bam BED/giab_150bp_merge.bed


