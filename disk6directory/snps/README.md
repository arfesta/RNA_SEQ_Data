# 1/29/18 - 2/1/18 SNP calling notes

*Raw tar files returned from GSL were extracted and names of folders were changed to correspond to a single pool (lane) run for each batch.  The .fastq.gz files for each library were QC'ed with fastqc and then adapters/contamination was removed using bbduk.  The filtered and trimmed reads are located at:

```linux
/media/disk6/ARF/RNASEQ/trimmedfiltreads
```

*The above path is the starting directory for processing SNPs.  Notes of lane extraction, fastqc, and filter/trim processing are on git at:

```linux
Breeding-Value-Prediction/disk6directory/RNASEQ/rawreads/processingnotes
```

## 1. Build index of 86k transcriptome

*Build bowtie2 index:

```linux
cd /media/disk6/ARF/RNASEQ
./software/bowtie/bowtie2-build -f ./references/Pita.86ktxptome.fasta Pita.86k
```

## 2. Align .fastq.gz techincal reps to 86K transcriptome with bowtie

*The following code aligns techincal rplicates to the 86k transcriptome and then converts to a compresed bam file

*EW samples were ran first.  A semicolon was used to string together processing of separate lanes.

```linux
while read file; do ./software/bowtie/bowtie2 -p 60 -x ./references/Pta.86k -U ./trimmedfiltreads/EW/lane01/$file | ./software/samtools-1.6/samtools view -@ 60 -b1 - -o ./snps/86k/EW/lane01/$file; done < ./rawreads/EWfasta/lane01/lane01.list;while read file; do ./software/bowtie/bowtie2 -p 60 -x ./references/Pta.86k -U ./trimmedfiltreads/EW/lane02/$file | ./software/samtools-1.6/samtools view -@ 60 -b1 - -o ./snps/86k/EW/lane02/$file; done < ./rawreads/EWfasta/lane02/lane02.list;

while read file; do ./software/bowtie/bowtie2 -p 60 -x ./references/Pta.86k -U ./trimmedfiltreads/EW/lane03/$file | ./software/samtools-1.6/samtools view -@ 60 -b1 - -o ./snps/86k/EW/lane03/$file; done < ./rawreads/EWfasta/lane03/lane03.list; while read file; do ./software/bowtie/bowtie2 -p 60 -x ./references/Pta.86k -U ./trimmedfiltreads/EW/lane04/$file | ./software/samtools-1.6/samtools view -@ 60 -b1 - -o ./snps/86k/EW/lane04/$file; done < ./rawreads/EWfasta/lane04/lane04.list; while read file; do ./software/bowtie/bowtie2 -p 60 -x ./references/Pta.86k -U ./trimmedfiltreads/EW/lane05/$file | ./software/samtools-1.6/samtools view -@ 60 -b1 - -o ./snps/86k/EW/lane05/$file; done < ./rawreads/EWfasta/lane05/lane05.list; while read file; do ./software/bowtie/bowtie2 -p 60 -x ./references/Pta.86k -U ./trimmedfiltreads/EW/lane06/$file | ./software/samtools-1.6/samtools view -@ 60 -b1 - -o ./snps/86k/EW/lane06/$file; done < ./rawreads/EWfasta/lane06/lane06.list;  while read file; do ./software/bowtie/bowtie2 -p 60 -x ./references/Pta.86k -U ./trimmedfiltreads/EW/lane07/$file | ./software/samtools-1.6/samtools view -@ 60 -b1 - -o ./snps/86k/EW/lane07/$file; done < ./rawreads/EWfasta/lane07/lane07.list;  while read file; do ./software/bowtie/bowtie2 -p 60 -x ./references/Pta.86k -U ./trimmedfiltreads/EW/lane08/$file | ./software/samtools-1.6/samtools view -@ 60 -b1 - -o ./snps/86k/EW/lane08/$file; done < ./rawreads/EWfasta/lane08/lane08.list;  while read file; do ./software/bowtie/bowtie2 -p 60 -x ./references/Pta.86k -U ./trimmedfiltreads/EW/lane09/$file | ./software/samtools-1.6/samtools view -@ 60 -b1 - -o ./snps/86k/EW/lane09/$file; done < ./rawreads/EWfasta/lane09/lane09.list;  while read file; do ./software/bowtie/bowtie2 -p 60 -x ./references/Pta.86k -U ./trimmedfiltreads/EW/lane10/$file | ./software/samtools-1.6/samtools view -@ 60 -b1 - -o ./snps/86k/EW/lane10/$file; done < ./rawreads/EWfasta/lane10/lane10.list;  while read file; do ./software/bowtie/bowtie2 -p 60 -x ./references/Pta.86k -U ./trimmedfiltreads/EW/lane11/$file | ./software/samtools-1.6/samtools view -@ 60 -b1 - -o ./snps/86k/EW/lane11/$file; done < ./rawreads/EWfasta/lane11/lane11.list;  while read file; do ./software/bowtie/bowtie2 -p 60 -x ./references/Pta.86k -U ./trimmedfiltreads/EW/lane12/$file | ./software/samtools-1.6/samtools view -@ 60 -b1 - -o ./snps/86k/EW/lane12/$file; done < ./rawreads/EWfasta/lane12/lane12.list

```

*The same process was done for the LGEP lanes:

```linux
while read file; do ./software/bowtie/bowtie2 -p 60 -x ./references/Pta.86k -U ./trimmedfiltreads/LGEP/lane01/$file | ./software/samtools-1.6/samtools view -@ 60 -b1 - -o ./snps/86k/LGEP/lane01/$file; done < ./rawreads/LGEPfasta/lane01/lane01.list;while read file; do ./software/bowtie/bowtie2 -p 60 -x ./references/Pta.86k -U ./trimmedfiltreads/LGEP/lane02/$file | ./software/samtools-1.6/samtools view -@ 60 -b1 - -o ./snps/86k/LGEP/lane02/$file; done < ./rawreads/LGEPfasta/lane02/lane02.list;while read file; do ./software/bowtie/bowtie2 -p 60 -x ./references/Pta.86k -U ./trimmedfiltreads/LGEP/lane03/$file | ./software/samtools-1.6/samtools view -@ 60 -b1 - -o ./snps/86k/LGEP/lane03/$file; done < ./rawreads/LGEPfasta/lane03/lane03.list; while read file; do ./software/bowtie/bowtie2 -p 60 -x ./references/Pta.86k -U ./trimmedfiltreads/LGEP/lane04/$file | ./software/samtools-1.6/samtools view -@ 60 -b1 - -o ./snps/86k/LGEP/lane04/$file; done < ./rawreads/LGEPfasta/lane04/lane04.list; while read file; do ./software/bowtie/bowtie2 -p 60 -x ./references/Pta.86k -U ./trimmedfiltreads/LGEP/lane05/$file | ./software/samtools-1.6/samtools view -@ 60 -b1 - -o ./snps/86k/LGEP/lane05/$file; done < ./rawreads/LGEPfasta/lane05/lane05.list; while read file; do ./software/bowtie/bowtie2 -p 60 -x ./references/Pta.86k -U ./trimmedfiltreads/LGEP/lane06/$file | ./software/samtools-1.6/samtools view -@ 60 -b1 - -o ./snps/86k/LGEP/lane06/$file; done < ./rawreads/LGEPfasta/lane06/lane06.list;  while read file; do ./software/bowtie/bowtie2 -p 60 -x ./references/Pta.86k -U ./trimmedfiltreads/LGEP/lane07/$file | ./software/samtools-1.6/samtools view -@ 60 -b1 - -o ./snps/86k/LGEP/lane07/$file; done < ./rawreads/LGEPfasta/lane07/lane07.list;  while read file; do ./software/bowtie/bowtie2 -p 60 -x ./references/Pta.86k -U ./trimmedfiltreads/LGEP/lane08/$file | ./software/samtools-1.6/samtools view -@ 60 -b1 - -o ./snps/86k/LGEP/lane08/$file; done < ./rawreads/LGEPfasta/lane08/lane08.list;  while read file; do ./software/bowtie/bowtie2 -p 60 -x ./references/Pta.86k -U ./trimmedfiltreads/LGEP/lane09/$file | ./software/samtools-1.6/samtools view -@ 60 -b1 - -o ./snps/86k/LGEP/lane09/$file; done < ./rawreads/LGEPfasta/lane09/lane09.list;  while read file; do ./software/bowtie/bowtie2 -p 60 -x ./references/Pta.86k -U ./trimmedfiltreads/LGEP/lane10/$file | ./software/samtools-1.6/samtools view -@ 60 -b1 - -o ./snps/86k/LGEP/lane10/$file; done < ./rawreads/LGEPfasta/lane10/lane10.list;  while read file; do ./software/bowtie/bowtie2 -p 60 -x ./references/Pta.86k -U ./trimmedfiltreads/LGEP/lane11/$file | ./software/samtools-1.6/samtools view -@ 60 -b1 - -o ./snps/86k/LGEP/lane11/$file; done < ./rawreads/LGEPfasta/lane11/lane11.list;  while read file; do ./software/bowtie/bowtie2 -p 60 -x ./references/Pta.86k -U ./trimmedfiltreads/LGEP/lane12/$file | ./software/samtools-1.6/samtools view -@ 60 -b1 - -o ./snps/86k/LGEP/lane12/$file; done < ./rawreads/LGEPfasta/lane12/lane12.list; 

while read file; do ./software/bowtie/bowtie2 -p 60 -x ./references/Pta.86k -U ./trimmedfiltreads/LGEP/lane13/$file | ./software/samtools-1.6/samtools view -@ 60 -b1 - -o ./snps/86k/LGEP/lane13/$file; done < ./rawreads/LGEPfasta/lane13/lane13.list;  while read file; do ./software/bowtie/bowtie2 -p 60 -x ./references/Pta.86k -U ./trimmedfiltreads/LGEP/lane14/$file | ./software/samtools-1.6/samtools view -@ 60 -b1 - -o ./snps/86k/LGEP/lane14/$file; done < ./rawreads/LGEPfasta/lane14/lane14.list;  while read file; do ./software/bowtie/bowtie2 -p 60 -x ./references/Pta.86k -U ./trimmedfiltreads/LGEP/lane15/$file | ./software/samtools-1.6/samtools view -@ 60 -b1 - -o ./snps/86k/LGEP/lane15/$file; done < ./rawreads/LGEPfasta/lane15/lane15.list;  while read file; do ./software/bowtie/bowtie2 -p 60 -x ./references/Pta.86k -U ./trimmedfiltreads/LGEP/lane16/$file | ./software/samtools-1.6/samtools view -@ 60 -b1 - -o ./snps/86k/LGEP/lane16/$file; done < ./rawreads/LGEPfasta/lane16/lane16.list;  while read file; do ./software/bowtie/bowtie2 -p 60 -x ./references/Pta.86k -U ./trimmedfiltreads/LGEP/lane17/$file | ./software/samtools-1.6/samtools view -@ 60 -b1 - -o ./snps/86k/LGEP/lane17/$file; done < ./rawreads/LGEPfasta/lane17/lane17.list;  while read file; do ./software/bowtie/bowtie2 -p 60 -x ./references/Pta.86k -U ./trimmedfiltreads/LGEP/lane18/$file | ./software/samtools-1.6/samtools view -@ 60 -b1 - -o ./snps/86k/LGEP/lane18/$file; done < ./rawreads/LGEPfasta/lane18/lane18.list; 

```

## 3. Sort the aligned techincal reps and mark/remove duplicates:

*Within each of the EW/lane(xx) directories a directory titled "sort.rmdup" was created to store finished sorted/rmdup bams

```linux
cd /media/disk6/ARF/RNASEQ/snps/86k/EW/lane01

mkdir sort.rmdup
```

*The above was done for every lane in both EW & LGEP SNP directories.  Now we can pipe together the sorting and removing of duplicates for each lane:

```linux
while read file; do ./software/samtools-1.6/samtools sort -m 5G -@ 20 ./snps/86k/EW/lane01/$file | ./software/samtools-1.6/samtools markdup -rs -@ 20 - ./snps/86k/EW/lane01/sort.rmdup/$file; done < ./rawreads/EWfasta/lane01/lane01.list

while read file; do ./software/samtools-1.6/samtools sort -m 5G -@ 20 ./snps/86k/EW/lane02/$file | ./software/samtools-1.6/samtools markdup -rs -@ 20 - ./snps/86k/EW/lane02/sort.rmdup/$file; done < ./rawreads/EWfasta/lane02/lane02.list; while read file; do ./software/samtools-1.6/samtools sort -m 5G -@ 20 ./snps/86k/EW/lane03/$file | ./software/samtools-1.6/samtools markdup -rs -@ 20 - ./snps/86k/EW/lane03/sort.rmdup/$file; done < ./rawreads/EWfasta/lane03/lane03.list; while read file; do ./software/samtools-1.6/samtools sort -m 5G -@ 20 ./snps/86k/EW/lane04/$file | ./software/samtools-1.6/samtools markdup -rs -@ 20 - ./snps/86k/EW/lane04/sort.rmdup/$file; done < ./rawreads/EWfasta/lane04/lane04.list; while read file; do ./software/samtools-1.6/samtools sort -m 5G -@ 20 ./snps/86k/EW/lane05/$file | ./software/samtools-1.6/samtools markdup -rs -@ 20 - ./snps/86k/EW/lane05/sort.rmdup/$file; done < ./rawreads/EWfasta/lane05/lane05.list; while read file; do ./software/samtools-1.6/samtools sort -m 5G -@ 20 ./snps/86k/EW/lane06/$file | ./software/samtools-1.6/samtools markdup -rs -@ 20 - ./snps/86k/EW/lane06/sort.rmdup/$file; done < ./rawreads/EWfasta/lane06/lane06.list; while read file; do ./software/samtools-1.6/samtools sort -m 5G -@ 20 ./snps/86k/EW/lane07/$file | ./software/samtools-1.6/samtools markdup -rs -@ 20 - ./snps/86k/EW/lane07/sort.rmdup/$file; done < ./rawreads/EWfasta/lane07/lane07.list; while read file; do ./software/samtools-1.6/samtools sort -m 5G -@ 20 ./snps/86k/EW/lane08/$file | ./software/samtools-1.6/samtools markdup -rs -@ 20 - ./snps/86k/EW/lane08/sort.rmdup/$file; done < ./rawreads/EWfasta/lane08/lane08.list; while read file; do ./software/samtools-1.6/samtools sort -m 5G -@ 20 ./snps/86k/EW/lane09/$file | ./software/samtools-1.6/samtools markdup -rs -@ 20 - ./snps/86k/EW/lane09/sort.rmdup/$file; done < ./rawreads/EWfasta/lane09/lane09.list; while read file; do ./software/samtools-1.6/samtools sort -m 5G -@ 20 ./snps/86k/EW/lane10/$file | ./software/samtools-1.6/samtools markdup -rs -@ 20 - ./snps/86k/EW/lane10/sort.rmdup/$file; done < ./rawreads/EWfasta/lane10/lane10.list; while read file; do ./software/samtools-1.6/samtools sort -m 5G -@ 20 ./snps/86k/EW/lane11/$file | ./software/samtools-1.6/samtools markdup -rs -@ 20 - ./snps/86k/EW/lane11/sort.rmdup/$file; done < ./rawreads/EWfasta/lane11/lane11.list; while read file; do ./software/samtools-1.6/samtools sort -m 5G -@ 20 ./snps/86k/EW/lane12/$file | ./software/samtools-1.6/samtools markdup -rs -@ 20 - ./snps/86k/EW/lane12/sort.rmdup/$file; done < ./rawreads/EWfasta/lane12/lane12.list

while read file; do ./software/samtools-1.6/samtools sort -m 5G -@ 20 ./snps/86k/LGEP/lane01/$file | ./software/samtools-1.6/samtools markdup -rs -@ 20 - ./snps/86k/LGEP/lane01/sort.rmdup/$file; done < ./rawreads/LGEPfasta/lane01/lane01.list

while read file; do ./software/samtools-1.6/samtools sort -m 5G -@ 20 ./snps/86k/LGEP/lane13/$file | ./software/samtools-1.6/samtools markdup -rs -@ 20 - ./snps/86k/LGEP/lane13/sort.rmdup/$file; done < ./rawreads/LGEPfasta/lane13/lane13.list

while read file; do ./software/samtools-1.6/samtools sort -m 5G -@ 20 ./snps/86k/LGEP/lane14/$file | ./software/samtools-1.6/samtools markdup -rs -@ 20 - ./snps/86k/LGEP/lane14/sort.rmdup/$file; done < ./rawreads/LGEPfasta/lane14/lane14.list

while read file; do ./software/samtools-1.6/samtools sort -m 5G -@ 20 ./snps/86k/LGEP/lane15/$file | ./software/samtools-1.6/samtools markdup -rs -@ 20 - ./snps/86k/LGEP/lane15/sort.rmdup/$file; done < ./rawreads/LGEPfasta/lane15/lane15.list

while read file; do ./software/samtools-1.6/samtools sort -m 5G -@ 20 ./snps/86k/LGEP/lane16/$file | ./software/samtools-1.6/samtools markdup -rs -@ 20 - ./snps/86k/LGEP/lane16/sort.rmdup/$file; done < ./rawreads/LGEPfasta/lane16/lane16.list

while read file; do ./software/samtools-1.6/samtools sort -m 5G -@ 20 ./snps/86k/LGEP/lane17/$file | ./software/samtools-1.6/samtools markdup -rs -@ 20 - ./snps/86k/LGEP/lane17/sort.rmdup/$file; done < ./rawreads/LGEPfasta/lane17/lane17.list

while read file; do ./software/samtools-1.6/samtools sort -m 5G -@ 20 ./snps/86k/LGEP/lane18/$file | ./software/samtools-1.6/samtools markdup -rs -@ 20 - ./snps/86k/LGEP/lane18/sort.rmdup/$file; done < ./rawreads/LGEPfasta/lane18/lane18.list

```
***Note: Once finished, unsorted and duplicated bam technical reps were deleted.  The newly created bam files within the ./sort.rmdup directory will be kept as source of techincal reps. 

## 3. Merge technical rep bam files

```linux
 cd /media/disk6/ARF/RNASEQ/snps/86k/EW/

 while read file; do ls -R */*/$file*.fastq.gz >out.txt;/media/disk6/ARF/RNASEQ/software/samtools-1.6/samtools merge -1 -b    out.txt -@ 60 ./$file.ew.bio.rep.merge.bam ;done < test

 cd /media/disk6/ARF/RNASEQ/snps/86k/LGEP/

 while read file; do ls -R */*/Sample_$file*.fastq.gz >out.txt;/media/disk6/ARF/RNASEQ/software/samtools-1.6/samtools merge -1 -b  out.txt -@ 60 ./$file.lgep.bio.rep.merge.bam ;done < test

 /media/disk6/ARF/RNASEQ/software/samtools-1.6/samtools merge -1 -b out.txt -@ 60 ./$file.lgep.bio.rep.merge.bam
```

## 4. Sort biological rep bam file & Remove duplicates

*Create new directory to hold sorted and de-duplicated biological reps & run samtools
 ```linux
 mkdir bio.sort.rmdup
 
 ls *.bam > test

cd ../../..
 
while read file; do ./software/samtools-1.6/samtools sort -m 5G -@ 40 ./snps/86k/EW/$file | ./software/samtools-1.6/samtools markdup -rs -@ 60 - ./snps/86k/EW/bio.sort.rmdup/$file; done < ./snps/86k/EW/test
 ```
*Do the same for LGEP
 ```linux
cd /media/disk6/ARF/RNASEQ/snps/86k/LGEP
 
mkdir bio.sort.rmdup

ls *.bam > test

cd ../../..
 
while read file; do ./software/samtools-1.6/samtools sort -m 5G -@ 60 ./snps/86k/LGEP/$file | ./software/samtools-1.6/samtools markdup -rs -@ 60 - ./snps/86k/LGEP/bio.sort.rmdup/$file; done < ./snps/86k/LGEP/test
 ```
*As before old files were deleted and files were kept in ./bio.sort.rmdup

## 5. Add RG & SM to each biological rep

```linux
cd /media/disk6/ARF/RNASEQ/snps/86k/EW/bio.sort.rmdup

mkdir msd.wrg

ls *.bam > add.rgs

cd /media/disk6/ARF/RNASEQ/

while read file; do java -jar ./software/picard-2.17.4/picard.jar AddOrReplaceReadGroups I=./snps/86k/EW/bio.sort.rmdup/$file O=snps/86k/EW/bio.sort.rmdup/msd.wrg/$file RGID=$file RGLB=$file RGPL=illumina RGPU=$file RGSM=$file ; done < ./snps/86k/EW/bio.sort.rmdup/add.rgs
```

```linux
cd /media/disk6/ARF/RNASEQ/snps/86k/LGEP/bio.sort.rmdup

mkdir msd.wrg

ls *.bam > add.rgs

cd /media/disk6/ARF/RNASEQ/

while read file; do java -jar ./software/picard-2.17.4/picard.jar AddOrReplaceReadGroups I=./snps/86k/LGEP/bio.sort.rmdup/$file O=snps/86k/LGEP/bio.sort.rmdup/msd.wrg/$file RGID=$file RGLB=$file RGPL=illumina RGPU=$file RGSM=$file ; done < ./snps/86k/LGEP/bio.sort.rmdup/add.rgs

```

## 6. Merge biological reps into final bam represnting EW & LGEP resepctively

```linux
cd /media/disk6/ARF/RNASEQ/snps/86k/EW/bio.sort.rmdup/msd.wrg

ls > merge.all.ew

/media/disk6/ARF/RNASEQ/software/samtools-1.6/samtools merge -1 -b ./merge.all.ew -@ 60 EW.snps.merged.bam

/media/disk6/ARF/RNASEQ/software/samtools-1.6/samtools index -@ 60 -b ./EW.snps.merged.bam
```

```linux
cd /media/disk6/ARF/RNASEQ/snps/86k/LGEP/bio.sort.rmdup/msd.wrg

ls > merge.all.lgep

/media/disk6/ARF/RNASEQ/software/samtools-1.6/samtools merge -1 -b ./merge.all.lgep -@ 60 LGEP.snps.merged.bam

```

## 7. Provide full ew and lgep bam to freebayes

Run Freebayes
started at 530pm ct 2/218
```linux
cd /RNASEQ/snps/86k/EW/bio.sort.rmdup/msd.wrg

/home/rosswhet/software/freebayes/scripts/freebayes-parallel <(/home/rosswhet/software/freebayes/scripts/fasta_generate_regions.py /media/disk6/ARF/RNASEQ/references/Pita.86ktxptome.fasta.fai 100000) 60 -f /media/disk6/ARF/RNASEQ/references/Pita.86ktxptome.fasta -b /media/disk6/ARF/RNASEQ/snps/86k/EW/bio.sort.rmdup/msd.wrg/EW.snps.merged.bam | bgzip -c > ew.bam.freebayes.out.vcf.gz
```

## 8. Tabix output vcf files, merge both vcf files, filter

```linux
cd /RNASEQ/snps/86k/EW/bio.sort.rmdup/msd.wrg

tabix ./ew.bam.freebayes.out.vcf.gz

cd /RNASEQ/snps/86k/LGEP/bio.sort.rmdup/msd.wrg

taxib ./lgep.bam.freebayes.out.vcf.gz
```

```linux
 bcftools merge --threads 60 ./EW/bio.sort.rmdup/msd.wrg/ew.bam.freebayes.out.vcf.gz ./LGEP/bio.sort.rmdup/msd.wrg/lgep.bam.freebayes.out.vcf.gz  | bgzip -c > ./192.biorep.vcf.bgzip.gz
 
 tabix 192.biorep.vcf.gz
 ```
 
``` linux
vcftools --gzvcf ./192.biorep.vcf.bgzip.gz --remove-indels --max-missing 1.0 --maf 0.05 --min-meanDP 10  --012 
```

```linux
vcftools --gzvcf ./192.biorep.vcf.bgzip.gz --remove-indels --max-missing 1.0 --maf 0.05 --min-meanDP 10  --relatedness2 
```
