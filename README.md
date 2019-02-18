# Collection of RNA-SEQ Data

- [Collection-of-RNA-SEQ-Data](#collection-of-rna-seq-data)
  * [Abstract](#abstract)
  * [Background of samples](#background-of-samples)
  * [Location of data](#location-of-data)
  * [Analyses](#analyses)
    + [Step 1 - Data Prep](#step-1---data-prep)
    + [Step 2 - Load Count Data](#step-2---load-count-data)

## Abstract

* Two batches of loblolly pine families were grown in Fall 2015 and Spring 2016 in order to test 

## Background of samples

* LGEP
   
   - 144 Biological Replicates x 3 technical replicates = 432 technical replicates

* EW

   - 80 Biological Replicates x 3 technical replicates = 240 technical replicates

## Location of data

Data Subject Type | Data File Type | Path | Notes
--- | --- | --- | ---
**raw read files**  | | `/media/disk6/ARF/RNASEQ/shared/rawreads/86kSalmon` | Raw files returned from GSL
| | *raw tar* | `./EWtarfiles or ./LGEPtarfiles`
| | *raw fasta* | `./EWfasta or ./LGEPfasta` 
**trimmed and filtered read files** | |`/media/disk6/ARF/RNASEQ/shared/trimmedfiltreads/86k` | Files post trim & adapater removal
|  |*EW* | `./EW/lane01 ... ./lane12` | 
|  |*LGEP* | `./LGEP/lane01 ... ./lane18` | 
**salmon count files** | |`/media/disk6/ARF/RNASEQ/shared/counts/86kSalmon` | Direcotries containing quant.sf files
|  |*EW tech reps* | `./EW/lane01 ... ./lane12` | 
|  |*LGEP tech reps* | `./LGEP/lane01 ... ./lane18` |
|  |*EW bio reps* | `./bio_EW/Sample_<animal_id>/` | 
|  |*LGEP bio reps* | `./bio_LGEP/Sample_<animal_id>/` | 
**experimental data resources**  | | `/media/disk6/ARF/RNASEQ/BV-Prediction/Breeding-Value-Prediction/disk6directory/resources` | Experiment information
|  |*sequencing* | `./exptdesign/sequencing` | `./EWtarfiles or ./LGEPtarfiles`
|  |*pedigree* | `./pedigree` | `./EWfasta or ./LGEPfasta`
|  |*phenotypes* | `./phenos` | `./EWfasta or ./LGEPfasta`


## Analyses

### Step 1 - Data Prep

   Data prep includes everything from unpacking the original tar files recieved by the GSL up to estimating transcript               
      abundance with Salmon. Additionally, this step includes identification of the indicies used within both batches and creates an experimental info matrix containing all meta data from both batches.
      
   See the [raw reads README](http://htmlpreview.github.com/?https://github.com/arfesta/RNA_SEQ_Data/blob/master/disk6directory/rawreads/012718raw_data_processing.html) for step by step processing of files.  
   
  Extra prep scripts: [sample index identification](http://htmlpreview.github.com/?https://github.com/arfesta/RNA_SEQ_Data/blob/master/disk6directory/analyses/step1.data_prep/identify_index_used.ouput.html) & [creation of experimental data](http://htmlpreview.github.com/?https://github.com/arfesta/RNA_SEQ_Data/blob/master/disk6directory/analyses/step1.data_prep/create_expt_data.html)

### Step 2 - Load Count Data
      
   Once counts have been estimated, the next step involves reading in the aligned technical, or biological, replicate counts using the tximport package.
      
   Additionally, the phenotype and other sample meta-data is constructed for normalization.
   
   To see this process for **biological reps**, navigate to: 
   
   [load counts bio rep html file](http://htmlpreview.github.com/?https://github.com/arfesta/RNA_SEQ_Data/blob/master/disk6directory/analyses/step2.loadcounts/load.counts.html) which contains the complete markdown and output.
   
  Origin information was added to the biological replicate matrix [Script](https://github.com/arfesta/RNA_SEQ_Data/blob/master/disk6directory/analyses/step2.load.counts.add_origin_dat_step2.R)

   To see this process for **technical reps**, navigate to: 
   
   [load counts tech rep html file](http://htmlpreview.github.com/?https://github.com/arfesta/RNA_SEQ_Data/blob/master/disk6directory/analyses/step2.loadcounts/load.counts_techreps.html) which contains the complete markdown and output.

