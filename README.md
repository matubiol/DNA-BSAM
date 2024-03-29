# DNA-BSAM
This is essentially a wrapper script that helps to automate some of the QIIME2 steps. The pipeline works as follows:

__1.__ Create directories and a proxy metadata file including just filenames and an arbitary metadata column.

__2.__ Import the data into a QIIME2.

__3.__ Use cutadapt to trim primers. Select primers from amplicons currently supported in the primer.tsv file, or add your own

__4.__ Truncate reads (requires manual input) and denoise with DADA2.

__5.__ Filter again for chimeras with UCHIME.

__6.__ Generate features.

__7.__ Cluster features.

__8.__ Filter based on the number of reads.

__9.__ Classify features with a user specified feature classifier.

__10.__ Filter (either include or exclude) based on taxa specific to each amplicon. Also generate a tabulated data output to the Metadata folder where feature ID, sequence, classification and confidence are recorded.

__11.__ Generate bar plots.

__12.__ Generate Shannon diversity outcomes.

__13.__ Calculate Cmin for SRS normalization.

__14.__ Calculate taxa richness.

__15.__ Constione with DNA-BSAM-R.R

__Installation__

This pipeline should be run in a QIIME2 Conda environment. Currently tested environments include QIIME2-2021.2.

Create a directory name QIIME2. Inside this, create a directory called Scripts. Copy the DNA-BSAM.py script into this directory.

__Getting started__

Create a directory for your project in the QIIME2 directory. Inside this new directory, create a raw_data directory and place your raw data there.

Classifiers Classifiers can be trained for QIIME2 using the instructions located here https://docs.qiime2.org/2023.2/tutorials/feature-classifier/

Example usage:
```
python DNA-BSAM-QIIME2.py --input 16S_project/raw_data/ --output 16S_project/ --q2_classifier silva-138-99-515-806-nb-classifier.qza --amplicon 16S --cutadapt_times 2 --filter_m_ie exclude --filter_list mitochondria,chloroplast,archaea --classify_threads 20 --classify_conf 0.7
```

Parameters:

```
usage: DNA-BSAM-QIIME2.py [-h] --input INPUT --output OUTPUT
                 [--q2_classifier Q2_CLASSIFIER]
                 [--amplicon {16S,ITS,18S,APTA,BF3BR2,mlCOIintF-jgHCO2198,rbcla,trnL,ITS6-5.8S-1R}]
                 [--primer_file PRIMER_FILE] [--create_dirs {Y,N}]
                 [--metadata {Y,N}] [--import_q2 {Y,N}] [--cutadapt {Y,N}]
                 [--dada2 {Y,N}] [--chimera {Y,N}] [--features {Y,N}] [--cluster {Y,N}]
                 [--filter {Y,N}] [--classify {Y,N}] [--barplots {Y,N}]
                 [--shannon {Y,N}] [--richness {Y,N}] [--srs {Y,N}]
                 [--cutadapt_times CUTADAPT_TIMES]
                 [--cutadapt_error CUTADAPT_ERROR] [--filter_n FILTER_N]
                 [--filter_list FILTER_LIST]
                 [--filter_m_ie {exclude,include,none}]
                 [--dada2_EE_rev DADA2_EE_REV] [--dada2_minfold DADA2_MINFOLD]
                 [--classifier {nb,vsearch}]
                 [--classify_threads CLASSIFY_THREADS]
                 [--classify_conf CLASSIFY_CONF] [--vsearch_db VSEARCH_DB]
                 [--vsearch_taxonomy VSEARCH_TAXONOMY]

Runs the Fera Qiime2 metabarcoding pipeline.

optional arguments:
  -h, --help            show this help message and exit
  --input INPUT         This is the location of the raw data directory.
  --output OUTPUT       This where the output data will be generated.
  --q2_classifier Q2_CLASSIFIER
                        This is the path to the classifier.
  --amplicon {16S,ITS,18S,APTA,BF3BR2,mlCOIintF-jgHCO2198,rbcla,trnL,ITS6-5.8S-1R}
                        This is a reference to the primers used for the marker
                        region. Please update the dictionary in the script for
                        more choices.
  --primer_file PRIMER_FILE
                        Tab separated list of the amplicon and the primers
                        used. Follow the example in the example provided.                       
  --create_dirs {Y,N}   Creates the directory structure. Default Y.
  --metadata {Y,N}      Create metadata file. Default Y.
  --import_q2 {Y,N}     Import data into QIIME2. Default Y.
  --cutadapt {Y,N}      Run Cutadapt. Default Y.
  --dada2 {Y,N}         Run DADA2. Default Y.
  --chimera {Y,N}       Run UCHIME for chimera removal. Default Y.
  --features {Y,N}      Generate features from chimera removed data. Default
                        Y.
  --cluster {Y,N}       Generate clustered features. Default Y.
  --filter {Y,N}        Filter samples. Default Y.
  --classify {Y,N}      Classify features. Default Y.
  --barplots {Y,N}      Generate barplots. Default Y.
  --shannon {Y,N}       Generate Shannon diversity outcomes. Default Y.
  --srs {Y,N}           Run SRS. Default Y.
  --richness {Y,N}      Calculate taxa richness. Default Y.
  --cutadapt_times CUTADAPT_TIMES
                        The number of times to look for primers to trim in the
                        reads. Default 2.
  --cutadapt_error CUTADAPT_ERROR
                        The number of errors to allow when detecting primers
                        to trim. Default 0.2.
  --filter_n FILTER_N   Remove sample IDs which have fewer than X reads.
                        Default 3000.
  --filter_list FILTER_LIST
                        Remove sample features based on comma separated list.
                        Default mitochondria,chloroplast,archaea,eukaryota.
  --filter_m_ie {exclude,include,none}
                        Include or exclude sample features.
  --dada2_EE_rev DADA2_EE_REV
                        The maximum number of expected errors in the reverse
                        read. Default 2.
  --dada2_minfold DADA2_MINFOLD
                        The minimum abundance of potential parents of a
                        sequence being tested as chimeric, expressed as a
                        fold-change versus the abundance of the sequence being
                        tested . Default 1.
  --classifier {nb,vsearch}
                        The classifier to run. Default nb classifier.
  --classify_threads CLASSIFY_THREADS
                        The number of threads to run the classification tool
                        with. Default 60.
  --classify_conf CLASSIFY_CONF
                        The confidence required for a successful
                        classification. Default 0.7.
  --vsearch_db VSEARCH_DB
                        The path to the database sequence file.
  --vsearch_taxonomy VSEARCH_TAXONOMY
                        The path to the database taxonomy file.
```
