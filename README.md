# DNA-BSAM
This is essentially a wrapper script that helps to automate some the Qiime2 steps. The pipeline works as follows:

__1.__ Create directories and a proxy metadata file including just filenames and an arbitary metadata column.

__2.__ Import the data into a Qiime2.

__3.__ Use cutadapt to trim primers. Select primers from amplicons currently supported in the primer.tsv file, or add your own

__4.__ Truncate reads (requires manual input) and denoise with Dada2.

__5.__ Filter again for chimeras with Uchime.

__6.__ Generate features.

__7.__ Filter based on the number of reads.

__8.__ Classify features with a user specified feature classifier.

__9.__ Filter (either include or exclude) based on taxa specific to each amplicon. Also generate a tabulated data output to the Metadata folder where feature ID, sequence, classification and confidence are recorded.

__10.__ Generate bar plots.

__11.__ Generate Shannon diversity indicies.

__12.__ Calculate Cmin for SRS rarefaction.

__13.__ Generate faith pd, pileou evenness and jaccard diversity outputs.

__Installation__

This pipeline should be run in a Qiime2 conda environment. Currently tested environments include Qiime2-2021.2.

Create a directory name Qiime2. Inside this, create a directory called Scripts. Copy the Qiime2.py script into this directory.

Getting started Create a directory for your project in the Qiime2 directory. Inside this new directory, create a raw_data directory and place your raw data there.

Classifiers Classifiers can be trained for Qiime2 using the instructions located here https://docs.qiime2.org/2020.8/tutorials/feature-classifier/.

Example usage:
```
python Scripts/Qiime2.py --input 16S_project/raw_data/ --output 16S_project/ --q2_classifier silva-138-99-515-806-nb-classifier.qza --amplicon 16S --cutadapt_times 2 --filter_m_ie exclude --filter_list mitochondria,chloroplast,archaea --classify_threads 20 --classify_conf 0.7
```
