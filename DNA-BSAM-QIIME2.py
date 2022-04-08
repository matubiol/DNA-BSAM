import os
import sys
import subprocess
import multiprocessing
import argparse
import re
import random
import string
from matplotlib import cbook

def main():

	# Qiime2 script version
	qiime2_script_version = 1.2

	### Input arguments
	options = parseArguments()

### Primers

	# Dictionary containing  primers as keys
	# Each key has a list attached with 4 elements
	# forward, reverse, forward reverse complement, reverse reverse complement
	# f_primer, r_primer, f_primer_rc, r_primer_rc 
	primer_dict = {}
	with open(options.primer_file) as primers:
		next(primers)
		for line in primers:
			line = line.strip()
			primer_dict[line.split("\t")[0]] = line.split("\t")[1:]

### Create project folder

	if(os.path.isdir(options.output)):
		print("Project folder already created.")
	else:
		print("Creating project folder.")
		os.mkdir(options.output)

### Intial directory setup

	if(options.create_dirs == "Y"):
		os.mkdir(f"{options.output}/Metadata")
		os.mkdir(f"{options.output}/Import")
		os.mkdir(f"{options.output}/Cutadapt")
		os.mkdir(f"{options.output}/Dada2")
		os.mkdir(f"{options.output}/Filtered_data")
		os.mkdir(f"{options.output}/Features")
		os.mkdir(f"{options.output}/Classify")
		os.mkdir(f"{options.output}/Barplots")
		os.mkdir(f"{options.output}/Shannon")
		os.mkdir(f"{options.output}/SRS")
		os.mkdir(f"{options.output}/Tree")
		os.mkdir(f"{options.output}/Diversity")

	else:
		print("Directory structure creation skipped.")

### Metadata setup

	if(options.metadata == "Y"):
		file_list = []
		raw_data_list = os.listdir(options.input)
		for file in raw_data_list:
			if("_R1_" in file):
				filesplit = file.split("_S")
				file_list.append(filesplit[0])

		with open(f"{options.output}/Metadata/metadata.tsv", "w") as metadata_out:
			metadata_out.write("#SampleID" + "	" + "Proxy" + "\n")
			for ID in file_list:
				metadata_out.write(ID + "	" + random.choice("ABCXYZ") + "\n")


### Import data

	if(options.import_q2 == "Y"):
		input_import = options.input
		output_import = f"{options.output}/Import/{options.amplicon}"
		print("Importing data into Qiime2.")
		
		# Import command
		subprocess.call(f"qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path {input_import} --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path {output_import}", shell = True)
		
		# Visualise imported data
		subprocess.call(f"qiime demux summarize --i-data {output_import}.qza --o-visualization {output_import}", shell = True)
		print(f"Importing into Qiime2 complete for {options.output}")

	else:
		print("Importing data into Qiime2 format skipped.")

### Run Cutadapt

	if(options.cutadapt == "Y"):
		input_cutadapt = f"{options.output}/Import/{options.amplicon}.qza"
		output_cutadapt = f"{options.output}/Cutadapt/{options.amplicon}"
		print(f"Running cutadapt with {options.amplicon} primers selected.")

		# Cutadapt command
		subprocess.call(f"qiime cutadapt trim-paired --i-demultiplexed-sequences {input_cutadapt} " + 
			f"--p-front-f {primer_dict[options.amplicon][0]} --p-adapter-f {primer_dict[options.amplicon][3]} --p-front-r {primer_dict[options.amplicon][1]} --p-adapter-r {primer_dict[options.amplicon][2]} " +
			f"--p-cores 60 --p-times {options.cutadapt_times} --p-error-rate {options.cutadapt_error} --p-minimum-length 50 --o-trimmed-sequences {output_cutadapt} --verbose > {output_cutadapt}.log", shell = True)

		# Visualise cudadapt results
		subprocess.call(f"qiime demux summarize --i-data {output_cutadapt}.qza --o-visualization {output_cutadapt}", shell=True)
		print(f"Cutadapt complete for {options.output}")

	else:
		print("Cutadapt skipped.")

### Run Dada2
	
	trunclen_f = ""
	trunclen_r = ""
	if(options.dada2 == "Y"):
		input_dada2 = f"{options.output}/Cutadapt/{options.amplicon}.qza"
		output_dada2_table = f"{options.output}/Dada2/{options.amplicon}-table"
		output_dada2_repseq = f"{options.output}/Dada2/{options.amplicon}-repseq"
		output_dada2_stats = f"{options.output}/Dada2/{options.amplicon}-stats"
		print("please inspect the qzv file using qiime2 - view https://view.qiime2.org/ - in the Cutadapt folder and select appropriate cutoffs.")
		trunclen_f = input("Please input length to truncate forward reads to: ")
		trunclen_r = input("Please input length to truncate reverse reads to: ")
		print("Running dada2.")

		# Dada2 command
		subprocess.call(f"qiime dada2 denoise-paired --i-demultiplexed-seqs {input_dada2} --p-trunc-len-f {trunclen_f} --p-trunc-len-r {trunclen_r} --p-n-threads 60 " + 
			f"--o-table {output_dada2_table} --p-max-ee-r {options.dada2_EE_rev} --o-representative-sequences {output_dada2_repseq} --p-min-fold-parent-over-abundance {options.dada2_minfold} --o-denoising-stats {output_dada2_stats}", shell = True)

		# Visualise Dada2 result
		subprocess.call(f"qiime metadata tabulate --m-input-file {output_dada2_stats}.qza --o-visualization {output_dada2_stats}", shell=True)
		print(f"Dada2 complete for {options.output}")

	else:
		print("Dada2 skipped.")

### Chimera removal

	if(options.chimera == "Y"):
		input_chimera_repseq = f"{options.output}/Dada2/{options.amplicon}-repseq.qza"
		input_chimera_table = f"{options.output}/Dada2/{options.amplicon}-table.qza"
		output_chiemra_dir = f"{options.output}/Chimera_removal"
		print("Removing chimeras with Uchime.")

		# Uchime command
		subprocess.call(f"qiime vsearch uchime-denovo --i-sequences {input_chimera_repseq} --i-table {input_chimera_table} --output-dir {output_chiemra_dir} --verbose", shell = True)
		subprocess.call(f"qiime metadata tabulate --m-input-file {output_chiemra_dir}/stats.qza --o-visualization {output_chiemra_dir}/stats.qzv", shell = True)

		print(f"Chimeral removal complete for {options.output}")

	else:
		print("Chimera removal skipped.")

### Generate features
	
	if(options.features == "Y"):
		input_features_repseq = f"{options.output}/Dada2/{options.amplicon}-repseq.qza"
		input_features_table = f"{options.output}/Dada2/{options.amplicon}-table.qza"
		input_chimeras = f"{options.output}/Chimera_removal/nonchimeras.qza"
		output_features_repseq = f"{options.output}/Features/{options.amplicon}-repseq"
		output_features_table = f"{options.output}/Features/{options.amplicon}-table"

		print("Generating features.")

		# Feature generation command
		subprocess.call(f"qiime feature-table filter-seqs --i-data {input_features_repseq} --m-metadata-file {input_chimeras} --o-filtered-data {output_features_repseq}.qza", shell = True)
		subprocess.call(f"qiime feature-table filter-features --i-table {input_features_table} --m-metadata-file {input_chimeras}  --o-filtered-table {output_features_table}.qza", shell = True)
		subprocess.call(f"qiime feature-table tabulate-seqs --i-data {output_features_repseq}.qza --o-visualization {output_features_repseq}.qzv", shell = True)
		subprocess.call(f"qiime feature-table summarize --i-table {output_features_table}.qza --o-visualization {output_features_table}.qzv", shell = True)

		print(f"Feature tables created for {options.output}")

	else:
		print("Feature table creation skipped.")

### Filter samples

	if(options.filter == "Y"):
		input_features_dir = f"{options.output}/Features/"
		output_features = f"{options.output}/Filtered_data/ID_filtered_table"
		print("Filtering sample IDs based on reads.")

		# Unzip archive and extract read information
		subprocess.call(f"unzip -p {input_features_dir}{options.amplicon}-table.qzv */data/sample-frequency-detail.csv > {input_features_dir}reads.csv", shell = True)

		# Get samples with fewer than X reads
		sample_list = []
		with open(f"{input_features_dir}reads.csv") as read_count:
			next(read_count)
			for line in read_count:
				line_list = line.split(",")
				if(float(line_list[1]) < int(options.filter_n)):
					print(f"Sample: {line_list[0]} has fewer than {options.filter_n} reads so will be removed. Number of reads = {line_list[1]}")
					sample_list.append(line_list[0])

		# Write samples2remove file
		with open(f"{options.output}/samples2remove.txt", "w") as samples2remove:
			samples2remove.write("#SampleID" + "\n")
			for sample in sample_list:
				samples2remove.write(sample + "\n")

		# Filter features if any samples removed, if not then copy the table to the Filtered data directory
		if not sample_list:
			print("No samples removed during filtering.")
			subprocess.call(f"cp {input_features_dir}{options.amplicon}-table.qza {output_features}.qza", shell = True)
		else:
			subprocess.call(f"qiime feature-table filter-samples --i-table {input_features_dir}{options.amplicon}-table.qza --m-metadata-file {options.output}/samples2remove.txt --p-exclude-ids --o-filtered-table {output_features}.qza", shell = True)

		subprocess.call(f"qiime feature-table summarize --i-table {output_features}.qza --o-visualization {output_features}.qzv", shell = True)

		print(f"Samples filtered {options.output}")

	else:
		print("Sample filtering skipped.")

### Classify features

	if(options.classify == "Y"):
		input_classifier = f"{options.output}/Features/{options.amplicon}-repseq.qza"
		output_classifier = f"{options.output}/Classify/{options.amplicon}-taxonomy"
		print("Classifying samples.")

		# Classification command
		if(options.classifier == "nb"):
			subprocess.call(f"qiime feature-classifier classify-sklearn --i-classifier {options.q2_classifier} --i-reads {input_classifier} --o-classification {output_classifier} --p-confidence {options.classify_conf} --p-n-jobs {options.classify_threads}", shell = True)
		elif(options.classifier == "vsearch"):
			subprocess.call(f"qiime feature-classifier classify-consensus-vsearch --i-query {input_classifier} --i-reference-reads {options.vsearch_db} --i-reference-taxonomy {options.vsearch_taxonomy} --o-classification {output_classifier} --p-threads {options.classify_threads}", shell = True)
		subprocess.call(f"qiime metadata tabulate --m-input-file {output_classifier}.qza --o-visualization {output_classifier}.qzv", shell = True)

		print(f"Samples classified {options.output}")

	else:
		print("Sample classification skipped.")

### Filter samples - taxonomy

	if(options.filter == "Y"):

		input_table = f"{options.output}/Filtered_data/ID_filtered_table.qza"
		input_taxonomy = f"{options.output}/Classify/{options.amplicon}-taxonomy.qza"
		output_table = f"{options.output}/Filtered_data/ID_taxonomy_filtered_table"
		print("Filtering samples based on taxonomy")

		# Filter samples for taxonomy command
		# If none is selected, then copy exisiting table to this location
		if(options.filter_m_ie != "none"):
			subprocess.call(f"qiime taxa filter-table --i-table {input_table} --i-taxonomy {input_taxonomy} --p-{options.filter_m_ie} {options.filter_list} --o-filtered-table {output_table}.qza", shell = True)
		else:
			subprocess.call(f"cp {input_table} {output_table}.qza", shell = True)
			print("No taxa filtered.")

		subprocess.call(f"qiime feature-table summarize --i-table {output_table}.qza --o-visualization {output_table}.qzv", shell = True)

		print(f"Samples filtered by metadata for {options.output}")

		# Generate tabulated data
		subprocess.call(f"qiime metadata tabulate --m-input-file {options.output}/Features/{options.amplicon}-repseq.qza --m-input-file {options.output}/Classify/{options.amplicon}-taxonomy.qza --o-visualization {options.output}/Metadata/{options.amplicon}-tabulated_data.qzv", shell = True)

	else:
		print("Taxonomy filtering skipped.")

### Generate Barplots

	if(options.barplots == "Y"):
		input_table = f"{options.output}/Filtered_data/ID_taxonomy_filtered_table.qza"
		input_taxonomy = f"{options.output}/Classify/{options.amplicon}-taxonomy.qza"
		input_metadata = f"{options.output}/Metadata/metadata.tsv"
		output_barplot = f"{options.output}/Barplots/{options.amplicon}-barplot.qzv"
		print("Generating barplots.")

		# Barplots commands
		subprocess.call(f"qiime taxa barplot --i-table {input_table} --i-taxonomy {input_taxonomy} --m-metadata-file {input_metadata} --o-visualization {output_barplot}", shell = True)

		print(f"Bar plots generated for {options.output}")

	else:
		print("Barplot generation skipped.")

### Generate Shannon Alpha Diversity

	if(options.shannon == "Y"):
		input_table = f"{options.output}/Filtered_data/ID_taxonomy_filtered_table.qza"
		output_shannon = f"{options.output}/Shannon/{options.amplicon}-shannon.qza"
		print("Generating Shannon diversity")

		# Shannon alpha diversity command
		subprocess.call(f"qiime diversity alpha --i-table {input_table} --p-metric shannon --o-alpha-diversity {output_shannon}", shell = True)

		print(f"Shannon diversity generated for {options.output}")

	else:
		print("Shannon diversity generation skipped.")

### SRS 
	whislo = ""
	if(options.srs == "Y"):
		input_table = f"{options.output}/Filtered_data/ID_taxonomy_filtered_table"
		output_features = f"{options.output}/SRS/features.csv"
		output_srs_table = f"{options.output}/SRS/SRS_table"
		print("Performing SRS subsampling")

		# Get Cmin value
		# Unzip archive and extract read information
		subprocess.call(f"unzip -p {input_table}.qzv */data/sample-frequency-detail.csv > {output_features}", shell = True)

		# Get samples with feature counts
		feature_count_list = []
		with open(output_features) as feature_count:
			next(feature_count)
			for line in feature_count:
				line = line.strip()
				line_list = line.split(",")
				feature_count_list.append(float(line_list[1]))

		# Get the lower whisker to use as the cutoff
		whislo = cbook.boxplot_stats(feature_count_list)
		whislo = int(whislo[0]["whislo"])

		# Run SRS
		subprocess.call(f"qiime srs SRS --i-table {input_table}.qza --p-c-min {whislo} --o-normalized-table {output_srs_table}.qza --verbose", shell = True)
		subprocess.call(f"qiime feature-table summarize --i-table {output_srs_table}.qza --o-visualization {output_srs_table}.qzv", shell = True)

		print(f"SRS table generated for {options.output}")

	else:
		print("SRS table generation skipped.")

### Diversity

	if(options.diversity == "Y"):
		
		# Tree variables
		input_repseq = f"{options.output}/Features/{options.amplicon}-repseq.qza"
		output_repseq_aln = f"{options.output}/Tree/{options.amplicon}-aligned-repseq.qza"
		output_repseq_mask_aln = f"{options.output}/Tree/{options.amplicon}-masked-aligned-repseq.qza"
		output_rooted_tree = f"{options.output}/Tree/rooted_tree.qza"
		output_unrooted_tree = f"{options.output}/Tree/unrooted_tree.qza"
		
		# Diversity variables
		input_table = f"{options.output}/SRS/SRS_table.qza"
		metadata = f"{options.output}/Metadata/metadata.tsv"
		output_faith = f"{options.output}/Diversity/faith_pd_vector"
		output_evenness = f"{options.output}/Diversity/evenness_vector"
		output_features = f"{options.output}/Diversity/observed_features"
		output_jaccard = f"{options.output}/Diversity/jaccard"
		print("Generating phylogenetic trees and running beta diversity metrics")

		# Tree
		subprocess.call(f"qiime phylogeny align-to-tree-mafft-fasttree --i-sequences {input_repseq} --o-alignment {output_repseq_aln} --o-masked-alignment {output_repseq_mask_aln} --o-tree {output_unrooted_tree} --o-rooted-tree {output_rooted_tree}", shell = True)

		# Diversity
		# Faith
		subprocess.call(f"qiime diversity-lib faith-pd --i-table {input_table} --i-phylogeny {output_rooted_tree} --o-vector {output_faith}.qza", shell = True)
		subprocess.call(f"qiime diversity alpha-group-significance --i-alpha-diversity {output_faith}.qza --m-metadata-file {metadata} --o-visualization {output_faith}-group-significance.qzv", shell = True)

		# Pielou evenness
		subprocess.call(f"qiime diversity-lib pielou-evenness --i-table {input_table} --o-vector {output_evenness}.qza", shell = True)
		subprocess.call(f"qiime diversity alpha-group-significance --i-alpha-diversity {output_evenness}.qza --m-metadata-file {metadata} --o-visualization {output_evenness}-group-significance.qzv", shell = True)

		# Observed features
		subprocess.call(f"qiime diversity alpha --i-table {input_table} --p-metric 'observed_features' --o-alpha-diversity {output_features}.qza", shell = True)
		subprocess.call(f"qiime diversity alpha-group-significance --i-alpha-diversity {output_features}.qza --m-metadata-file {options.output}/Metadata/metadata.tsv --o-visualization {output_features}.qzv", shell = True)

		# Jaccard
		subprocess.call(f"qiime diversity-lib jaccard --i-table {input_table} --p-n-jobs auto --o-distance-matrix {output_jaccard}.qza", shell = True)
		subprocess.call(f"qiime diversity pcoa --i-distance-matrix {output_jaccard}.qza --o-pcoa {output_jaccard}-pcoa.qza", shell = True)
		subprocess.call(f"qiime emperor plot --i-pcoa {output_jaccard}-pcoa.qza --m-metadata-file {metadata} --o-visualization {output_jaccard}.qzv", shell = True)

		print(f"Beta diversity generated for {options.output}")

	else:
		print("Beta diversity generation skipped.")

### Program Details

	print("Printing Qiime2 pipeline version information")
	with open(f"{options.output}/Qiime2Pipeline_params.txt", "w") as qiime2_log:
		qiime2_log.write(f"Qiime2 Pipeline Version: {qiime2_script_version}\n")
		qiime2_log.write(f"Classifier:  {options.classifier}\n")
		qiime2_log.write(f"Vsearch classifier database:  {options.vsearch_db}\n")
		qiime2_log.write(f"cVsearch classifier taxonomy:  {options.vsearch_taxonomy}\n")
		qiime2_log.write(f"NB classifier file:  {options.q2_classifier}\n")
		qiime2_log.write(f"Amplicon: {options.amplicon}\n")
		qiime2_log.write(f"Dada2 truncation lengths: forward - {trunclen_f} reverse - {trunclen_r}\n")
		qiime2_log.write(f"Dada2 expected errors: forward - 2 reverse - {options.dada2_EE_rev}\n")
		qiime2_log.write(f"Number of reads required for a sample to pass filter: {options.filter_n}\n")
		qiime2_log.write(f"Amplicon taxonomy filtered for: {options.filter_m_ie} - {options.filter_list}\n")
		qiime2_log.write(f"Classification confidence: {options.classify_conf}\n")
		qiime2_log.write(f"SRS Cmin: {whislo}\n")
		subprocess.call(f"conda list > {options.output}/Qiime2Pipeline_env.txt", shell = True)

	print("Qiime2 pipeline completed!")


################################################################################################################################ Functions ################################################################################################################################

#################################### Get Arguments Function ####################################
def parseArguments():
	parser = argparse.ArgumentParser(description = "Runs the Fera Qiime2 metabarcoding pipeline.")

	# Main arguments
	parser.add_argument("--input", help = "This is the location of the raw data directory.", required = True)
	parser.add_argument("--output", help = "This where the output data will be generated.", required = True)

	# Main arguments
	parser.add_argument("--q2_classifier", help = "This is the path to the classifier.", default = "/biostore/bigbio_00/databases/classifiers/Silva_138_16S_qiime2.20.6/silva-138-99-515-806-nb-classifier.qza")
	parser.add_argument("--amplicon", help = "This is a reference to the primers used for the marker region. Please update the dictionary in the script for more choices.", choices = ["16S", "ITS", "18S", "APTA", "BF3BR2", "mlCOIintF-jgHCO2198", "rbcla", "trnL", "ITS6-5.8S-1R"])
	parser.add_argument("--primer_file", help = "Tab separated list of the amplicon and the primers used. Follow the example in the example provided. Default /home/smcgreig/Scripts/Qiime2/primers.tsv", default = "/home/smcgreig/Scripts/Qiime2/primers.tsv")

	# Extra arguments, useful for if a specific job has failed and you don't want to start from scratch
	parser.add_argument("--create_dirs", help = "Creates the directory structure. Default Y.", choices = ["Y", "N"], default = "Y")
	parser.add_argument("--metadata", help = "Create metadata file. Default Y.", choices = ["Y", "N"], default = "Y")
	parser.add_argument("--import_q2", help = "Import data into qiime2. Default Y.", choices = ["Y", "N"], default = "Y")
	parser.add_argument("--cutadapt", help = "Run cutadapt. Default Y.", choices = ["Y", "N"], default = "Y")
	parser.add_argument("--dada2", help = "Run dada2. Default Y.", choices = ["Y", "N"], default = "Y")
	parser.add_argument("--chimera", help = "Run uchime for chimera removal. Default Y.", choices = ["Y", "N"], default = "Y")	
	parser.add_argument("--features", help = "Generate features from chimera removed data. Default Y.", choices = ["Y", "N"], default = "Y")
	parser.add_argument("--filter", help = "Filter samples. Default Y.", choices = ["Y", "N"], default = "Y")
	parser.add_argument("--classify", help = "Classify features. Default Y.", choices = ["Y", "N"], default = "Y")
	parser.add_argument("--barplots", help = "Generate barplots. Default Y.", choices = ["Y", "N"], default = "Y")
	parser.add_argument("--shannon", help = "Generate shannon alpha diversity stats. Default Y.", choices = ["Y", "N"], default = "Y")
	parser.add_argument("--srs", help = "Run SRS. Default Y.", choices = ["Y", "N"], default = "Y")
	parser.add_argument("--diversity", help = "Run diversity metrics. Default Y.", choices = ["Y", "N"], default = "Y")

	# Cutadapt arguments 
	parser.add_argument("--cutadapt_times", help = "The number of times to look for primers to trim in the reads. Default 2.", default = "2")
	parser.add_argument("--cutadapt_error", help = "The number of errors to allow when detecting primers to trim. Default 0.2.", default = "0.2")

	# Filter arguments
	parser.add_argument("--filter_n", help = "Remove sample IDs which have fewer than X reads. Default 3000.", default = "3000")
	parser.add_argument("--filter_list", help = "Remove sample features based on comma separated list. Default mitochondria,chloroplast,archaea,eukaryota.", default = "mitochondria,chloroplast,archaea,eukaryota")
	parser.add_argument("--filter_m_ie", help = "Include or exclude sample features.", choices = ["exclude", "include", "none"], default = "exclude")

	# Dada2 arguments
	parser.add_argument("--dada2_EE_rev", help = "The maximum number of expected errors in the reverse read. Default 2.", default = "2")
	parser.add_argument("--dada2_minfold", help = "The minimum abundance of potential parents of a sequence being tested as chimeric, expressed as a fold-change versus the abundance of the sequence being tested . Default 1.", default = "1")

	# General classifier arguments
	parser.add_argument("--classifier", help = "The classifier to run. Default nb classifier.", default = "nb", choices = ["nb", "vsearch"])
	parser.add_argument("--classify_threads", help = "The number of threads to run the classification tool with. Default 60.", default = 60)

	# NB classifier arguments
	parser.add_argument("--classify_conf", help = "The confidence required for a successful classification. Default 0.7.", default = "0.7")

	# Vsearch classification arguments
	parser.add_argument("--vsearch_db", help = "The path to the database sequence file. Default /data/bigbio_00/smcgreig/18S_nematode/18S_nematode_full_ncbi/18S_nhmmer_final_seqs_uniq.qza", default = "/data/bigbio_00/smcgreig/18S_nematode/18S_nematode_full_ncbi/18S_nhmmer_final_seqs_uniq.qza")
	parser.add_argument("--vsearch_taxonomy", help = "The path to the database taxonomy file. Default /data/bigbio_00/smcgreig/18S_nematode/18S_nematode_full_ncbi/18S_nhmmer_final_taxa_uniq.qza", default = "/data/bigbio_00/smcgreig/18S_nematode/18S_nematode_full_ncbi/18S_nhmmer_final_taxa_uniq.qza")


	return parser.parse_args()

################################################################################################

############################################################################################################################## Functions End ###############################################################################################################################

if __name__ == '__main__':
	main()
