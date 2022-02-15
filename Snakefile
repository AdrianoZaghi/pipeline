
import json

configfile: "config.yaml"

dati = json.load(open("data_links.json"))

#The names that lables data_links.json entries are used to compose the report and should match the run_accession
names = []
for s in dati.keys():
	names.append(s+"_1")
	names.append(s+"_2")
	names.append(s+"_1_trimP")
	names.append(s+"_2_trimP")
		

rule main:
	"""This rule perform the alignment of the contigs to Card and WildCard using rgi main. In the config.yaml file can be set the alignemnt tool
	Input file	are the databases references, which are set up in the /rgi repo during installation, and the contig file resulting from the assembly
	Output file	is a json file whic contains all the information resulting from the assembly. Also a txt version is produced. Those files are found in the working directory at the end of execution
	"""
	input:
		"{sample}_assembly/contigs.fasta",
		"rgi/card_annotation.log",
		"rgi/wildcard_annotation.log"
	output:
		"{sample}_main.json"
	conda:
		"required_env.yaml"
	shell:
		"""
		cd rgi
		rgi main -i ../{input[0]} -o {wildcards.sample}_main -t contig -a {config[rgi_main_alignment_tool]} --clean -n {config[threads]} -d wgs --split_prodigal_jobs
		mv {wildcards.sample}_main.txt ../{wildcards.sample}_main.txt
		mv {wildcards.sample}_main.json ../{wildcards.sample}_main.json
		cd ..
		"""

rule assembly:
	"""This rule compose the assembly of the trimmed read libraries using spades. In the config.yaml it can be set a limit for RAM usage and the K paramether
	Input file	are the 2 pair end read libraries, both trimmed
	Output file	are contained in the directory {run_accession}_assembly. It is created ad hoc in the working direcotry. Beside some information about assembly process, the most important file is contigs.fasta
	"""
	input:
		"{sample}_1_trimP.fastq.gz",
		"{sample}_2_trimP.fastq.gz"
	output:
		"{sample}_assembly/contigs.fasta"
	conda:
		"required_env.yaml"
	shell:
		"spades.py -1 {input[0]} -2 {input[1]} -t {config[threads]} -m {config[ram]} -k {config[assembly_k_parameter]} -o {wildcards.sample}_assembly"

rule rgi_bwt:
	"""This rule perform the alignement of the trimmed reads to Card and WildCard databases using rgi bwt
	Input file	are the databases references, which are set up during installation and the 2 pair end read libraries, both trimmed
	Output file	is a json file which contains all the information resulting from the alignment. Also a txt version is produced.
	"""
	input:
		"{sample}_1_trimP.fastq.gz",
		"{sample}_2_trimP.fastq.gz",
		"rgi/card_annotation.log",
		"rgi/wildcard_annotation.log"
	output:
		"{sample}_bwt.allele_mapping_data.txt",
		"{sample}_bwt.allele_mapping_data.json"
	conda:
		"required_env.yaml"
	shell:
		"""
		cd rgi
		rgi bwt -1 ../{input[0]} -2 ../{input[1]} -a bowtie2 -n {config[threads]} --clean --include_wildcard -o ../{wildcards.sample}_bwt --local
		cd ..
		rm {wildcards.sample}_bwt.sorted.length_100.*
		rm {wildcards.sample}_bwt.reference_mapping_stats.txt
		rm {wildcards.sample}_bwt.gene_mapping_data.txt
		rm {wildcards.sample}_bwt.artifacts_mapping_stats.txt
		"""

rule rgi_env_setup:
	"""This rule download the references to the databases Card and WildCard and load both of them in rgi, so that those can be used for the alignment
	It don't require any input file since this operation have to be done in order to set up the envoirement for the pipeline to work
	Output file	are contained in the new created folder rgi/ , those are the references to the two databases
	"""
	output:
		"rgi/card_annotation.log",
		"rgi/wildcard_annotation.log"
	conda:
		"required_env.yaml"
	shell:
		"""
		mkdir -p rgi
		cd rgi
		#CARD DATABASE
		wget -O data https://card.mcmaster.ca/latest/data
		tar -xvf data ./card.json
		rgi load --card_json card.json --local
		rgi card_annotation -i card.json > card_annotation.log 2>&1
		rgi load -i card.json --card_annotation card_database* --local
		#WILDCARD DATABASE
		wget -O wildcard_data.tar.bz2 https://card.mcmaster.ca/latest/variants
		mkdir -p wildcard
		tar -xjf wildcard_data.tar.bz2 -C wildcard
		gunzip wildcard/*.gz
		rgi wildcard_annotation -i wildcard --card_json card.json -v_non_specificata > wildcard_annotation.log 2>&1
		rgi load --wildcard_annotation wildcard_database* --wildcard_index wildcard/index-for-model-sequences.txt --card_annotation card_database* --local
		cd ..
		"""

rule trimmomatic:
	"""This rule perform the pair end trimming of the read libraries eliminating the NextFles-96 adapters. Is used trimmomatic and all of it's parameters can be set in the config.yaml file.
	Input files	are the non trimmed pair ended read libraries
	Output files	are the trimmed pair ended read libraries
	"""
	input:
		"{sample}_1.fastq.gz",
		"{sample}_2.fastq.gz"

	output:
		"{sample}_1_trimP.fastq.gz",
		"{sample}_2_trimP.fastq.gz",
	conda:
		"required_env.yaml"
	shell:
		"""
		trimmomatic PE -threads {config[threads]} -phred33 -trimlog trimlog.txt {input} {output[0]} {wildcards.sample}_1_trimS.fastq.gz {output[1]} {wildcards.sample}_2_trimS.fastq.gz ILLUMINACLIP:adapters/NEXTflex_96.fa:{config[trim_seed_mismatch]}:{config[trim_palindrome_clip_threshold]}:{config[trim_simple_clip_threshold]} LEADING:{config[trim_leading_quality]} TRAILING:{config[trim_trailing_quality]} SLIDINGWINDOW:{config[trim_window_size]}:{config[trim_quality_required]} MINLEN:{config[trim_minlen]}
		rm {wildcards.sample}_1_trimS.fastq.gz
		rm {wildcards.sample}_2_trimS.fastq.gz
		"""

rule reportrule:
	"""This rule create a report about the read libraries mentioned in the data_links.json file. It regards both the trimmed and untrimmed version
	Input files	are the read libraries
	Output files	are an empthy file report.out and the file report.html, which contain important stathistics about the read libraries, provided by FastQC
	"""
	input:
		[name + ".fastq.gz" for name in names]
	output:
		"report.out",
		[report("FastQC_report/" + name + "_fastqc/Images/per_base_quality.png", category = name) for name in names],
		[report("FastQC_report/" + name + "_fastqc/Images/adapter_content.png", category = name) for name in names],
		[report("FastQC_report/" + name + "_fastqc/Images/sequence_length_distribution.png", category = name) for name in names]
	conda:
		"required_env.yaml"
	shell:
		"""
		touch report.out
		mkdir -p FastQC_report
		fastqc {input} --outdir=FastQC_report 
		cd FastQC_report
		unzip -o \*.zip
		cd ..
		"""

rule get_data:
	"""This rule download all the reads libraries mentioned in the data_links.json file
	Even if the data_links.json file is required for this rule execution, it is not an input file. In this the rule is not executed again each time the links are modified, but just when a non existing sample is required
	Output data are the two pair end read libraries
	"""
	output:
		"{sample}_1.fastq.gz",
		"{sample}_2.fastq.gz"
	run:
		f = open("data_links.json")
		dati = json.load(f)
		shell("wget --tries=10 " + dati[wildcards.sample][0])
		shell("wget --tries=10 " + dati[wildcards.sample][1])
