import json
from urllib import request

configfile: "config.yaml"

rule rgi_bwt:
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
		"""

rule rgi_env_setup:
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
	"""
	trimming rule
	in questa parte del workflow vengono rimossi gli adapter dalle reads 
	"""
	input:
		"{sample}_1.fastq.gz",
		"{sample}_2.fastq.gz"

	output:
		"{sample}_1_trimP.fastq.gz",
		"{sample}_2_trimP.fastq.gz"
	conda:
		"required_env.yaml"
	shell:
		"""
		trimmomatic PE -threads {config[threads]} -phred33 -trimlog trimlog.txt {input} {output[0]} {wildcards.sample}_1_trimS.fastq.gz {output[1]} {wildcards.sample}_2_trimS.fastq.gz ILLUMINACLIP:adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
		rm {wildcards.sample}_1_trimS.fastq.gz
		rm {wildcards.sample}_2_trimS.fastq.gz
		"""
rule FastQC:
	"""
	rule per il controllo di qualità delle reads non trimmate
	i dati che produce vanno nel report
	"""

	input:
		"{sample}_1.fastq.gz",
		"{sample}_2.fastq.gz"
	output:
		report("../FastQC_report/{sample}_1_fastqc/Images/per_base_quality.png", category = "original reads"),
		report("../FastQC_report/{sample}_2_fastqc/Images/per_base_quality.png", category = "original reads")
	conda:
		"required_env.yaml"
	shell:
		"""
		mkdir -p FastQC_report
		fastqc {input[1]} --outdir=FastQC_report
		fastqc {input[2]} --outdir=FastQC_report
		unzip FastQC_report/{wildcards.sample}_2_fastqc.zip
		unzip FastQC_report/{wildcards.sample}_1_fastqc.zip
		"""

rule get_data:
	"""
	scarica i dati che sono presenti nel config file
	"""
	input:
		"data_links.json"
	output:
		"{sample}_1.fastq.gz",
		"{sample}_2.fastq.gz"
	run:
		f = open("data_links.json")
		dati = json.load(f)
		shell("wget " + dati[wildcards.sample][0])
		shell("wget " + dati[wildcards.sample][1])
#-------------------------------------------------------------------------------

#capire come ottenere il report
#implementa un modo di scartare le contigs troppo corte dopo l'assembly
