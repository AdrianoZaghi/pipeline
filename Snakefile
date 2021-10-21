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
		rm {wildcards.sample}_bwt.sorted.length_100.*
		rm {wildcards.sample}_bwt.reference_mapping_stats.txt
		rm {wildcards.sample}_bwt.gene_mapping_data.txt
		rm {wildcards.sample}_bwt.artifacts_mapping_stats.txt
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
		"{sample}_2_trimP.fastq.gz",
	conda:
		"required_env.yaml"
	shell:
		"""
		trimmomatic PE -threads {config[threads]} -phred33 -trimlog trimlog.txt {input} {output[0]} {wildcards.sample}_1_trimS.fastq.gz {output[1]} {wildcards.sample}_2_trimS.fastq.gz ILLUMINACLIP:adapters/NEXTflex_96.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
		rm {wildcards.sample}_1_trimS.fastq.gz
		rm {wildcards.sample}_2_trimS.fastq.gz
		"""


dati = json.load(open("data_links.json"))
names = []
for s in dati.keys():
	names.append(s+"_1")
	names.append(s+"_2")
	names.append(s+"_1_trimP")
	names.append(s+"_2_trimP")
print(names)

rule reportrule:
#segui questo link per vedere come si usa la list comprehension https://www.w3schools.com/python/python_lists_comprehension.asp
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

#devo ordinare bene le voci di questo report in modo tale che si possa capire da che read_lib vengono, forse dovrò modificare il nome del file
#agiungere sequence length e overrepresented sequences?

rule get_data:
	"""
	scarica i dati che sono presenti nel config file
	tecnicamente richiede come input il file data_links.json, ma è meglio limitarsi ad usarlo solo nella linea di comando
	altrimenti ogni volta che modifico il file, visto che l'ultima modifica all'input è più recente dell'output, anche tutti i campioni già scaricati vengono riscaricati
	"""
#	input:
#		"data_links.json"
	output:
		"{sample}_1.fastq.gz",
		"{sample}_2.fastq.gz"
	run:
		f = open("data_links.json")
		dati = json.load(f)
		shell("wget --tries=10 " + dati[wildcards.sample][0])
		shell("wget --tries=10 " + dati[wildcards.sample][1])
#-------------------------------------------------------------------------------

#capire come ottenere il report
#implementa un modo di scartare le contigs troppo corte dopo l'assembly
