import pytest
import json
import os
import yaml
from trimmomatic_test_randreads import get_adapters, get_trimmomatic_test_data

def uniform_quality(x, parameter = 0):
        return 0.1*x

def adapter_che_avanza(x):
        return 20

#esempio di reads file
#nome
#GTGTGGCGATGCTAGCATGCATCGATGACGACTGCTAGCATGCATCGATCGATCGATGACT
#+
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

def into_dict(stream):
	dictt = {}
	nome = stream.readline()
	while nome != "":
		seq = stream.readline()
		stream.readline()
		qual = stream.readline()
		dictt[nome[0:len(nome)-1]] = {"seq" : seq[0:len(seq)-1], "qual" : qual[0:len(qual)-1]}
		nome = stream.readline()
	return dictt

#esempio di contig file
#>NODE_2_length_213416_cov_13.619424
#GCAAATCTACAGTTCTGACAACCTGCTTCAATGGACAAAAGAGAGTTCCTTTGGAGCCGA
#GTATGGCTCACACGAAGGTGTTTGGGAATGTCCTGACCTTATAAAATTGCCTATTAGGGG

def contigs_into_dict(stream):
	r = {}
	name = stream.readline()
	sequence = ""
	linea = ""
	while not (name == ""):
		linea = stream.readline()
		while not (linea == "" or linea[0] == ">"):
			sequence += linea[0:len(linea)-2]
			linea = stream.readline()
		r[name] = sequence
		name = linea
	return r

def quality_back_casting(letter):
	letters = """!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"""
	return letters.find(letter)

class Trim_data:
	"""Quando viene costruito un oggetto con questa classe vengono creati i dati per i test di trimmomatic, nel caso non fossero già presenti, con le caratteristiche richieste
	I dati vengono sottoposti al trimming qualora non fosse già presente la loro versione trimmata
	Si noti che anche cambiando modificando i dati di partenza non verranno creati nuovi file trimmati, a meno che i vecchi non siano rimossi
	Questa classe espone 4 metodi pubblici che sono dictionary aventi struttura {nome_read : {"seq" : sequenza_delle_basi, "qual" : qualità_a_loro_associata}}:
		d_1 e d_2 sono i dati non trimmati
		tr_1 e tr_2 sono i dati risultanti dal trimming pair end di d_1 e d_2
	Il costruttore accetta come argomento il nome filename dei dati di partenza
	"""
	def __init__(self, data_name = "trim_test_data"):
		self.name = data_name
		self.d_1 = into_dict(open("test_dir/" + data_name + "_1.fastq"))
		self.d_2 = into_dict(open("test_dir/" + data_name + "_2.fastq"))
		os.system("gzip -k " + "test_dir/" + data_name + "_1.fastq")
		os.system("gzip -k " + "test_dir/" + data_name + "_2.fastq")
		os.system("snakemake " + "test_dir/" + data_name + "_1_trimP.fastq.gz --core 10 --use-conda")
		os.system("gzip -d " + "test_dir/" + data_name + "_1_trimP.fastq.gz")
		os.system("gzip -d " + "test_dir/" + data_name + "_2_trimP.fastq.gz")
		self.tr_1 = into_dict(open("test_dir/" + data_name + "_1_trimP.fastq"))
		self.tr_2 = into_dict(open("test_dir/" + data_name + "_2_trimP.fastq"))

	def __del__(self):
		os.system("rm " + "test_dir/" + self.name + "_1.fastq.gz " + "test_dir/" + self.name + "_2.fastq.gz")
		os.system("rm " + "test_dir/" + self.name + "_1_trimP.fastq " + "test_dir/" + self.name + "_2_trimP.fastq")

class Bwt_data:
	def __init__(self, data_name = "bwt_test_data"):
		self.name = data_name
		self.d_1 = into_dict(open("test_dir/" + data_name + "_1_trimP.fastq"))
		self.d_2 = into_dict(open("test_dir/" + data_name + "_2_trimP.fastq"))
		os.system("gzip -k " + "test_dir/" + data_name + "_1_trimP.fastq")
		os.system("gzip -k " + "test_dir/" + data_name + "_2_trimP.fastq")
		os.system("snakemake test_dir/" + data_name + "_bwt.allele_mapping_data.json --core 10 --use-conda")
		self.mapping = json.load(open("test_dir/" + data_name + "_bwt.allele_mapping_data.json", "r"))
		#mapping qui è una lista

	def __del__(self):
		os.system("rm test_dir/" + self.name + "_1_trimP.fastq.gz")
		os.system("rm test_dir/" + self.name + "_2_trimP.fastq.gz")
		os.system("rm test_dir/" + self.name + "_bwt.overall_mapping_stats.txt")
		os.system("rm test_dir/" + self.name + "_bwt.allele_mapping_data.txt")
		os.system("rm test_dir/" + self.name + "_bwt.allele_mapping_data.json")

class Spades_data:
	def __init__(self, data_name = "reads_for_test_spades"):
		self.name = data_name
		self.d_1 = into_dict(open("test_dir/" + data_name + "_1_trimP.fastq"))
		self.d_2 = into_dict(open("test_dir/" + data_name + "_2_trimP.fastq"))
		os.system("gzip -k " + "test_dir/" + data_name + "_1_trimP.fastq")
		os.system("gzip -k " + "test_dir/" + data_name + "_2_trimP.fastq")
		os.system("snakemake test_dir/" + data_name + "_assembly/contigs.fasta --core 10 --use-conda")
		self.contigs = contigs_into_dict(open("test_dir/" + data_name + "_assembly/contigs.fasta"))

	def __del__(self):
		os.system("rm test_dir/" + self.name + "_1_trimP.fastq.gz")
		os.system("rm test_dir/" + self.name + "_2_trimP.fastq.gz")
		os.system("rm -r test_dir/" + self.name + "_assembly")

class Main_data:
	def __init__(self, data_name = "main_test_data"):
		self.name = data_name
		os.system("snakemake test_dir/" + data_name + "_main.json --core 10 --use-conda")
		self.contigs = contigs_into_dict(open("test_dir/" + data_name + "_assembly/contigs.fasta"))
		self.mapping = json.load(open("test_dir/" + data_name + "_main.json", "r"))

	def __del__(self):
		os.system("rm " + self.name + "_main.*")


@pytest.fixture
def trim_dati():
	return Trim_data()

@pytest.fixture
def bwt_dati():
	return Bwt_data()

@pytest.fixture
def spades_dati():
	return Spades_data()

@pytest.fixture
def main_dati():
	return Main_data()

@pytest.fixture
def config():
	return yaml.safe_load(open("config.yaml", "r"))

#---------------------------------------------------------------------------------------------------------------------------------------------------------------
#	TEST FOR TRIMMOMATIC RULE
#---------------------------------------------------------------------------------------------------------------------------------------------------------------

@pytest.mark.trimmomatic
def test_len_trimmed_data(trim_dati):
	assert len(trim_dati.tr_1) == len(trim_dati.tr_2)

@pytest.mark.trimmomatic
def test_lower_len_for_trimmed_1(trim_dati):
	A = True
	for r in trim_dati.tr_1:
		A = A and (len(trim_dati.tr_1[r]["seq"]) <= len(trim_dati.d_1[r]["seq"]))
	assert A

@pytest.mark.trimmomatic
def test_lower_len_for_trimmed_2(trim_dati):
	A = True
	for r in trim_dati.tr_2:
		A = A and (len(trim_dati.tr_2[r]["seq"]) <= len(trim_dati.d_2[r]["seq"]))
	assert A

@pytest.mark.trimmomatic
def test_base_letters(trim_dati):
	A = True
	for r in trim_dati.tr_1:
		a = trim_dati.tr_1[r]["seq"].count('A')
		t = trim_dati.tr_1[r]["seq"].count("T")
		c = trim_dati.tr_1[r]["seq"].count("C")
		g = trim_dati.tr_1[r]["seq"].count("G")
		n = trim_dati.tr_1[r]["seq"].count("N")
		A = A and (len(trim_dati.tr_1[r]["seq"]) == a+t+c+g+n)
	for r in trim_dati.tr_2:
		a = trim_dati.tr_2[r]["seq"].count('A')
		t = trim_dati.tr_2[r]["seq"].count("T")
		c = trim_dati.tr_2[r]["seq"].count("C")
		g = trim_dati.tr_2[r]["seq"].count("G")
		n = trim_dati.tr_2[r]["seq"].count("N")
		A = A and (len(trim_dati.tr_2[r]["seq"]) == a+t+c+g+n)
	assert A

@pytest.mark.trimmomatic
def test_minlen(trim_dati, config):
	A = True
	for r in trim_dati.tr_1:
		A = A and len(trim_dati.tr_1[r]["seq"]) >= config["trim_minlen"]
	for r in trim_dati.tr_2:
		A = A and len(trim_dati.tr_2[r]["seq"]) >= config["trim_minlen"]
	assert A

@pytest.mark.trimmomatic
def test_leading(trim_dati, config):
	A = True
	for r in trim_dati.tr_1:
		A = A and quality_back_casting(trim_dati.tr_1[r]["qual"][0]) >= config["trim_leading_quality"]
	for r in trim_dati.tr_2:
		A = A and quality_back_casting(trim_dati.tr_2[r]["qual"][0]) >= config["trim_leading_quality"]
	assert A

@pytest.mark.trimmomatic
def test_trailing(trim_dati, config):
	A = True
	for r in trim_dati.tr_1:
		A = A and quality_back_casting(trim_dati.tr_1[r]["qual"][len(trim_dati.tr_1[r]["qual"])-1]) >= config["trim_trailing_quality"]
	for r in trim_dati.tr_2:
		A = A and quality_back_casting(trim_dati.tr_2[r]["qual"][len(trim_dati.tr_2[r]["qual"])-1]) >= config["trim_trailing_quality"]
	assert A

#--------------------------------------------------------------------------------------------------------------------------------------
#	TEST FOR BWT RULE
#--------------------------------------------------------------------------------------------------------------------------------------
@pytest.mark.bwt
def test_bwt_1(bwt_dati, config):
	assert len(bwt_dati.mapping) > 0

@pytest.mark.bwt
def test_mapped_reads_number(bwt_dati):
	A = True
	for r in bwt_dati.mapping:
		A = A and r["reads"]["all"] == r["reads"]["mapped"] + r["reads"]["unmapped"]
	assert A

@pytest.mark.bwt
def test_reference_length(bwt_dati):
	A = True
	for r in bwt_dati.mapping:
		A = A and r["reference"]["sequence_length"] == r["length_coverage"]["covered"] + r["length_coverage"]["uncovered"]
	assert A

@pytest.mark.bwt
def test_aro_accession_bwt(bwt_dati):
	A = True
	for r in bwt_dati.mapping:
		A = A and r["aro_accession"].isdigit()
	assert A

#-----------------------------------------------------------------------------------------------------------------------------------
#	TEST FOR SPADES RULE
#-------------------------------------------------------------------------------------------------------------------------------------
@pytest.mark.spades
def test_spades_non_vuoto(spades_dati, config):
	assert len(spades_dati.contigs) > 0

@pytest.mark.spades
def test_correct_data_in_contig_name(spades_dati):
	A = True
	for r in spades_dati.contigs:
		start = 8 + r.find("_length_")
		end = r.find("_cov_")
		try:
			int(r[start:end])
			float(r[end+5:])
		except ErrorValue:
			A = False
	assert A

@pytest.mark.spades
def test_base_letters_contigs(spades_dati):
	A = True
	for r in spades_dati.contigs:
		for char in spades_dati.contigs[r]:
			A = A and (char == 'A' or char == 'T' or char == 'C' or char == 'G' or char == 'N')
		if not A:
			break
	assert A

#-----------------------------------------------------------------------------------------------------------------------------------
#       TEST FOR MAIN RULE
#-------------------------------------------------------------------------------------------------------------------------------------
@pytest.mark.main
def test_main_1(main_dati, config):
	assert len(main_dati.mapping) > 0

@pytest.mark.main
def test_coverage(main_dati):
	A = True
	for r in main_dati.mapping:
		start = 5 + r.find("_cov_")
		end = r.find("_")
		try:
			float(r[start:end])
		except ValueError:
			A = False
	return A

@pytest.mark.main
def test_aro_accession(main_dati):
	A = True
	for r in main_dati.mapping:
		for f in main_dati.mapping[r]:
			A = A and main_dati.mapping[r][f]["ARO_accession"].isdigit()
	assert A
