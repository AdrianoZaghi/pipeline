import pytest
import json
import os
import yaml
from trimmomatic_test_randreads import get_adapters, get_trimmomatic_test_data

#-----------------------------------------------------
#	ROUTINE FUNCTIONS
#-----------------------------------------------------

#example of entry of .fastq file
#name
#GTGTGGCGATGCTAGCATGCATCGATGACGACTGCTAGCATGCATCGATCGATCGATGACT
#+
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

def into_dict(stream):
	"""This function is used to organize in a dictionariey the informations in a .fastq file
	It's argument is an Istream object, opened on a .fastq file
	It returns a dictionary with the following structure:
		{<read name> : {"seq" : <string with the DNA sequence>, "qual" : <string with quality base score>}}
	"""
	dictt = {}
	name = stream.readline()
	while name != "":
		seq = stream.readline()
		stream.readline()
		qual = stream.readline()
		dictt[name[0:len(name)-1]] = {"seq" : seq[0:len(seq)-1], "qual" : qual[0:len(qual)-1]}
		name = stream.readline()
	return dictt

#example of entry of contig.fasta file
#>NODE_2_length_213416_cov_13.619424
#GCAAATCTACAGTTCTGACAACCTGCTTCAATGGACAAAAGAGAGTTCCTTTGGAGCCGA
#GTATGGCTCACACGAAGGTGTTTGGGAATGTCCTGACCTTATAAAATTGCCTATTAGGGG

def contigs_into_dict(stream):
	"""This function is used to organize in a dictionariey the informations in a contigs.fasta file
        It's argument is an Istream object, opened on a contig.fastq file
        It returns a dictionary with the following structure:
                {<contig name> : <contig sequence>}
	"""
	r = {}
	name = stream.readline()
	sequence = ""
	line = ""
	while not (name == ""):
		line = stream.readline()
		while not (line == "" or line[0] == ">"):
			sequence += line[0:len(line)-2]
			line = stream.readline()
		r[name] = sequence
		name = line
	return r

def quality_back_casting(letter):
	"""This function is used to convert the quality characters used in .fasta files in numeric score
	"""
	letters = """!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"""
	return letters.find(letter)

#-------------------------------------------------------------------------------
#	CLASSES TO ORGAINZE RULES INPUT AND OUTPUT
#-------------------------------------------------------------------------------

class Trim_data:
	"""When an object of this class is constructed, is executed the trimming of the data characterized my "data_name" run accession, if trimmed versions are not already present in the working directory
	The data of trimmed and untrimmed read libraries are collected in 4 public metods. Those are dictionaryes in which each element have the following structure {<read_name>: {"seq" : >base sequence>, "qual" : <quality score> }}
		.d_1 and .d_2 contains reads from untrimmed data
		.tr_1 and .tr_2 contains reads from trimmed data
	"data_name" standard value is the run accession of the data created for the thest, present in the test_dir/
	"""
	def __init__(self, data_name = "trim_test_data"):
		self.name = data_name
		os.system("gzip -k " + "test_dir/" + data_name + "_1.fastq")
		os.system("gzip -k " + "test_dir/" + data_name + "_2.fastq")
		os.system("snakemake " + "test_dir/" + data_name + "_1_trimP.fastq.gz --core 10 --use-conda")
		os.system("gzip -d " + "test_dir/" + data_name + "_1_trimP.fastq.gz")
		os.system("gzip -d " + "test_dir/" + data_name + "_2_trimP.fastq.gz")
		with open("test_dir/" + data_name + "_1.fastq", mode = 'r', encoding = "utf-8") as dati_1, open("test_dir/" + data_name + "_2.fastq", mode = 'r', encoding = "utf-8") as dati_2, open("test_dir/" + data_name + "_1_trimP.fastq", mode = 'r', encoding = "utf-8") as trimmed_1, open("test_dir/" + data_name + "_2_trimP.fastq", mode = 'r', encoding = "utf-8") as trimmed_2:
			self.d_1 = into_dict(dati_1)
			self.d_2 = into_dict(dati_2)
			self.tr_1 = into_dict(trimmed_1)
			self.tr_2 = into_dict(trimmed_2)

#	def __del__(self):
#		os.system("rm " + "test_dir/" + self.name + "_1.fastq.gz " + "test_dir/" + self.name + "_2.fastq.gz")
#		os.system("rm " + "test_dir/" + self.name + "_1_trimP.fastq " + "test_dir/" + self.name + "_2_trimP.fastq")

class Bwt_data:
	"""When an object of this class is constructed, the data characterized by "data_name" run accession are aligned with bowtie2. This won't happen if there is already, in the working directory, a file that contains the alignement result
	The aligment have to be performed on trimmed data, and if there arent _trimP versions of the given run accession, those will be created performing the trimming also
	The standard value of "data_name" corresponds to the run accession of the data created for the testing in test_dir/
	The created object have 3 public methods, in which are organized input and output data of bwt alignment:
		.d_1 and .d_2 contains the information about trimmed reads used in the alignemnt. Those have structure: {<read name> : {"seq" : <read sequence>, "qual" : <read quality sequence>}}
		.mapping is a json list that contains mapping results
	"""
	def __init__(self, data_name = "bwt_test_data"):
		self.name = data_name
		os.system("gzip -k " + "test_dir/" + data_name + "_1_trimP.fastq")
		os.system("gzip -k " + "test_dir/" + data_name + "_2_trimP.fastq")
		os.system("snakemake test_dir/" + data_name + "_bwt.allele_mapping_data.json --core 10 --use-conda")
		with open("test_dir/" + data_name + "_1_trimP.fastq", mode = 'r', encoding = "utf-8") as dati_1, open("test_dir/" + data_name + "_2_trimP.fastq", mode = 'r', encoding = "utf-8") as dati_2, open("test_dir/" + data_name + "_bwt.allele_mapping_data.json", mode = 'r', encoding = "utf-8") as bwt:
			self.d_1 = into_dict(dati_1)
			self.d_2 = into_dict(dati_2)
			self.mapping = json.load(bwt)
		#mapping qui è una lista

#	def __del__(self):
#		os.system("rm test_dir/" + self.name + "_bwt.overall_mapping_stats.txt")
#		os.system("rm test_dir/" + self.name + "_bwt.allele_mapping_data.txt")
#		os.system("rm test_dir/" + self.name + "_bwt.allele_mapping_data.json")

class Spades_data:
	"""When an object of this class is constructed, the data characterized by "data_name" run accession are assembled with Spades_rule, if is not already present an aligment for that run accession in the working directory
	The data used for the assembly have to be alredy trimmed, or with a _trimP filename, otherwise also the trimming is executed.
	The standard value of "data_name" corresponds to the run accession of the data created for the testing in test_dir/
	This class contains 3 public methods in which are orgainzed input and output data about assembly:
		.d_1 and .d_2 contains the information about trimmed reads used in the assembly. Those have structure: {<read name> : {"seq" : <read sequence>, "qual" : <read quality sequence>}}
		.contigs contains informations about the contigs obtained in the assembly. Is a dict with structure {<contig name> : <contig sequence>}
	"""
	def __init__(self, data_name = "reads_for_test_spades"):
		self.name = data_name
		os.system("gzip -k " + "test_dir/" + data_name + "_1_trimP.fastq")
		os.system("gzip -k " + "test_dir/" + data_name + "_2_trimP.fastq")
		os.system("snakemake test_dir/" + data_name + "_assembly/contigs.fasta --core 10 --use-conda")
		with open("test_dir/" + data_name + "_1_trimP.fastq", mode = 'r', encoding = "utf-8") as dati_1, open("test_dir/" + data_name + "_2_trimP.fastq", mode = 'r', encoding = "utf-8") as dati_2, open("test_dir/" + data_name + "_assembly/contigs.fasta", mode = 'r', encoding = "utf-8") as contigs:
			self.contigs = contigs_into_dict(contigs)
			self.d_1 = into_dict(dati_1)
			self.d_2 = into_dict(dati_2)
#	def __del__(self):
#		os.system("rm -r test_dir/" + self.name + "_assembly")

class Main_data:
	"""When an object of this class is constructed, the contigs characterized by "data_name" run accession are aligned, if there is not already in the working directory a file with alignemnt output
	If there isn't any <run accession>_assembly/contigs.fasta that matchs "data_name", the pipeline is executed in order to obtain them, if possible
	The standard value of "data_name" refers to the data created for the testing and those are located in the test_dir/ folder
	Those objects have 2 public methods:
		.contigs ontains informations about the contigs used for the alignemnt. Is a dict with structure {<contig name> : <contig sequence>}
		.mapping is a json list that contains mapping results
	"""
	def __init__(self, data_name = "main_test_data"):
		self.name = data_name
		os.system("snakemake test_dir/" + data_name + "_main.json --core 10 --use-conda")
		with open("test_dir/" + data_name + "_assembly/contigs.fasta", mode = "r", encoding = "utf-8") as contigs, open("test_dir/" + data_name + "_main.json", mode = "r", encoding = "utf-8") as main:
			self.contigs = contigs_into_dict(contigs)
			self.mapping = json.load(main)

#	def __del__(self):
#		os.system("rm " + self.name + "_main.*")

#----------------------------------------
#	FIXTURES
#-------------------------------------------------------------

@pytest.fixture
def trim_dati(name = "trim_test_data"):
	yield Trim_data(name)
	os.system("rm " + "test_dir/" + name + "_1.fastq.gz " + "test_dir/" + name + "_2.fastq.gz")
	os.system("rm " + "test_dir/" + name + "_1_trimP.fastq " + "test_dir/" + name + "_2_trimP.fastq")
	os.system("rm " + "test_dir/" + name + "_1_trimP.fastq.gz " + "test_dir/" + name + "_2_trimP.fastq.gz")

@pytest.fixture
def bwt_dati(name = "bwt_test_data"):
	yield Bwt_data(name)
	os.system("rm test_dir/" + name + "_bwt.overall_mapping_stats.txt")
	os.system("rm test_dir/" + name + "_bwt.allele_mapping_data.txt")
	os.system("rm test_dir/" + name + "_bwt.allele_mapping_data.json")

@pytest.fixture
def spades_dati(name = "reads_for_test_spades"):
	yield Spades_data(name)
	os.system("rm -r test_dir/" + name + "_assembly")

@pytest.fixture
def main_dati(name = "main_test_data"):
	yield Main_data(name)
	os.system("rm " + name + "_main.*")

@pytest.fixture
def config():
	with open("config.yaml", mode = 'r', encoding = "utf-8") as stream:
		yield yaml.safe_load(stream)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------
#	TESTS FOR TRIMMOMATIC RULE
#---------------------------------------------------------------------------------------------------------------------------------------------------------------

@pytest.mark.trimmomatic
def test_len_trimmed_data(trim_dati):
	"""Ensure the pair end read libraries to have the same read number after the trimming"""
	assert len(trim_dati.tr_1) == len(trim_dati.tr_2)

@pytest.mark.trimmomatic
def test_lower_len_for_trimmed(trim_dati):
	"""Ensure trimmed sequences to have length <= respect original sequences lengt"""
	assert all(len(trim_dati.tr_1[r]["seq"]) <= len(trim_dati.d_1[r]["seq"]) for r in trim_dati.tr_1) and all(len(trim_dati.tr_2[r]["seq"]) <= len(trim_dati.d_2[r]["seq"]) for r in trim_dati.tr_2)

@pytest.mark.trimmomatic
def test_base_letters(trim_dati):
	"""Ensure only DNA base letters are contained in read sequences"""
	assert all(set(trim_dati.tr_1[r]["seq"]) < {'A','a','T','t','C','c','G','g','T','t','N','n'} for r in trim_dati.tr_1) and all(set(trim_dati.tr_2[r]["seq"]) < {'A','a','T','t','C','c','G','g','T','t','N','n'} for r in trim_dati.tr_2)

@pytest.mark.trimmomatic
def test_minlen(trim_dati, config):
	"""Ensure all trimmed sequences to have length >= to the minlen value specified in config.jaml file"""
	assert all(len(trim_dati.tr_2[r]["seq"]) >= config["trim_minlen"] for r in trim_dati.tr_2) and all(len(trim_dati.tr_1[r]["seq"]) >= config["trim_minlen"] for r in trim_dati.tr_1)

@pytest.mark.trimmomatic
def test_leading(trim_dati, config):
	"""Ensure first base of trimmed sequences to have quality >= then leading value specified in config.jaml file"""
	assert all(quality_back_casting(trim_dati.tr_2[r]["qual"][0]) >= config["trim_leading_quality"] for r in trim_dati.tr_2) and all(quality_back_casting(trim_dati.tr_1[r]["qual"][0]) >= config["trim_leading_quality"] for r in trim_dati.tr_1)

@pytest.mark.trimmomatic
def test_trailing(trim_dati, config):
	"""Ensure last base of trimmed sequences to have quality >= then trailing value specified in config.jaml file"""
	assert all(quality_back_casting(trim_dati.tr_2[r]["qual"][len(trim_dati.tr_2[r]["qual"])-1]) >= config["trim_trailing_quality"] for r in trim_dati.tr_2) and all(quality_back_casting(trim_dati.tr_1[r]["qual"][len(trim_dati.tr_1[r]["qual"])-1]) >= config["trim_trailing_quality"] for r in trim_dati.tr_1)

#--------------------------------------------------------------------------------------------------------------------------------------
#	TESTS FOR BWT RULE
#--------------------------------------------------------------------------------------------------------------------------------------
@pytest.mark.bwt
def test_bwt_1(bwt_dati, config):
	"""Ensure rgi bwt output not to be empthy"""
	assert len(bwt_dati.mapping) > 0

@pytest.mark.bwt
def test_mapped_reads_number(bwt_dati):
	"""Ensure coherence between total read number and the the sum of perfectly mapped and flanking mapped reads"""
	assert all(int(r["reads"]["all"]) == int(r["reads"]["mapped"]) + int(r["reads"]["unmapped"]) for r in bwt_dati.mapping)

@pytest.mark.bwt
def test_reference_length(bwt_dati):
	"""Ensure coherence between total length of reference and the sum fo covered and uncovered refeence length"""
	assert all(int(r["reference"]["sequence_length"]) == int(r["length_coverage"]["covered"]) + int(r["length_coverage"]["uncovered"]) for r in bwt_dati.mapping)

@pytest.mark.bwt
def test_aro_accession_bwt(bwt_dati):
	"""Ensure the aro accession to be noly digit value"""
	assert all(r["aro_accession"].isdigit() for r in bwt_dati.mapping)

#-----------------------------------------------------------------------------------------------------------------------------------
#	TESTS FOR SPADES RULE
#-------------------------------------------------------------------------------------------------------------------------------------
@pytest.mark.spades
def test_spades_non_vuoto(spades_dati, config):
	"""Ensure spades output not to be empthy"""
	assert len(spades_dati.contigs) > 0

@pytest.mark.spades
def test_correct_data_in_contig_name(spades_dati):
	"""Ensure the coverage value to be contained in contig names"""
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
	"""Ensure all contigs letters to be DNA letters"""
	assert all(set(spades_dati.contigs[r]) <= {'A','a','T','t','C','c','G','g','T','t','N','n'} for r in spades_dati.contigs)

#-----------------------------------------------------------------------------------------------------------------------------------
#       TESTS FOR MAIN RULE
#-------------------------------------------------------------------------------------------------------------------------------------
@pytest.mark.main
def test_main_1(main_dati, config):
	"""Ensure rgi main output not to be empthy"""
	assert len(main_dati.mapping) > 0

@pytest.mark.main
def test_coverage(main_dati):
	"""Ensure coverage value to be present in mapped contigs name"""
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
def test_aro_accession_main(main_dati):
	"""Ensure the aro accession to be noly digit value"""
	assert all(((main_dati.mapping[r][f]["ARO_accession"].isdigit() for f in main_dati.mapping[r]) for r in main_dati.mapping))
