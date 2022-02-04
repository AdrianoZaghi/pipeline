import pytest
import json
import os
import yaml
from trimmomatic_test_randreads import get_adapters, get_trimmomatic_test_data

#-----------------------------------------------------
#	ROUTINE FUNCTIONS
#-----------------------------------------------------

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

#-------------------------------------------------------------------------------
#	CLASSI PER TESTATE UNA RULE DELLA PIPELINE
#-------------------------------------------------------------------------------

class Trim_data:
	"""Quando viene costruito un oggetto con questa classe i dati predisposti per il test vengono sottoposti al trimming, qualora non fosse già presente la loro versione trimmata
	Questa classe espone 4 metodi pubblici che sono dictionary di dictionary, ogni elemento ha la struttura {"nome_read" : {"seq" : sequenza_delle_basi, "qual" : qualità_a_loro_associata}}
		d_1 e d_2 contengono i dati non trimmati
		tr_1 e tr_2 contengono i dati risultanti dal trimming pair end di d_1 e d_2
	Il costruttore accetta come argomento il nome filename dei dati di partenza, ma questo ha un valore standard
	Il distruttore è programmato per distruggere file prodotti nella lavorazione e anche gli output del trimming alla fine del test NON SO SE VADA BENE DISTRUGGERLI
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
	"""Quando vienne costruito un oggetto con questa classe i dati predisposti per il test vengono allineati con la bwt_rule, qualora non fosse già presente un allienamento riferito a quei dati
	Devono essere presenti delle libreie di reads già timmate (o contrassegnate da trimP), altrimenti verranno prodotte utilizzando la trimmomatic rule da reads non trimmate coerenti col filename
	Questa classe espone 3 metodi pubblici.
		Due dictionary di dictionary dalla struttura {"nome_read" : {"seq" : sequenza_delle_basi, "qual" : qualità_a_loro_associata}}, contenenti di dati di input (d_1, d_2)
		Una json list che contiene i risultati del mapping (mapping)
	Il costruttore accetta come argomento il nome filename dei dati di partenza, ma questo ha un valore standard
	Il distruttore è programmato per eliminare i dati dell'allineamento alla fine del test NON SO SE VADA BENE DISTRUGGERLI
	"""
	def __init__(self, data_name = "bwt_test_data"):
		self.name = data_name
		self.d_1 = into_dict(open("test_dir/" + data_name + "_1_trimP.fastq"))
		self.d_2 = into_dict(open("test_dir/" + data_name + "_2_trimP.fastq"))
		os.system("gzip -k " + "test_dir/" + data_name + "_1_trimP.fastq")
		os.system("gzip -k " + "test_dir/" + data_name + "_2_trimP.fastq")
		os.system("snakemake test_dir/" + data_name + "_bwt.allele_mapping_data.json --core 10 --use-conda")
		self.mapping = json.load(open("test_dir/" + data_name + "_bwt.allele_mapping_data.json", "r"))
		#mapping qui è una lista

#	def __del__(self):
#		os.system("rm test_dir/" + self.name + "_bwt.overall_mapping_stats.txt")
#		os.system("rm test_dir/" + self.name + "_bwt.allele_mapping_data.txt")
#		os.system("rm test_dir/" + self.name + "_bwt.allele_mapping_data.json")

class Spades_data:
	"""Quando vienne costruito un oggetto con questa classe i dati predisposti per il test vengono assemblati con Spades_rule, qualora non fosse già presente un allienamento riferito a quei dati
	Devono essere presenti delle libreie di reads già timmate (o contrassegnate da trimP), altrimenti verranno prodotte utilizzando la trimmomatic rule da reads non trimmate coerenti col filename
	Questa classe contiene 3 metodi pubblici:
		Due dictionary di dictionary dalla struttura {{"nome read" : {"seq" : sequenza_delle_basi, "qual" : qualità_a_loro_associata}}}, contenenti di dati di input (d_1, d_2)
		Una dictionary dalla struttura {"nome contig": sequenza}, contenente i dati sulle contigs prodotte dall'assembly
	Il costruttore accetta come argomento il nome filename dei dati di partenza, ma questo ha un valore standard
	Il distruttore è programmato per eliminare i dati dell'assembly alla fine del test NON SO SE VADA BENE DISTRUGGERLI
	"""
	def __init__(self, data_name = "reads_for_test_spades"):
		self.name = data_name
		self.d_1 = into_dict(open("test_dir/" + data_name + "_1_trimP.fastq"))
		self.d_2 = into_dict(open("test_dir/" + data_name + "_2_trimP.fastq"))
		os.system("gzip -k " + "test_dir/" + data_name + "_1_trimP.fastq")
		os.system("gzip -k " + "test_dir/" + data_name + "_2_trimP.fastq")
		os.system("snakemake test_dir/" + data_name + "_assembly/contigs.fasta --core 10 --use-conda")
		self.contigs = contigs_into_dict(open("test_dir/" + data_name + "_assembly/contigs.fasta"))

#	def __del__(self):
#		os.system("rm -r test_dir/" + self.name + "_assembly")

class Main_data:
	"""Quando vienne costruito un oggetto con questa classe i dati predisposti per il test vengono allineati con rule_main, qualora non fosse già presente un allienamento riferito a quei dati
	Devono essere già presenti delle contigs nella cartella di assembly, altrimenti verranno prodotte utilizzando la spades_rule da reads trimmate coerenti con il filename
	Questa classe contiene 2 metodi pubblici:
		Una dictionary dalla struttura {"nome contig": sequenza}, contenente i dati sulle contigs usate per l'allineamento
		Un json dataset contentente i risultati del mapping
	Il costruttore accetta come argomento il nome filename dei dati di partenza, ma questo ha un valore standard
	Il distruttore è programmato per eliminare i dati dell' allineamento alla fine del test NON SO SE VADA BENE DISTRUGGERLI
	"""
	def __init__(self, data_name = "main_test_data"):
		self.name = data_name
		os.system("snakemake test_dir/" + data_name + "_main.json --core 10 --use-conda")
		self.contigs = contigs_into_dict(open("test_dir/" + data_name + "_assembly/contigs.fasta"))
		self.mapping = json.load(open("test_dir/" + data_name + "_main.json", "r"))

#	def __del__(self):
#		os.system("rm " + self.name + "_main.*")

#----------------------------------------
#	FIXTURES
#-------------------------------------------------------------

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
	"""Assicura che le read libraries in pair and abbiano lo stesso numero di reads alla fine del trimming"""
	assert len(trim_dati.tr_1) == len(trim_dati.tr_2)

@pytest.mark.trimmomatic
def test_lower_len_for_trimmed_1(trim_dati):
	"""Controlla che tutte le sequenze trimmate abbiano una lunghezza <= a quella delle reads originali (forward sequence)"""
	A = True
	for r in trim_dati.tr_1:
		A = A and (len(trim_dati.tr_1[r]["seq"]) <= len(trim_dati.d_1[r]["seq"]))
	assert A

@pytest.mark.trimmomatic
def test_lower_len_for_trimmed_2(trim_dati):
	"""Come sopra, ma reverse sequence"""
	A = True
	for r in trim_dati.tr_2:
		A = A and (len(trim_dati.tr_2[r]["seq"]) <= len(trim_dati.d_2[r]["seq"]))
	assert A

@pytest.mark.trimmomatic
def test_base_letters(trim_dati):
	"""Controlla che nelle sequenze siano contenute solo le lettere relative alle basi del DNA"""
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
	"""Controlla che, dopo il trimming, siano state escluse le sequenze aventi una lunghezza inferiore a quella specificata nelle configurazioni di trimmomatic"""
	A = True
	for r in trim_dati.tr_1:
		A = A and len(trim_dati.tr_1[r]["seq"]) >= config["trim_minlen"]
	for r in trim_dati.tr_2:
		A = A and len(trim_dati.tr_2[r]["seq"]) >= config["trim_minlen"]
	assert A

@pytest.mark.trimmomatic
def test_leading(trim_dati, config):
	"""Controlla che la prima base delle sequenze trimmate abbia una qulità >= a quella specificata nelle configurazioni di trimmomatic"""
	A = True
	for r in trim_dati.tr_1:
		A = A and quality_back_casting(trim_dati.tr_1[r]["qual"][0]) >= config["trim_leading_quality"]
	for r in trim_dati.tr_2:
		A = A and quality_back_casting(trim_dati.tr_2[r]["qual"][0]) >= config["trim_leading_quality"]
	assert A

@pytest.mark.trimmomatic
def test_trailing(trim_dati, config):
	"""Controlla che l'ultima base delle sequenze trimmate abbia una qualità >= a quella specificata nelle configurazioni di trimmomatic"""
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
	"""Controlla che l'output di rgi_bwt non sia vuoto"""
	assert len(bwt_dati.mapping) > 0

@pytest.mark.bwt
def test_mapped_reads_number(bwt_dati):
	"""Controlla la coerenza della quantità di reads mappate perfettamente e mappate con flanking, con il numero di reads totali mappate"""
	A = True
	for r in bwt_dati.mapping:
		A = A and int(r["reads"]["all"]) == int(r["reads"]["mapped"]) + int(r["reads"]["unmapped"])
	assert A

@pytest.mark.bwt
def test_reference_length(bwt_dati):
	"""Controlla che la lunghezza del riferimento coperta dalla mappatura, più quella non coperta, sia uguale alla lunghezza totale del riferimento"""
	A = True
	for r in bwt_dati.mapping:
		A = A and int(r["reference"]["sequence_length"]) == int(r["length_coverage"]["covered"]) + int(r["length_coverage"]["uncovered"])
	assert A

@pytest.mark.bwt
def test_aro_accession_bwt(bwt_dati):
	"""Controlla che il campo contrassegnato da aro accession sia numerico"""
	A = True
	for r in bwt_dati.mapping:
		A = A and r["aro_accession"].isdigit()
	assert A

#-----------------------------------------------------------------------------------------------------------------------------------
#	TEST FOR SPADES RULE
#-------------------------------------------------------------------------------------------------------------------------------------
@pytest.mark.spades
def test_spades_non_vuoto(spades_dati, config):
	"""Controlla che l'output di Spades non sia vuoto"""
	assert len(spades_dati.contigs) > 0

@pytest.mark.spades
def test_correct_data_in_contig_name(spades_dati):
	"""Controlla che nel nome della contig in output sia presente un valore di coverage"""
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
	"""Controlla che tutte le lettere nella contig siano riconducibili alle basi del DNA"""
	A = True
	for r in spades_dati.contigs:
		for char in spades_dati.contigs[r]:
			A = A and (char == 'A' or char == 'T' or char == 'C' or char == 'G' or char == 'N' or char == 'a' or char == 't' or char == 'c' or char == 'g' or char == 'n')
		if not A:
			break
	assert A

#-----------------------------------------------------------------------------------------------------------------------------------
#       TEST FOR MAIN RULE
#-------------------------------------------------------------------------------------------------------------------------------------
@pytest.mark.main
def test_main_1(main_dati, config):
	"""Controlla che l'output di mani non sia vuoto"""
	assert len(main_dati.mapping) > 0

@pytest.mark.main
def test_coverage(main_dati):
	"""Controlla che nel nome della contig mappata sia presente un valore di coverage"""
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
	"""Controlla che il campo contrassegnato come aro_accession sia numerico"""
	A = True
	for r in main_dati.mapping:
		for f in main_dati.mapping[r]:
			A = A and main_dati.mapping[r][f]["ARO_accession"].isdigit()
	assert A
