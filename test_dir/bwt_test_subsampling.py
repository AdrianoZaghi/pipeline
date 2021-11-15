import os
from numpy import random

seed = 3

def line_number_is_not_ok(filename, stop):
	"""Produce un report nel caso il subsampling si spinga troppo oltre ed esaurisca il fastq file
	"""
	f = open(filename)
	count = 0
	linea = "pasta"
	while (count < stop):
		count +=1
		if(f.readline()) == "":
			print("sono andato troppo in la nel file a cercare le reads")
			return True
	return False

def get_bwt_test_data(test_file = "bwt_test_data", origin = "ERR2241634", how_many_reads = 100):
	"""Produce due fastq file che sono il risultato di un random subsampling del fastq "orogon"
	Contengono how_many_reads reads ciascunno e sono prodotti in modo tale da essere paired eand
	La distribuzione con cui si fa il campionamento non è uniforme
	"""
	random.seed(seed)
	os.system("gzip -dk ../" + origin + "_1_trimP.fastq.gz")
	os.system("gzip -dk ../" + origin + "_2_trimP.fastq.gz")
	selected_steps = random.randint(0, 100, how_many_reads)
	while line_number_is_not_ok("../" + origin + "_1_trimP.fastq", sum(selected_steps)*4+how_many_reads):
		selected_steps = random.randint(0, 100, how_many_reads)
	f_1 = open(test_file + "_1_trimP.fastq", 'w')
	f_2 = open(test_file + "_2_trimP.fastq", 'w')
	origin_1 = open("../" + origin + "_1_trimP.fastq")
	origin_2 = open("../" + origin + "_2_trimP.fastq")
	for n in selected_steps:
		for _ in range(n*4):
			origin_1.readline()
			origin_2.readline()
		for _ in range(4):
			print(origin_1.readline(), file = f_1, end = '')
			print(origin_2.readline(), file = f_2, end = '')
	print("subsampling concluso")
	os.system("rm ../" + origin + "_1_trimP.fastq")
	os.system("rm ../" + origin + "_2_trimP.fastq")

get_bwt_test_data()
