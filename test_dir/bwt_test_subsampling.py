import os
from numpy import random

seed = 3

def count_reads(filename):
	"""Return the number of the text based file passed as argument, devided by 4
	"""
	with open(filename, mode = 'r', encoding = "utf-8") as f:
		count = 0
		line = "pasta"
		while (line != ""):
			count +=1
			line = f.readline()
	return count/4

def get_bwt_test_data(test_file = "bwt_test_data", origin = "ERR2241634", how_many_reads = 100):
	"""This function don't returns anything, but create in the working ditectory two .fastq files.
	Those contains reads from a random (same indexes for both files) subsampling of the files with run accession "ogigin"
	Is expected that the "origin" read libraries are already trimmed, or at least that have a filename like trimmed reads
	The function works iven if "origin" files are compressed (.gz)
	"origin" file is suposed to be located in the external directory respect test_dir/, where should be also all other sample you are working on
	If origin = "example", this function sample from the file ../example_1_trimP.fastq(.gz) and ../example_2_trimP.fastq(.gz)
	Subsampled read libraries contains "how_many_reads" reads
	The run accession for resulting files will be "test_file" and will have a trimmed like name (example_1_trimP.fastq)
	"""
	random.seed(seed)
	os.system("gzip -dk ../" + origin + "_1_trimP.fastq.gz")
	os.system("gzip -dk ../" + origin + "_2_trimP.fastq.gz")
	selected_steps = random.randint(0, 100, how_many_reads)
	l = count_reads(origin + "_1_trimP.fastq.gz")
	while sum(selected_steps)*4 + how_many_reads) > l:
		selected_steps = random.randint(0, 100, how_many_reads)
	with open(test_file + "_1_trimP.fastq", mode = 'w', encoding = "utf-8") as f_1, open(test_file + "_2_trimP.fastq", mode = 'w', encoding = "utf-8") as f_2, open("../" + origin + "_1_trimP.fastq", mode = 'r', encoding = "utf-8") as orifin_1, open("../" + origin + "_2_trimP.fastq", mode = 'r', encoding = "utf-8") as origin_2:
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

#	COMANDO CONSIGLIATO

#	get_bwt_test_data()
