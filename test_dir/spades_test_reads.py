import math
import random

def quality_casting(i):
	if i >= 1:
#		print("qualità troppo alta, clipping effettuato")
		return "~"
	if i <= 0:
#		print("qualità troppo bassa, clipping effettuato")
		return "!"
	letters = """!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"""
	return letters[int(i*94)]

#SEQUENZE POSSIBILI

def random_seq(len):
	seq = ""
	for i in range(len):
		seq = seq + random.choice(["A", "T", "C", "G"])
	return seq

def banal_seq(len):
	seq = ""
	for i in range(len):
		if i < len/2:
			seq = seq + "A"
		else:
			seq = seq + "T"
	return seq

#QUALITÀ DISPONIBILI

def sinusoidal_qual(len, coefs):
	seq = ""
	norm = 0
	for i in coefs:
		norm += i[0]
	for i in range(len):
		val = 0
		for c in coefs:
			val += c[0] * math.sin(i * c[1] * math.pi/len)
		seq = seq + quality_casting(val/norm)
	return seq

def const_qual(len, coef):
	r = ""
	for i in range(len):
		r = r + quality_casting(coef)
	return r

def write_read(name, len, seq, qual, parameters_of_qual):
	read =[]
	read.append("@" + name)
	read.append(seq(len))
	read.append("+")
	read.append(qual(len, parameters_of_qual))
	return read

def print_data(reads, file_name = "spades_test_data", append_or_overwrite = "w"):
	f_1 = open(file_name + "_1_trimP.fastq", append_or_overwrite)
	f_2 = open(file_name + "_2_trimP.fastq", append_or_overwrite)
	if not (len(reads[0]) == len(reads[1])):
		print("errore_1")
		return
	for i in reads[0]:
		print(i[0], file = f_1)
		print(i[1], file = f_1)
		print(i[2], file = f_1)
		print(i[3], file = f_1)
#		print(i)
	for i in reads[1]:
		print(i[0], file = f_2)
		print(i[1], file = f_2)
		print(i[2], file = f_2)
		print(i[3], file = f_2)
#		print(i)

#scrivere delle reads con una qualità sempre che ha due sinusoidi sommate (sempre le stesse, con gli stessi coefficienti) e delle basi a caso, in modo che sovrasctiva l'output precedente

def get_data_1(quante = 100, len = 150, parameters = [[3,4], [5,6]]):
	data_0 = []
	data_1 = []
	for i in range(quante):
		data_0.append(write_read("get_data_1_" + str(i), len, random_seq, sinusoidal_qual, parameters))
		data_1.append(write_read("get_data_2_" + str(i), len, random_seq, sinusoidal_qual, parameters)) 
	return [data_0, data_1]

def get_data_2(quante = 100, len = 150, parameters = 1):
	data_0 = []
	data_1 = []
	for i in range(quante):
		data_0.append(write_read("get_data_1_" + str(i), len, banal_seq, const_qual, parameters))
		data_1.append(write_read("get_data_2_" + str(i), len, banal_seq, const_qual, parameters))
	return [data_0, data_1]

print_data(get_data_2())

