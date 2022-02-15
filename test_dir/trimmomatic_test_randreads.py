from numpy import random
import os

def quality_casting(i):
	"""Convert numerichal values (between 0 and 1) in the quality codification of .fastq format
	The only argument is the numerichal value to convert and is returned the character resulting from the convertion
	In case i > 1 or i < 0, the returned dalue will be the higher ('~') or the lower ('!') quality character
	"""
	if i >= 1:
		return "~"
	if i <= 0:
		return "!"
	letters = """!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"""
	return letters[int(i*94)]

def const_adapter_end(i):
	"""Example of a distribution that can be used to produce the values of "adapter_end_position" for various reads
	"i" is suposed to be the number of the generated read
	"""
	return 33 + i * 0

def growing_quality(i, parameter):
	"""Example of a distribution of quality for a read sequence
	"i" is the number index of the base and parameter should be the number of generated reads
	"""
	return 0.4 + i * 0.01 - parameter * 0.001

def write_read(name, read_length, adapter, adapter_end_position, quality_distribution, parameter_of_quality_distribution = 0):
	"""This function returns "read", a 4 string list that contains the 4 fields that a .fastq file dedicates to a read:
	read[0]		Read name, it will be equal to the function argument "name"
	read[1]		Base sequence of the read. basis are randomly generated and are in number "read_length".
				"adapter" sequences will be overwritten in the final part of them. A number of "adapter_end_position" basis will vìbe overwritten, starting from the end of the read
				If "adapter_end_position" is <= 0 the adapter is not introduced
				If "adapter_end_position" is > len(adapter), it will be overwritten as "adapter_end_position" = len(adapter)
	read[2]		"+"
	read[3]		Sequence of base quality value. It will follow the "quality_distribution"
				The distribution passed as argumenmust have a parameter, so that it can vary from read to read when those are generated in serie. It is not necesssary that this parameter is actually used in the function.
				Values of the distribution are translated into corresponding quality characters as long as those are >= 0 and <= 1; otherwise the quality value is set to the lowest or the highest possible value
	"""
	read = []
	read.append("@" + name)
	sequence = []
	quality = ""
	for i in range(read_length):
		sequence.append(random.choice(["A", "T", "C", "G"]))
		quality = quality + quality_casting(quality_distribution(i, parameter_of_quality_distribution))
	if (adapter_end_position > len(adapter) - 1):
		adapter_end_position = len(adapter) - 1
	if (adapter_end_position >= 0):
		sequence[len(sequence)-adapter_end_position-1 : len(sequence)] = adapter[0 : adapter_end_position+1]
	read.append("".join(sequence))
	read.append("+")
	read.append(quality)
	return read


def get_trimmomatic_test_data(adapter_1, adapter_2, file_name = "trim_test_data", read_number = 300, read_length = 150, adapter_end_distribution = const_adapter_end, quality_distribution = growing_quality):
	"""This function don't return anything, but creates 2 files in the working directory called "file_name"_1.fastq and "file_name"_1.fastq
	Each of them contains "read_number" reads of length "read_length", composed of ransom bases.
	Such situation reproduces a pair end read database with null overlap, similar to [Munk et all. 2018] data
	The two "adapter" (_1 for the _1 file and _2 for the _2 file) sequences are then overwritten to part of the sequences, to one of the margin.
	"adapter_end_distribution" should be a function that specify a relation between the adapter end position and the read index.
		In this way, when generating iteratively the reads, those can have different adapter end position
	"quality_distribution" specify the distribution of reads base quality. This distribution can depend on an additional parameter (other than the read base number)
		This parameter is suposed to be the index of the read, so that the quality can change from a read to an other
	"""
	with open(file_name + "_1.fastq", mode = 'w', encoding = "utf-8") as f_1, open(file_name + "_2.fastq", mode = 'w', encoding = "utf-8") as f_2:
		for i in range(read_number):
			read = write_read("read " + str(i) + "_1", read_length, adapter_1, adapter_end_distribution(i), quality_distribution, i)
			print(read[0], file = f_1)
			print(read[1], file = f_1)
			print(read[2], file = f_1)
			print(read[3], file = f_1)
			opposit = write_read("read " + str(i) + "_2", read_length, adapter_2, adapter_end_distribution(i), quality_distribution, i)
			print(opposit[0], file = f_2)
			print(opposit[1], file = f_2)
			print(opposit[2], file = f_2)
			print(opposit[3], file = f_2)


def get_adapters(adapters_file_path = "../adapters/NEXTflex_96.fa", adapter_number = 4):
	"""This function returns a list of two strings. Those strings contain the 2 adapts indexed by the number "adapter_number": [adapter_1, adapter_2], where 1 is the adapter for forward reads and 2 is the adapter for reverse reads
	Those adapteres are taken from the file "adapters_file_path". This file should have all the lables of the adapter such that those should end with "<adapter_number>/1" "or <adapter_number>/2"
	"""
	with open(adapters_file_path, mode = 'r', encoding = "utf-8") as f:
		look_for = str(adapter_number) + "/1"
		linea = ""
		while linea[len(linea)-len(look_for)-1:len(linea)-1] != look_for:
			linea = f.readline()
		adapter_1 = f.readline()
		adapter_1 = adapter_1[0:len(adapter_1)-1]
		f.close()
	with open(adapters_file_path, mode = 'r', encoding = "utf-8") as f:
		look_for = str(adapter_number) + "/2"
		while linea[len(linea)-len(look_for)-1:len(linea)-1] != look_for:
			linea = f.readline()
		adapter_2 = f.readline()
		adapter_2 = adapter_2[0:len(adapter_2)-1]
	return [adapter_1, adapter_2]

#	COMANDO CONSIGLIATO

#	adapters = get_adapters()
#	get_trimmomatic_test_data(adapters[0], adapters[1])
