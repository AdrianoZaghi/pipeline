def split_contigs(how_many, inputname = "original_contigs.fasta"):
	"""This function creates a .fasta file containing the first "how_maiy" contigs which are found in the "inputname" file.
	Parameters of this function are the number of contigs required and the filename of the .fasta file from which the contigs are taken from
	This function don't return anything, but creates in the working directory a file named contig.fasta
	"""
	with open(inputname, mode = 'r', encoding = "utf-8") as input, open("contigs.fasta", mode = 'w', encoding = "utf-8") as output:
		line = input.readline()
		for n in range(how_many):
			output.write(line)
			while True:
				line = input.readline()
				if line[0] == '>':
					break
				output.write(line)

#	SUGGESTED COMAND

#	split_contigs(30)
