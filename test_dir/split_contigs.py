def split_contigs(quante, inputname = "original_contigs.fasta"):
	"""Produce un fasta file che contiene le prime "quante" contigs del file "inputname"
	Il nome del file creato è contigs.fasta
	"""
	f = open(filename)
	output = open("contigs.fasta", 'w')
	line = f.readline()
	for n in range(quante):
		output.write(line)
		while True:
			line = f.readline()
			if line[0] == '>':
				break
			output.write(line)

#	COMANDO CONSIGLIATO
#	split_contigs(30)
