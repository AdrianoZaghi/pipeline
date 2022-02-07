from numpy import random
import os

def quality_casting(i):
	"""Converte valori numerici (tra 0 e 1) della nel codice di qualità del formato fastq
	Nel caso i > 1 o i < 0 viene effettuato un clipping e si assegna il valore massimo o minimo di qualità, rispettivamente 
	"""
	if i >= 1:
#		print("qualità troppo alta, clipping effettuato")
		return "~"
	if i <= 0:
#		print("qualità troppo bassa, clipping effettuato")
		return "!"
	letters = """!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"""
	return letters[int(i*94)]

def const_adapter_start(i):
	return 33

def const_quality(i):
	return 0.4

def write_read(name, read_length, adapter, adapter_end_position, quality_distribution, parameter_of_quality_distribution = 0):
	"""Returna una lista di 4 stringhe che contengono:
		Nome della read, che sarà uguale a "name"
		Sequenza di basi della read: le basi sono estratte a random, è lunga "read_length" e ci viene sovrascritta una la sequenza dell'"adapter"
			La sovrascrittura dell'adapter comincia da una base distante "adapter_end_position" dalla fine della read
			"adapter_end_position" può essere anche <=0: in questo caso l'adapter non sarà introdotto
			"adapter_end_position" può essere anche > len(adapter), in questo caso l'adapter sarà posizionato come se "adapter_end_position" = len(adapter)
		+
		Sequenza dei valori di qualità delle basi: segue la "quality_distribution"
			questa quality distribution deve avere valori tra 0 e 1 per x da 0 a "read_length", altrimenti viene fatto il clipping
			queste distribuzioni sono munite di un parametro per poter essere diversificata quando le reads sono prodotte in serie
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


def get_trimmomatic_test_data(adapter_1, adapter_2, file_name = "trim_test_data", read_number = 300, read_length = 150, adapter_end_distribution = const_adapter_start, quality_distribution = const_quality):
	"""Produce due fastq file pair end dal nome "file_name_1.fastq" e "file_name_2.fastq".
	Ciascuno contiene "read_number" reads di lunghezza "read_length" con basi random, dunque come se l'overlap tra le due sequenze pair ended fosse nullo (caso di Munk)
	Si introduce in queste reads anche la sequenze di degli "adapter" specificate in argomento (_1 e _2 nei rispettivi file)
		Queste sequenze sono sovrascritte alle basi random, non per intero, ma solo come se fossero legate a un estremo della sequenza significativa e fossero in parte tagliate
		Adapter_end_distribution specifica una relazione tra l'indice della read e la posizione di collegamento adapter-sequenza_significativa
	Va specificata in argomento anche una distribbuzione dei valori di qualità associati alle basi (quality_distribution)
		questa distribuzione, a seconda di come sia implementaata, può dipenderere dall' indice della read nel campione.
	"""
	f_1 = open(file_name + "_1.fastq", 'w')
	f_2 = open(file_name + "_2.fastq", 'w')
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
	"""cerca nel file dato in argomento la coppia di adapter selezionata da "adapter_number" nel file "adapters_file_path".
	si presuppone che la riga del nome delgli adapters termini con il loro numero, prima di avere il consueto "/1" o "/2", necessario a distinguere i due adapter pair ended
	"""
	f = open(adapters_file_path)
	look_for = str(adapter_number) + "/1"
	linea = ""
	while linea[len(linea)-len(look_for)-1:len(linea)-1] != look_for:
		linea = f.readline()
	adapter_1 = f.readline()
	adapter_1 = adapter_1[0:len(adapter_1)-1]
	look_for = str(adapter_number) + "/2"
	while linea[len(linea)-len(look_for)-1:len(linea)-1] != look_for:
		linea = f.readline()
	adapter_2 = f.readline()
	adapter_2 = adapter_2[0:len(adapter_2)-1]
	print("adapter_1 = " + adapter_1)
	print("adapter_2 = " + adapter_2)
	return [adapter_1, adapter_2]

#	COMANDO CONSIGLIATO

#	adapters = get_adapters()
#	get_trimmomatic_test_data(adapters[0], adapters[1])
