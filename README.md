# INSTALLAZIONE DELLA PIPELINE

Per installare i tool impiegati in questa pipeline è sufficiente clonare questa git repo.
Unici requisiti fondamentali sono la precedente intallazione di SNAKEMAKE e ANACONDA.
La prima volta che si eesegue un allineamento, saranno installati tutti i programmi richiesti e sarà predisposto l'envoirment per il funzionamento di rgi.
L'installazione è locale e al termine sarà presente nella working directory una repository *rgi/*, contenente i riferimenti ai database Card e WildCard. 
Per eseguire esclusivamente l'installazione è possibile lanciare il comando

	snakemake rgi/card_annotation.log --core 10 --use-conda

# ESECUZIONIE DELLA PIPELINE

La pipeline è predisposta per allineare read libraries ai database sull'antibiotico resistenza CARD e wildcard.
I file di inpput devono essere due read libraries pair end e devono essere compresse (.gz).
Sono disponibili due modalità per completare questo compito: si può sfruttare l'algoritmo bowtie2 per allineare direttamente le reads ai resistomi di riferimento, oppure comporre un assembly delle reads utilizzando Spades e poi allineare le contigs ottenute ai resistomi (utilizzando Dimond o Blast)
Per eseguire interamente uno dei due rami della pipeline è necessario aggiungere al file *data_links.json* i link delle reads library che si intende analizzare, con etichetta uguale alla run accession del campione.

	"<run_accession>" : ["link_1", "link_2"]

Si osservi il file citato per un esempio, i link attualmente presenti sono riferiti ai dati di Munk et all. [2018].
Fatto questo si può richiedere a snakemake il file di output di uno dei due rami, in modo che vengano eseguite tutte le rules che sono necessarie ad ottenerlo
Se si vuole allineare un campione con bowtie2:
	
	snakemake <run_accesion>_bwt.allele_mapping_data.json --core <cores> --use-conda

Al termine dell'esecuzione saranno presenti i file:
- <run_accession>_bwt.allele_mapping_data.json
- <run_accession>_bwt.allele_mapping_data.txt
- <run_accession>_bwt.overall_mapping_stats.json

I primi due contengono entrambi il risultato dell'assembly e il terzo alcune statistiche relative al processo di allineamento.
Si faccia riferimento alla documentazione di rgi (https://github.com/arpcard/rgi#running-rgi-main-with-genome-or-assembly-dna-sequences) per il significato delle voci di questi file.
A seconda di quale sia il punto di partenza del processo, saranno presenti anche tutti i file generati dalle altre rules utilizzate per ottenere l'allineamento.
Si faccia riferimento alla struttura della piepline per capire quali.

Se si vuole allineare un campione componendo prima un assembly con Spades:

	snakemake <run accesion>_main.json --core <cores> --use-conda

Al termine dell'esecuzione saranno presenti i file:
- <run_accession>_main.json
- <run_accession>_main.txt

Entrambi contengono il risultato dell'assembly.
Si faccia riferimento alla documentazione di rgi [link](https://github.com/arpcard/rgi#running-rgi-main-with-genome-or-assembly-dna-sequences) per il significato delle voci di questi file.
A seconda di quale sia il punto di partenza del processo, saranno presenti anche tutti i file generati dalle altre rules utilizzate per ottenere l'allineamento.
Si faccia riferimento alla struttura della piepline per capire quali.

È possibile anche eseguire singolarmente ogni step della pipeline:
- Scaricare i file dai link riportati nel file *data_links.json*
Al termine saranno presenti nella working directory le read libraries scaricate

	snakemake <run_accession>_1.fastq.gz --core <core> --use-conda

- Eseguire il trimming pair end delle librerie. Per questa operazione sono utilizzate le sequenze degli adapter contenute nella directory *adapters/* e basta aggiungere in questa un altri .fa file per includere altri adapter nel trimming.
Al termine saranno presenti nella working directory le read libraries trimmate.

	snakemake <run_accession>_1_trimP.fastq.gz --core <core> --use-conda

- Eseguire l'assembly.
Al termine sarà presente nella working directory una directory chiamata *<run_accessio>_assembly/*, in essa sono contenuti i file risultanti dall'assembly

	snakemake <run_accession>_assembly/contigs.fasta --core <core> --use-conda

Al posto di <core> andrà inserito il numero di threads che si intende dedicare alla pipeline.
Nel file *config.yaml* si può specificare qunati threads dedicare all'elaborazione di un singolo campione (8 di default), qulora si decidesse di allineare più campioni contemporaneamente.
Nello stesso file sono stabiliti altri parametri per l'esecuzione della pipeline, la cui funzione è riassunta nello stesso, ma per informazioni più dettagliate si rimanda alla documentazione dei tools utilizzati:

- Ttimmomatic: [link](http://www.usadellab.org/cms/?page=trimmomatic)
- Spades: [link](https://cab.spbu.ru/files/release3.15.2/manual.html)
- Rgi: [link](https://github.com/arpcard/rgi#running-rgi-main-with-genome-or-assembly-dna-sequences)


## REPORT

Utilizzando il comando 

	snakemake report.out --report

Viene prodotto il report (*report.html*) riguardante tutti i campioni che sono contenuti nel file *data_links.json*
Questo file contiene, per ogni read library, trimmata e non:
- La distribuzione delle lunghezze delle reads
- La distribuzione della qualità media delle basi in ogni read
- Loccorrenza degli adapters in ogni read library
Questi grafici sono ottenuti da FastQC.
Si tenga presente che, se non si è già eseguito il trimming, questo verrà eseguito prima di ottenere l'output

# TESTING

Per compiere il testing della pipeline è utilizzato pytest.
Il testing viene eseguito utilizzando dei dati prodotti appositamente e che si trovano nella directory *test_dir*.
In essa si trovano anche gli script atti a generare detti dati:
- trimmomatic_test_reds.py
- bwt_test_subsampling.py
- get_spades_test_data.sh
- split_contigs.py
Nel file data_description.txt è spiegato come sono stati usati questi script per produrre i dati allegati e utilizzati per il testing.
Quando si abbia installato questo programma, è possibile testare tutta la pipeline con il comando (lanciato dalla directory principale).

	pytest test_dir/test_pipeline.py

Nel caso si voglia testare una sola rule della pipeline è possibile il sottoinseieme di tests dedicati attraverso i marker di pytest

	pytest test_dir/test_pipeline.py -m [marker]

I possibili valori di [marker] sono:
- trimmomatic
- bwt
- spades
- main

Si consiglia di utilizzare questi ultimi e di non eseguire tutti i tests assieme, potrebbe richiedere molto tempo.
