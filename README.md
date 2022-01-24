# ESECUZIONIE DELLA PIPELINE

Allo stato attuale, la pipeline è predisposta per completare un alignment ai database CARD e wildcard attraverso bowtie2.
Per eseguirla è necessario aggiungere al file *data_links.json* i link delle reads library che si intende analizzare, con etichetta uguale alla run accession del campione.
Dopo di che si può eseguire il comando
	
	snakemake {run accesion}_bwt.allele_mapping_data.json --core [cores] --use-conda

Se si intende analizzare tutti i campioni a cui ci si riferisce in data_links si può eseguire il comando

	snakemake --core [cores] --use-conda

Si noti che al posto di [core] andrà inserito il numero di threads che si intende dedicare alla pipeline.
Nel file config.yaml si può specificare qunati threads dedicare all'elaborazione di un singolo campione (8 di default).
Nel caso ci si appresti a un analisi in parallelo di più campioni sarà necessario assegnare alla pipeline un numero di threads multiplo di questo parametro.

Per il funzionamento della pipeline dovrebbe bastare uno snakemake envoirement ed aver installato anaconda.
Questa infatti è predisposta per installare automaticamente tutte le dipendenze richieste.

## REPORT

Utilizzando il comando 

	snakemake report.out --report

Viene prodotto il report riguardante tutti i campioni che sono contenuti nel file *data_links.json*
Questo report allo stato attuale contiene solo alcune voci salienti dell analisi di FastQC sulle read libraries prima e dopo il trimming

# TESTING

Per compiere il testing della pipeline è utilizzato pytest.
Il testing viene eseguito utilizzando dei dati prodotti appositamente e che si trovano nella directory *test_dir*.
In essa si trovano anche gli script atti a generare detti dati.
Quando si abbia installato questo programma, è possibile testare tutta la pipeline con il comando (lanciato dalla directory principale)

	pytest

Nel caso si voglia testare una sola rule della pipeline è possibile specificarla attraverso i marker di pytest

	pytest -m [marker]

Allo stato attuale i possibili valori di [marker], oltre a quelli dei default, sono:
- trimmomatic
- bwt
- spades
