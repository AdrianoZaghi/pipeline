Allo stato attuale, la pipeline è predisposta per completare un alignment ai database CARD e wildcard attraverso bowtie2.
Per eseguirla è necessario aggiungere al file json data_links i link delle reads library che si intende analizzare, con etichetta uguale alla run accession del campione.
Dopo di che si può eseguire il comando
	
	snakemake {run accesion}_bwt.allele_mapping_data.json --core 10 --use-conda

Per il funzionamento della pipeline dovrebbe bastare uno snakemake envoirement ed aver installato anaconda.
Questa infatti è predisposta per installare automaticamente tutte le dipendenze richieste.

Utilizzando il comando 

	snakemake report.out --report

Viene prodotto il report riguardante tutti i campioni che sono contenuti nel file "data_links.json"
Questo report allo stato attuale contiene solo alcune voci salienti dell analisi di FastQC sulle read libraries prima e dopo il trimming
