#PIPELINE INSTALLING

In order to install all the tools used you should just clone this repo and run the pipeline for the first time. The only requirements are the previous installation of SNAKEAMAKE and ANACONDA (or MAMBA). During the installation will be set up the environment required for rgi to work. The installation is local and at the end of it a new *rgi/* repo will be found. It contains all the references to Card and WildCard  database. If you intend to execute the installation only, you can launch the command

	snakemake rgi/card_annotation.log --core 10 --use-conda

#PIPELINE SET UP

The pipeline is made to align read libraries to the antibiotic resistance databases Card and WildCard. Input files should be pair end read libraries compressed (.gz). File names should be provided in the form of *<run_accession>*_n.fastq.gz, where “n” should be 1 or 2 depending if the reads in the library are the forward or the backward once of the sample. To perform the alignment can be used two different methods:
- reads can be aligned directly to the references exploiting bowtie2
- otherwise can be composed an assembly of those reads using Spades and than the resulting contigs are aligned to the references (using Dimond or Blast, depending on the configuration of the pipeline)
To run entirely the pipeline (one of the two branches) you should add to the *data_links.json* file the links to the the online read libraries you intend to analyze, with the following sintax.

	"<run_accession>" : ["link_1", "link_2"]

At the beginning of the execution, data will be downloaded in the working directory. Their filename should match the standard required, but is surely worth a check. The links actually present in the *data_links.json* file refers to data from a study performed by Munk ef all. [2018]. In case the data ere already on your device and have a correct filename, the pipeline will skip the downloads, but is still useful to fill the “<run_accession>” part of the *data_links.json* file. This will enable the report.

#PIPELINE EXECUTION

After completing the setting up, the pipeline can be executed requesting to snakemake the output file of one fo the two branches. In this way are executed all the rules required to obtain it. In case you want to perform an alignment with bowtie2:

	snakemake <run_accesion>_bwt.allele_mapping_data.json --core <cores> --use-conda

At the end of the execution should be found, in the working directory, the following files:
- <run_accession>_bwt.allele_mapping_data.json
- <run_accession>_bwt.allele_mapping_data.txt
- <run_accession>_bwt.overall_mapping_stats.json

First two files contain the results of the alignment and the third contains some statistics about the alignment process, provided by rgi. You should refer to rgi documentation [link](https://github.com/arpcard/rgi#running-rgi-main-with-genome-or-assembly-dna-sequences) for the explanation of such output content. Depending from what the starting point of the analysis is, in the working directory will be found also all the files produced by the rules used in the alignment. Look at the picture of pipeline structure to figure out what they are.

If you want to align a sample composing the assembly and then using Dimond or Blast:

	snakemake <run accesion>_main.json --core "cores" --use-conda

At the end of the process will be found, in the working directory, the following output files:
- <run_accession>_main.json
- <run_accession>_main.txt

Both contain the results of the alignment. You should again refer to tgi documentation [link](https://github.com/arpcard/rgi#running-rgi-main-with-genome-or-assembly-dna-sequences) for the explanation of the content of those files. Depending from what the starting point of the analysis is, in the working directory will be found also all the files produced by the rules used in the alignment. Look at the picture of pipeline structure to figure out what they are.

It is possible to run each step of the pipeline singularly:

- Download the data from links in *data_links.json*. At the end you will find the data in the working directory.

	 snakemake <run_accession>_1.fastq.gz --core "cores" --use-conda

- Perform read library trimming trimming (pair end). For this task are used the adapter sequences contained in the *adaprers/* dir. To include new adapters in the trimming other *.fa* file should be added to this directory. At the end of this operation trimmed read libraries will be found in the working directory.

	snakemake <run_accession>_1_trimP.fastq.gz --core "cores" --use-conda

- Perform the assembly. At the end of this task in the working directory will be found *<run_accession>_assembly/*. In it are contained assembly results

	snakemake <run_accession>_assembly/contigs.fasta --core "cores" --use-conda

Is quite hard to estimate the computational resources required for each step of the pipeline, but all of them are quite expensive. Mainly the assembly. I suggest to execute reach step singularly, at the first usage.

In the place of *"cores"* should be put the number of threads dedicated to pipeline execution. In the *config.yaml* file you can specify how many threads dedicate to the processing of a single sample (8 default), in case you are aligning more then one sample at once. In the same file are sett some other parameters for pipeline execution. Their meaning is explained in that file, but for more information you should refer to the documentation of those specific tools.
- Ttimmomatic: [link](http://www.usadellab.org/cms/?page=trimmomatic)
- Spades: [link](https://cab.spbu.ru/files/release3.15.2/manual.html)
- Rgi: [link](https://github.com/arpcard/rgi#running-rgi-main-with-genome-or-assembly-dna-sequences)



#REPORT

Using the command:

	snakemake report.out --report

Is produced the report (*report.html*). It regards all the samples referenced in the file *data_links.json*. This report contains, for each read library, trimmed and non trimmed:
- length distribution of reads
- average quality distribution for the n-th base of all the reads
- occurrence of adapters in each read library

Those graphs are obtained from FastQC [link](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). Notice that, if the trimming is not already been executed, it will be done before producing the report.


#TESTING

To test the functionality of the pipeline is used *pytest*. Ad hoc data are used to perform tests. Those can be found in the *test_dir* directory.

- *trim_test_data* are generated using the script *trimmomatic_test_rendreads.py*. Those are .fastq file containing pair end reads randomly generated. Those have null superposition. In those reads are added the NextFlex-96 adapters.

- *bwt_test_data* is obtained subsampling the sample ERR2241634 (already trimmed). To do this have been used the script *bwt_test_subsampling.py*

- *reads_for_test_spades* have been created using ART_illumia. The command used is reported in *get_spades_test_data.sh* and require to have an *original_contig.fasta* file in *test_dir/*. Because of usual dimensions of .fasta file, original_contig.fasta is not included in this directory .Anyway whatever .fasta file can be used.

- In the directory *main_test_data/* can be found the output of the script *split_contigs.py*. This function compose, from the mentioned *original_contig.fasta*, a new .fasta file, containing the first n contigs of the original

In order perform the tests, you should run, in the working directory, the command

	pytest test_dir/test_pipeline.py

In case you want to test just a single rule of the pipeline is possible to specify the required subset of tests to perform, using pytest markers.

	pytest test_dir/test_pipeline.py -m [marker]

Possivle [marker] values are:
- trimmomatic
- bwt
- spades
- main


