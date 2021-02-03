# TANDEMNEO: Analysis of neofunctionalization following gene tandem duplication in vertebrate evolution

## 1. Crea un ambiente virtuale con i pacchetti necessari (specificati nel file 'spec-file.txt')
###### nella Shell digita (necessario solo al primo utilizzo):

`$ conda create --name TandemNeo --file spec-file.txt`

## 2. Modifica il file config.py in base ai parametri di interesse

``` python
threads = 8 # define available threads for multiprocessing processes 
seed    = 2020 # define a seed for random sampling

########## 1. Species list ########## 2. Download

have_list         = False # True if you already have a species list, None if you want the script to write one 
species_file_name = 'species_list.csv' # Species list follows the syntax: [['class', 'order', 'genre', 'specie', 'publications', 'taxid', 'assembly'], ['class2', 'order2', etc...]] # write None as seat holder if you miss any field 
mg = 1 # define max number of genres, None if you don't want a limit 
mo = 12 # define max number of orders, None if you don't want a limit rf = ['Mammalia', 'Sauropsida', 'Actinopteri'] # define a reference classes list, None if you don't want to choose one 
rs = ['Homo_sapiens', 'Gallus_gallus', 'Danio_rerio'] # define a reference species list

########## 3. Main isoforms and BlastP

num_threads     = 20 # define blast params 
evalue          = '10e-6' 
max_hsps        = 1 
max_target_seqs = 5

########## 4. Duplications 
########## 5. Orthology

dups = ['tandem', 'divergent', 'convergent'] # define duplication classes

########## 6. Database 
########## 7. FASTAs for alignments 
########## 8. Alignments

matrix      = 'BLOSUM62' # define matrix that will be used whithin alignment phases 
srefs       = ['ENSP0', 'ENSP0'] # define reference specie ENSEMBL tags that will be used for extract alignment positions 
th_coverage = 60 # define minimum accepted coverage between orthogroups sequence lenghts

########## 9. Features

features_ref_specie = 'Homo_sapiens' # define reference specie for extract Uniprot 
features pos_range  = 1 # define max distance between alignment residue positions and feature positions 
alignment_threshold = 1 # define minumum accepted position score

######### 10. 3D Clustering

clust_ref_specie  = 'Homo_sapiens' 
clust_min_samples = 3 # minimum cluster dimension (?) 
clust_eps         = 10 # Angstrom -------
```

## 3. Tre modi di utilizzo del programma
### - Jupyter-Notebook:
###### Aggiunta dell'ambiente virtuale in jupyter (necessario solo al primo utilizzo)

`$ python -m ipykernel install --user --name=TandemNeo`

###### attivare prima la cella per l'import dei moduli e del config, poi attivare le celle di interesse

### - Linea di comando:
###### attivare l'ambiente virtuale (necessario ogni volta che si desidera utilizzare il programma)

`$ conda activate TandemNeo`

###### lasciare il sys.argv[1] == None per far andare tutto il programma

`$ python3 TandemNeo`

###### oppure specificare fasi in particolare come lista:

`$ python3 TandemNeo species_list`

```
usage: TandemNeo [-h] [-t THREADS] [-s SEED] [-l HAVE_SPECIES_LIST]
                 [-ll SPECIES_LIST] [-mg MAX_GENRES] [-mo MAX_ORDERS]
                 [-rf REF_CLASSES] [-rs REF_SPECIES] [-bt BLASTP_THREADS]
                 [-ev BLASTP_EVALUE] [-mh BLASTP_MAX_HSPS]
                 [-mt BLASTP_MAX_TARGET_SEQS] [-dc DUPLICATION_CLASSES]
                 [-m MATRIX] [-ta TAGS] [-tc THRESHOLD_COVERAGE]
                 [-fr FEATURES_REF_SPECIE] [-pr POS_RANGE]
                 [-at ALIGNMENT_THRESHOLD] [-cr CLUST_REF_SPECIE]
                 [-ms CLUST_MIN_SAMPLES] [-ep CLUST_EPS]
                 [species_list] [download] [main_isoforms] [duplications]
                 [orthology] [database] [fa_for_alignments] [alignments]
                 [features] [clustering]

positional arguments:
  species_list          Generate a species list file
  download              Download FASTA and GTF from Ensembl ftp servers
  main_isoforms         Looking for principal isoforms through Appris database
  duplications          Find duplicated genes along the genome
  orthology             Looking for orthology informations in Compara database
  database              Generate a database containing accessions found
                        informations for speed up next processes
  fa_for_alignments     Generate FASTA file containing orthology columns to
                        align
  alignments            Generate aligned files (FASTA format)
  features              Looking for features through the Uniprot database
  clustering

optional arguments:
  -h, --help            show this help message and exit
  -t THREADS, --threads THREADS
                        Define available threads for multiprocessing processes
  -s SEED, --seed SEED  Define a seed for random sampling processes
  -l HAVE_SPECIES_LIST, --have_species_list HAVE_SPECIES_LIST
                        Specify (True or False) if you already have a species
                        list
  -ll SPECIES_LIST, --species_list SPECIES_LIST
                        Species list follows the syntax: [[class, order,
                        genus, specie, publications, taxid, assembly],
                        [class2, order2, ...]]. Write None as seat holder if
                        you miss any field
  -mg MAX_GENRES, --max_genres MAX_GENRES
                        Define max number of genres, None if you don't want a
                        limit
  -mo MAX_ORDERS, --max_orders MAX_ORDERS
                        Define max number of orders, None if you don't want a
                        limit
  -rf REF_CLASSES, --ref_classes REF_CLASSES
                        Define a reference classes list (python list format),
                        None if you don't want to choose one
  -rs REF_SPECIES, --ref_species REF_SPECIES
                        Define a reference species list (python list format)
  -bt BLASTP_THREADS, --blastp_threads BLASTP_THREADS
                        Define blast params
  -ev BLASTP_EVALUE, --blastp_evalue BLASTP_EVALUE
                        Define blast params
  -mh BLASTP_MAX_HSPS, --blastp_max_hsps BLASTP_MAX_HSPS
                        Define blast params
  -mt BLASTP_MAX_TARGET_SEQS, --blastp_max_target_seqs BLASTP_MAX_TARGET_SEQS
                        Define blast params
  -dc DUPLICATION_CLASSES, --duplication_classes DUPLICATION_CLASSES
                        Define duplication classes list (python list format)
  -m MATRIX, --matrix MATRIX
                        Define matrix that will be used whithin alignment
                        phases
  -ta TAGS, --tags TAGS
                        Define reference specie ENSEMBL tags that will be used
                        for extract alignment positions
  -tc THRESHOLD_COVERAGE, --threshold_coverage THRESHOLD_COVERAGE
                        Define minimum accepted coverage between orthogroups
                        sequence lenghts
  -fr FEATURES_REF_SPECIE, --features_ref_specie FEATURES_REF_SPECIE
                        Define reference specie for extract Uniprot features
  -pr POS_RANGE, --pos_range POS_RANGE
                        Define max distance between alignment residue
                        positions and feature positions
  -at ALIGNMENT_THRESHOLD, --alignment_threshold ALIGNMENT_THRESHOLD
                        Define minumum accepted position score
  -cr CLUST_REF_SPECIE, --clust_ref_specie CLUST_REF_SPECIE
  -ms CLUST_MIN_SAMPLES, --clust_min_samples CLUST_MIN_SAMPLES
  -ep CLUST_EPS, --clust_eps CLUST_EPS
  ```

### - Installando il pacchetto:






# script 1 genera file input contenente lista specie "filtrata", genera albero cartelle

    # 

# script 2 fa il download

    # 
    
# script 3 isoforme principali e blastp

    # 
    
# script 4 duplicazioni

    #

# script 5 ortologia compara

    #
    
# script 6 database

    #
    
# script 7 FASTAS allineamento

    #
    
# script 8 Allineamenti

    # 
    
# script 9 Features (UNIPROT)