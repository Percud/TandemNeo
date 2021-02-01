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

## 3. Due modi di utilizzo del programma
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

`$ python3 TandemNeo ['species_list', 'download']`

###### Le fasi possono essere: 'species_list', 'download', 'main_isoforms', 'duplications', 'orthology', 'database', 'fa_for_alignments', 'alignments', 'features', 'clustering'







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