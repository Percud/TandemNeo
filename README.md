# TANDEMNEO
#### Analysis of neofunctionalization following gene tandem duplication in vertebrate evolution

https://github.com/Percud/TandemNeo

# 1. Download:
Download dal server GitHub della repository del progetto 
``` sh
$ git clone https://github.com/Percud/TandemNeo
```

# 2. Ambiente virtuale:
Creazione (necessaria solo al primo utilizzo) e attivazione dell'ambiente di lavoro virtuale
#### 2.1 Conda virtual environments
Consigliato se non si hanno già installati blastp e clustalo, entrambi necessari ($ conda install blastp, $ conda install clustalo)
``` sh
$ conda create -n TandemNeo python=3.7
$ conda activate TandemNeo
(TandemNeo) ... $
```
oppure 
#### 2.2 Python virtual environments
``` sh
$ python3.7 -m venv TandemNeo
$ source TandemNeo/bin/activate
(TandemNeo) ... $
```

# 3. Installazione pacchetto:
Setup.py installerà il pacchetto TandemNeo con i suoi moduli. Gli altri pacchetti necessari saranno installati automaticamente. Spostarsi nella cartella in cui è presente il setup.py e digitare:
``` sh
(TandemNeo) ... TandemNeo $ python3 setup.py install
```
Verificare se l'installazione è avvenuta correttamente:
``` sh
(TandemNeo) ... TandemNeo $ tandemneo --help
```

# 4. Utilizzo:
Raccomandazione: 
#### 4.1 Jupyter-Notebook:
Aggiunta dell'ambiente virtuale in jupyter (necessario solo al primo utilizzo), previa installazione dell'ipykernel
```sh
(TandemNeo) ... TandemNeo $ pip install ipykernel
(TandemNeo) ... TandemNeo $ ipython kernel install --user --name=TandemNeo
```
Attivando tutte le celle per avviare tutte le fasi: 
```sh
run --> run all cells
```
Oppure per una fase alla volta, prima la cella con gli import dei moduli e del config:
```python
import os, sys, gzip, json
...
from config import *
```
Poi la cella per la scrittura delle cartelle:
```python
os.mkdir(cwd + '/appris')
...
os.mkdir(cwd + '/clustering')
```
Infine le celle di interesse fornendo i file nel formato corretto (vedi cartella esempio)

#### 4.2 Linea di comando:
Oppure mediante linea di comando, lasciando il primo argomento vuoto per far eseguire tutte le fasi in successione: 
``` sh
(TandemNeo) ... TandemNeo $ tandemneo
```
oppure scrivendo la fase come primo argomento:
``` sh
(TandemNeo) ... TandemNeo $ tandemneo download
```
le fasi possibili sono: 
``` sh
  species_list          Generate a species list file
  download              Download FASTA and GTF from Ensembl ftp servers
  main_isoforms         Looking for principal isoforms through Appris database
  duplications          Find duplicated genes along the genome
  orthology             Looking for orthology informations in Compara database
  database              Generate a database containing accessions found informations for speed up next processes
  fa_for_alignments     Generate FASTA file containing orthology columns to align
  alignments            Generate aligned files (FASTA format)
  features              Looking for features through the Uniprot database
  clustering
  ```

# 5. Parametri:
#### 5.1 Config:
I parametri possono essere modificati agendo nel file config.py:
``` python
threads = 8 # define available threads for multiprocessing processes 
seed = 2020 # define a seed for random sampling

########## 1. Species list 
########## 2. Download

have_list = False                                           # True if you already have a species list, None if you want the script to write one 
species_file_name = 'species_list.csv'                      # Species list follows the syntax: [['class', 'order', 'genre', 'specie', 'publications', 'taxid', 'assembly'], ['class2', 'order2', etc...]]                                                    
                                                            # write None as seat holder if you miss any field 
mg = 1                                                      # define max number of genres, None if you don't want a limit 
mo = 12                                                     # define max number of orders, None if you don't want a limit 
rf = ['Mammalia', 'Sauropsida', 'Actinopteri']              # define a reference classes list, None if you don't want to choose one 
rs = ['Homo_sapiens', 'Gallus_gallus', 'Danio_rerio']       # define a reference species list

########## 3. Main isoforms and BlastP

num_threads = 20    # define blast params 
evalue = '10e-6' 
max_hsps = 1 
max_target_seqs = 5

########## 4. Duplications 
########## 5. Orthology

dups = ['tandem', 'divergent', 'convergent'] # define duplication classes

########## 6. Database 
########## 7. FASTAs for alignments 
########## 8. Alignments

matrix = 'BLOSUM62'         # define matrix that will be used whithin alignment phases 
srefs = ['ENSP0', 'ENSP0']  # define reference specie ENSEMBL tags that will be used for extract alignment positions 
th_coverage = 60            # define minimum accepted coverage between orthogroups sequence lenghts

########## 9. Features

features_ref_specie = 'Homo_sapiens'    # define reference specie for extract Uniprot 
features pos_range = 1                  # define max distance between alignment residue positions and feature positions 
alignment_threshold = 1                 # define minumum accepted position score

######### 10. 3D Clustering

clust_ref_specie = 'Homo_sapiens' 
clust_min_samples = 3                   # minimum cluster dimension (?) 
clust_eps = 10 # Angstrom -------
```

#### 4.2 Flags (linea di comando):
Oppure mediante l'utilizzo di flag appositi:
``` sh
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
                        Define max number of genres, None if you don t want a
                        limit
  -mo MAX_ORDERS, --max_orders MAX_ORDERS
                        Define max number of orders, None if you don t want a
                        limit
  -rf REF_CLASSES, --ref_classes REF_CLASSES
                        Define a reference classes list (python list format),
                        None if you don t want to choose one
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
Ad esempio:
``` sh
(TandemNeo) ... TandemNeo $ tandemneo features -fr 'Homo_sapiens'
```

# 6. Risultati
Nel notebook TandemNeo.ipynb è possibile trovare...

# 7. Extra
E' possibile utilizzare i moduli anche al di fuori del progetto, ad esempio:
``` sh
(base) ... TandemNeo $ conda activate TandemNeo
(TandemNeo) ... TandemNeo $ python
>>> from tandemneo.specieinfo import specieinfo
>>> h = specieinfo('Homo_sapiens')

>>> h.allinfo
['Mammalia', 'Primates', 'Homo', 'sapiens', '19139154', '9606', 'GRCh38']

>>> h.taxlist
      TaxId        ScientificName          Rank
0    131567    cellular organisms       no rank
1      2759             Eukaryota  superkingdom
2     33154          Opisthokonta         clade
3     33208               Metazoa       kingdom
4      6072             Eumetazoa         clade
5     33213             Bilateria         clade
6     33511         Deuterostomia         clade
7      7711              Chordata        phylum
8     89593              Craniata     subphylum
9      7742            Vertebrata         clade
10     7776         Gnathostomata         clade
11   117570            Teleostomi         clade
12   117571          Euteleostomi         clade
13     8287         Sarcopterygii    superclass
14  1338369  Dipnotetrapodomorpha         clade
15    32523             Tetrapoda         clade
16    32524               Amniota         clade
17    40674              Mammalia         class
18    32525                Theria         clade
19     9347              Eutheria         clade
20  1437010         Boreoeutheria         clade
21   314146      Euarchontoglires    superorder
22     9443              Primates         order
23   376913           Haplorrhini      suborder
24   314293           Simiiformes    infraorder
25     9526            Catarrhini     parvorder
26   314295            Hominoidea   superfamily
27     9604             Hominidae        family
28   207598             Homininae     subfamily
29     9605                  Homo         genus
```


# TANDEMNEO: le fasi
### 0. Assunto iniziale
Testo qui
### 1. Lista delle specie
Testo qui
### 2. Download
Testo qui
### 3. Definizione isoforme principali e BLASTP
Testo qui
### 4. Ricerca duplicazioni
Testo qui
### 5. Ricerca ortologia
Testo qui
### 6. Scrittura di un database
Testo qui
### 7. Scrittura dei FASTA per gli allineamenti
Testo qui
### 8. Allineamenti
Testo qui
### 9. Ricerca delle features
Testo qui
### 10. Analisi di clustering 3D
Testo qui
### 11. Analisi dei risultati
Testo qui

