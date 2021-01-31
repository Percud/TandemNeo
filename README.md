# TANDEMNEO: Analysis of neofunctionalization following gene tandem duplication in vertebrate evolution

## 1. Crea un ambiente virtuale con i pacchetti necessari (specificati nel file 'spec-file.txt')
###### nella Shell digita (necessario solo al primo utilizzo):

    $ conda create --name TandemNeo --file spec-file.txt

## 2. Modifica il file config.py in base ai parametri di interesse

## 3. Due modi di utilizzo del programma
### - Jupyter-Notebook:
###### Aggiunta dell'ambiente virtuale in jupyter (necessario solo al primo utilizzo)

    $ python -m ipykernel install --user --name=TandemNeo

###### attivare prima la cella per l'import dei moduli e del config poi
###### attivare le celle di interesse

### - Linea di comando:
###### attivare l'ambiente virtuale (necessario ogni volta che si desidera utilizzare il programma)

    $ conda activate TandemNeo

###### lasciare il sys.argv[1] == None per far andare tutto il programma

    $ python3 TandemNeo

###### oppure specificare fasi in particolare come lista:

    $ python3 TandemNeo ['species_list', 'download']

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