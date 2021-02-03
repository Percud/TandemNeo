import os
import pandas as pd

cwd = os.path.dirname(os.getcwd())

class duplications:
    
    def blastfile(s):
        """return blast output file as dataframe"""

        path = cwd + '/blast_queries/'
        file = path + s + '_main_blast.txt'

        df = pd.read_table(file, 
                          comment='#', 
                          header=None, 
                          names=['Protein_id', 
                                 'Subject', 
                                 'Identity', 
                                 'A.lenght', 
                                 'Mismatches', 
                                 'Gap', 
                                 'Q.start', 
                                 'Q.end', 
                                 'S.start', 
                                 'S.end', 
                                 'Evalue', 
                                 'Bits'])

        return df
    
        #duplications.blastfile(s)
    
    def isbrh(s):
        """check which accessions are best reciprocal hit"""

        df = duplications.blastfile(s)

        df = df[df['Protein_id'] != df['Subject']]
        df = df.drop_duplicates('Protein_id')

        df['Isbrh'] = (df['Protein_id'] + df['Subject']).isin(
                       df['Subject'] + df['Protein_id'])
        
        df = df[df['Isbrh'] == True]
        df = df.drop(columns='Isbrh')

        return df
    
        #duplications.isbrh(s)
        
    def duplicationsdf(s):
        """filter for main isoforms and 
        return a dataframe"""

        brh = duplications.isbrh(s)

        main = cwd + '/main/' + s + '.tsv'
        main = pd.read_table(main)

        df = pd.merge(main, brh, on='Protein_id')
        df = df[(df['Protein_id'] == df['Subject'].shift(-1)) | 
                (df['Subject'] == df['Protein_id'].shift(+1))]
        
        return df

        #duplications.duplicationsdf(s)
        
    def isduplication(s, a, df):
        """looks for duplications and return 
        a dataframe. Verified if genes are 
        consecutive, have same direction and 
        are not overlapping"""

        if not any(df):
            df = duplications.duplicationsdf(s)[['Protein_id', 
                                                 'Subject', 
                                                 'Strand', 
                                                 'Start', 
                                                 'Stop']]
        else:
            df = df[['Protein_id',
                     'Subject',
                     'Strand', 
                     'Start', 
                     'Stop',
                     'Gene_id',
                     'Evalue', 
                     'Identity']]

        i = df[df['Protein_id'] == a].index[0]
        if i+1 in df.index.tolist():
            l = df.loc[[i, i+1]].values.tolist()
            que, sub = l[0], l[1]
            pid1, pid2 = que[0], sub[1]
            gen1, gen2 = que[5], sub[5]
            str1, str2 = que[2], sub[2]
            beg1, beg2 = que[3], sub[3]
            end1, end2 = que[4], sub[4]
            eva1, ide1 = que[6], que[7]
            

            if pid1 == pid2:
                if str1 == str2:
                    if str1 == '+' and end1 < beg2 or str1 == '-' and beg1 < end2:
                        return [gen1, pid1, gen2, que[1], ide1, eva1, 'tandem']
                if str1 != str2:
                    if str1 == '+' and end1 < end2:
                        return [gen1, pid1, gen2, que[1], ide1, eva1, 'convergent']
                    if str1 == '-' and beg1 < beg2:
                        return [gen1, pid1, gen2, que[1], ide1, eva1, 'divergent']
        else:
            return None
        
        #duplications.isduplication(s, accession, df (or None))
        
    def duplist(s, k):
        """return a duplications list for a given 
        specie and kind of duplications"""

        df = duplications.duplicationsdf(s)

        kl = []
        for a in df['Protein_id'].tolist():
            dup = duplications.isduplication(s, a, df)
            if dup and dup[6] == k:
                kl.append(dup[:6])

        return kl

        #duplications.duplist('Homo_sapiens', 'convergent')
        
    def tocsv(s, k):
        """store duplications containing 
        dataframe in a tsv format file"""

        dups = duplications.duplist(s, k)
        for l in dups:
            path = cwd + '/duplications/'
            file = path + s + '_' + k + '.tsv'
            print(*l, sep='\t', file=open(file, 'a'))

        #duplications.tocsv('Homo_sapiens', 'tandem')