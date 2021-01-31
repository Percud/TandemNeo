import requests
import pandas as pd
import os
from multiprocessing import Pool
import natsort as ns
from itertools import product
import re
import sys

cwd = os.path.dirname(os.getcwd())
sys.path.append(cwd)

class ortho:
    
    def mergedups(rs, k):
        """return a dataframe containing all 
        given reference species gene informations"""

        merged = []
        for s in rs:
            df = pd.read_table(
                cwd + '/duplications/' + s + '_' + k + '.tsv', 
                header=None)
            merged.append(df)
        df = pd.concat(merged)
        df = df.reset_index(drop=True)
        df = df.reset_index()
        df[6] = df['index']
        df = df.drop(columns='index')

        """add a tag for each duplicated genes"""

        df[6] = df[6].astype(str)
        df[7] = df[0] + '-' + df[1] + '-' + df[6] + 'A' + ',' + df[2] + '-' + df[3] + '-' + df[6] + 'B'
        df = df[[4,5,6,7]]
        df = (df.set_index([4,5,6])
              .apply(lambda x: x.str.split(',').explode())
              .reset_index())  
        df[8] = df[7].apply(lambda x: x.split('-')[0])
        df[9] = df[7].apply(lambda x: x.split('-')[1])
        df[10] = df[7].apply(lambda x: x.split('-')[2])
        df = df[[8,9,4,5,10]]

        return df
    
        #ortho.mergedups(rs, 'tandem')

    def df(q, havespecielist):
        """return a df containing found 
        orthologues (COMPARA server)"""

        rest = "https://rest.ensembl.org/homology/id/"
        typ = "?type=orthologues;"
        args = "aligned=0;cigar_line=0"
        head = {"Content-Type": "application/json"}
        req = rest + q + typ + args
        try:
            r = requests.get(req, headers=head).json()
        except:
            r = None
        
        ndf = pd.DataFrame({'id': None, 
                          'protein_id': None, 
                          'perc_id': None, 
                          'species': None, 
                          'specie_query': None}, 
                        index=[0])
        
        if not 'error' in r:
            df = pd.DataFrame(r['data'][0]['homologies'])
            
            if 'source' in df.columns.tolist():
                df = df[['source', 'target']]
                dfsource = df['source'].apply(pd.Series).loc[0]
                dftarget = df['target'].apply(pd.Series)
                df = dftarget.append(dfsource)
                df = df.sort_values(['species', 'perc_id'], 
                                    ascending=[True,False])
                df = df.drop_duplicates('species')
                df['species'] = df['species'].apply(lambda x: x.capitalize())
                df['specie_query'] = dfsource['species'].capitalize()

                """if specified, filter the database based on a species list"""

                if havespecielist:
                    slist = pd.read_csv(cwd + '/species/species_list.csv')
                    slist = [l[3] + '_' + l[4] for l in slist.values.tolist()]

                    df = df[df['species'].isin(slist)]

                df = df[['id', 'protein_id', 'perc_id', 'species', 'specie_query']]

                return df
            
            else:
                return ndf
        else:
            return ndf

        #ortho.df(q, True)
        
    def brh(path_to_no_filter):
        """perform a best reciprocal hit 
        along all the orthologues found"""

        slist = pd.read_csv(cwd + '/species/species_list.csv')
        taxa = {}
        for l in slist.values.tolist():
            taxa.update({l[3]: [l[1], l[2]]})

        df = pd.read_table(path_to_no_filter, 
                           sep=';', 
                           names=[
                                'gene_id', 
                                'protein_id', 
                                'perc_id', 
                                'Species', 
                                'specie_query', 
                                'orthogroup'])

        df['orthogroup'] = pd.Categorical(
            df['orthogroup'], 
            ordered=True, 
            categories=ns.natsorted(
                df['orthogroup'].unique()))

        df = df.sort_values(['orthogroup', 'Species'], ascending=(True,False))
        df = df.reset_index(drop=True)

        df['Class'] = df['Species'].apply(lambda x: x if x == 'None' else taxa.get(x.split('_')[0])[0])
        df['Order'] = df['Species'].apply(lambda x: x if x == 'None' else taxa.get(x.split('_')[0])[1])

        df = df[df['protein_id'] != 'None']

        df['class+orthogroup'] = df['specie_query'].apply(lambda x: x[0])
        df['class+orthogroup'] = df['class+orthogroup'] + '.' + df['orthogroup'].astype(str)
        df['pair_number'] = df['orthogroup'].apply(lambda x: re.split('A|B', x)[0]).astype(int)
        df['pair_letter'] = df['orthogroup'].apply(lambda x: x[-1])

        # number compara hit for each query of each specie (one-to-many)
        # if 20 orthologues found of a Gallus gene, number from 1 to 20
        df = df.sort_values([
            'specie_query',
            'pair_number',
            'pair_letter',
            'Species',
            'perc_id'], ascending=(
            False, 
            True, 
            True, 
            True, 
            False))

        df['count'] = df.groupby([
            'Species', 
            'class+orthogroup']).cumcount()+1

        # drop duplicates if same percentage identity for 
        # each query of each specie to easier sorting
        df = df.drop_duplicates([
            'Species',
            'perc_id',
            'class+orthogroup'
        ])

        # sort based on global percentage identity, drop duplicates gene but the first (one-to-one),
        # drop duplicates same specie and same group ---> "best-hit"
        df = df.sort_values(
            'perc_id', ascending=False)
        df = df.drop_duplicates('gene_id')
        df = df.drop_duplicates([
            'Species',
            'class+orthogroup'])

        # save only first one-to-one hits
        # e.g. return an orthologue as first but it was the third ---> "reciprocal"
        df = df[df['count'] == 1]
        df = df.sort_values([
            'specie_query', 
            'pair_number', 
            'pair_letter'], ascending=(
            False, 
            True, 
            True))
        df['count2'] = df.groupby(
            'class+orthogroup', sort=False)[
            'class+orthogroup'].transform("count")

        counts = df[['class+orthogroup', 'count2']]
        counts = counts.drop_duplicates(
            'class+orthogroup').values.tolist()
        pmt = [] # select pairs more than 3 
        for x in range(len(counts)-1):
            x1 = re.split('A|B', counts[x][0])[0]
            x1 = x1.split('.')[1]
            x2 = re.split('A|B', counts[x+1][0])[0]
            x2 = x2.split('.')[1]
            if counts[x][1] and counts[x+1][1] > 3 and x1 == x2:
                pmt.append(counts[x][0])
                pmt.append(counts[x+1][0])

        df = df[df['class+orthogroup'].isin(pmt)]
        df = df.reset_index()
        df = df.drop(columns=['index'])
        df['specie_query'] = df.groupby([
            'pair_number', 'specie_query'], 
            sort=False, 
            as_index=False).ngroup()

        # pivot table
        df = df.pivot(index='Species', 
                      columns='orthogroup', 
                      values='protein_id')
        
        df = df.rename(columns=str).reset_index()
        df.columns = pd.Index(list(df.columns))
        df['Class'] = df['Species'].apply(lambda x: x if x == 'None' else taxa.get(x.split('_')[0])[0])
        df['Order'] = df['Species'].apply(lambda x: x if x == 'None' else taxa.get(x.split('_')[0])[1])
        df = df.set_index([
            'Class', 
            'Order', 
            'Species']).sort_index()
        df.to_csv(cwd + '/orthologues/' + path_to_no_filter.split('/')[-1].split('_no_')[0] + '.csv', sep=';')