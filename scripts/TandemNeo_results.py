import os
import pandas as pd
import re
import gzip
from Bio import SeqIO
import itertools
from functools import reduce
import json
import numpy as np

cwd = os.path.dirname(os.getcwd())

class ort:
    
    def whichdup(acc):

        def read(kind):
            path = cwd + '/orthologues/' + kind + '.csv'
            return pd.read_table(path, sep=';')

        t = read('Tandem').isin([acc]).any()
        d = read('Divergent').isin([acc]).any()
        c = read('Convergent').isin([acc]).any()

        if any(t):
            return 'Tandem'
        elif any(d):
            return 'Divergent'
        elif any(c):
            return 'Convergent'
        else:
            return None

        #ort.whichdup('ENSGALP00000072350')
        
    def dataframe(acc):
        kind = ort.whichdup(acc)
        path = cwd + '/orthologues/' + kind + '.csv'
        df = pd.read_table(path, sep=';')
        return df.set_index(['Class', 'Order', 'Species'])
        
        #ort.dataframe('ENSGALP00000072350')

    def pairs(acc):
        df = ort.dataframe(acc)
        p = df.columns[df.isin([acc]).any()][0]
        p_num = re.split('A|B', p)[0]

        if 'B' in p:
            pairs = [p_num + 'A', p]
        else:
            pairs = [p, p_num + 'B']

        return pairs

        #ort.pairs('ENSGALP00000072350')

    def ortholist(acc, ref_classes, includenone):
        df = ort.dataframe(acc)
        df = df.loc[ref_classes]
        df = df[ort.pairs(acc)].values.tolist()

        if not includenone:
            dfA = [l[0] for l in df 
                   if str(l[0]) != 'nan']
            dfB = [l[1] for l in df 
                   if str(l[1]) != 'nan']
        else:
            dfA = [l[0] for l in df]
            dfB = [l[1] for l in df]
            
        return [dfA, dfB]

        #ort.ortholist('ENSGALP00000072350', 'Sauropsida', False)

    def orthodict(acc, ref_classes, includenone):
        df = ort.dataframe(acc)
        c = {}
        for line in ref_classes:
            df2 = df[ort.pairs(acc)].loc[line].values.tolist()

            if not includenone:
                dfA = [l[0] for l in df2 
                       if str(l[0]) != 'nan']
                dfB = [l[1] for l in df2 
                       if str(l[1]) != 'nan']
            else:
                dfA = [l[0] for l in df2]
                dfB = [l[1] for l in df2]
            c.update({line: [{'A': dfA}, {'B': dfB}]})
        return c
    
        #ort.orthodict('ENSGALP00000072350', rf, True)

    def isclass(acc, cla):
        df = ort.dataframe(acc)
        df = df[ort.pairs(acc)].loc[cla]
        numna = len(df.values.tolist())
        num = len(df.iloc[:,0].dropna().tolist())
        if num > numna * 0.15:
            return True
        else:
            return False

        #ort.isclass('ENSGALP00000072350', 'Sauropsida')
        
    def whereisdup(ac, rf, t):

        path = cwd + '/orthologues/' + ort.whichdup(ac) + '.csv'
        df = pd.read_table(path, sep=';')
        df = df.set_index(['Class', 'Order', 'Species'])

        m = df[ort.pairs(ac)].loc[rf[0]].count()
        s = df[ort.pairs(ac)].loc[rf[1]].count()
        a = df[ort.pairs(ac)].loc[rf[2]].count()

        tsa = len(df.loc[rf[0]].index.tolist())
        tsm = len(df.loc[rf[1]].index.tolist())
        tss = len(df.loc[rf[2]].index.tolist())

        if m.tolist()[0] > tsm*t and m.tolist()[1] > tsm*t:
            ma = rf[0]
        else:
            ma = None

        if s.tolist()[0] > tss*t and s.tolist()[1] > tss*t:
            sa = rf[1]
        else:
            sa = None

        if a.tolist()[0] > tsa*t and a.tolist()[1] > tsa*t:
            ac = rf[2]
        else:
            ac = None

        mat = {
                str(['Mammalia', None, None]): 'Only Mammalia',
                str([None, 'Sauropsida', None]): 'Only Sauropsida',
                str([None, None, 'Actinopteri']): 'Only Actinopteri',
                str(['Mammalia', 'Sauropsida', None]): 'Mammalia and Sauropsida',
                str(['Mammalia', None, 'Actinopteri']): 'Mammalia and Actinopteri',
                str([None, 'Sauropsida', 'Actinopteri']): 'Sauropsida and Actinopteri',
                str(['Mammalia', 'Sauropsida', 'Actinopteri']): 'Mammalia, Sauropsida and Actinopteri',
                str([None, None, None]): 'Have not accessions',
                }

        return mat.get(str([ma, sa, ac]))

        #ort.whereisdup('ENSGALP00000072350', rf, 0.15)


class stats():

    # species classes count
    def ts():
        path = 'species_list.txt'
        s = pd.read_table(path, sep=' ', header=None)
        df = pd.DataFrame(s.groupby(0).count()[1]).reset_index().T
        df.columns = df.iloc[0]
        return df[1:]

        # stats.ts()

    # all isoforms count
    def ai(s):
        a = []
        path = cwd + '/fa/' + s + '.fa.gz'
        fa = SeqIO.parse(gzip.open(path, 'rt'), 'fasta')
        for fasta in fa:
            if 'gene_biotype:protein_coding' in fasta.description:
                a.append(fasta.id)
        return len(a)

        # stats.ai('Homo_sapiens')

    # main isoforms count
    def mi(s):
        path = cwd + '/main/' + s + '.tsv'
        df = pd.read_table(path, sep='\t')
        count = df['Protein_id'].count().tolist()
        return count

        # stats.mi('Homo_sapiens')
        
    # all and main isoforms ratio
    def aim(rs):
        mi = [stats.mi(x) for x in rs]
        ai = [stats.ai(x) for x in rs]
        index=['main', 'all']
        df = pd.DataFrame([mi, ai], index=index, columns=rs).T
        df['% main'] = df['main'] * 100 / df['all']
        return df
    
        # stats.aim('Homo_sapiens')

    # duplicated genes count
    def tc(s, k):
        path = cwd + '/duplications/' + s + '_' + k + '.tsv'
        df = pd.read_table(path, sep='\t', header=None)
        df = df.values.tolist()
        return len(df)*2
    
        # stats.tc('Homo_sapiens', 'tandem')
        
    def tcdf(rs, du):
        tcdf = {
                rs[0]: {
                    du[0]: stats.tc(rs[0], du[0]), 
                    du[1]: stats.tc(rs[0], du[1]), 
                    du[2]: stats.tc(rs[0], du[2]),
                    },
                rs[1]: {
                    du[0]: stats.tc(rs[1], du[0]), 
                    du[1]: stats.tc(rs[1], du[1]), 
                    du[2]: stats.tc(rs[1], du[2]),
                    },
                rs[2]: {
                    du[0]: stats.tc(rs[2], du[0]), 
                    du[1]: stats.tc(rs[2], du[1]), 
                    du[2]: stats.tc(rs[2], du[2]),
                    },
        }

        return pd.DataFrame(tcdf).T
    
        #stats.tcdf(rs, du)
        
    # return numbers and percentage about main species
    def allinfo(df1, df2):
        df = pd.merge(df1, df2, on='index')
        df['% tandem'] = df['tandem'] * 100 / df['main']
        df['% convergent'] = df['convergent'] * 100 / df['main']
        df['% divergent'] = df['divergent'] * 100 / df['main']
        df = df.set_index('index')
        df.index.name = None
        
        return df
    
        # stats.allinfo(stats.aim(rs).reset_index(), stats.tcdf(rs, du).reset_index())

    # blast comparisons count
    def bcf(specie):
        path = cwd + '/blast_queries/' + specie + '_main_blast.txt'
        x = pd.read_table(path, comment='#', header=None)
        x = x[x[2] != 100]
        q = x.drop_duplicates(0)[0].count().tolist()
        s = x[1].count().tolist()
        return [q, s]
    
    def bc(rs):
        df = pd.DataFrame({
        rs[0]: {'queries': stats.bcf(rs[0])[0], 
                'subjects': stats.bcf(rs[0])[1]},
        rs[1]: {'queries': stats.bcf(rs[1])[0], 
                'subjects': stats.bcf(rs[1])[1]},
        rs[2]: {'queries': stats.bcf(rs[2])[0], 
                'subjects': stats.bcf(rs[2])[1]}})
        return df

        #stats.bc(rs)
        
    def wherdup(kind, rf, threshold):

        def isdup(kind, c):
            path = cwd + '/orthologues/' + kind + '.csv'
            df = pd.read_table(path, sep=';')
            df = df.set_index(['Class', 'Order', 'Species'])

            # total number of accessions
            it = itertools.chain.from_iterable(df.values.tolist())
            accessions = len([l for l in list(it) if not str(l) == 'nan'])

            # if both columns sum is > 20%
            df = df.applymap(lambda x: 1 if str(x) != 'nan' else 0).loc[c]
            df = df.reset_index().drop(columns=['Order', 'Species']).T
            ts = int(df.columns.tolist()[-1])
            df['sum'] = df.sum(axis=1)
            df['Isdup'] = df['sum'].apply(lambda x: 1 if x > ts*threshold else 0) 

            def pairwise(iterable):
                "s -> (s0, s1), (s2, s3), (s4, s5), ..."
                a = iter(iterable)
                return zip(a, a)
            lp = list(pairwise(df.index))

            x = pd.DataFrame([df.loc[list(l)]['Isdup'].tolist() for l in lp])
            x['Isdup_' + c] = x.sum(axis=1)  
            x = x['Isdup_' + c].reset_index()

            return x

            #stats.isdup('Convergent', 'Actinopteri')

        df_list = []
        for x in range(len(rf)):
            df_list.append(isdup(kind, rf[x]))

        msa = reduce(lambda df1,df2: pd.merge(df1,df2,on='index'), df_list)
        msa = msa.drop(columns='index')

        return msa
    
        # stats.wherdup('Tandem', rf, 0.15)

    def whichdup(kind, rf, threshold):

        df = stats.wherdup(kind, rf, threshold)
        df['m'] = df['Isdup_Mammalia'].apply(lambda x: 'Mammalia' if x == 2 else None)
        df['s'] = df['Isdup_Sauropsida'].apply(lambda x: 'Sauropsida' if x == 2 else None)
        df['a'] = df['Isdup_Actinopteri'].apply(lambda x: 'Actinopteri' if x == 2 else None)
        df = df[['m', 's', 'a']]

        mat = {
                    str(['Mammalia', None, None]): 'Only Mammalia',
                    str([None, 'Sauropsida', None]): 'Only Sauropsida',
                    str([None, None, 'Actinopteri']): 'Only Actinopteri',
                    str(['Mammalia', 'Sauropsida', None]): 'Mammalia and Sauropsida',
                    str(['Mammalia', None, 'Actinopteri']): 'Mammalia and Actinopteri',
                    str([None, 'Sauropsida', 'Actinopteri']): 'Sauropsida and Actinopteri',
                    str(['Mammalia', 'Sauropsida', 'Actinopteri']): 'Mammalia, Sauropsida and Actinopteri',
                    str([None, None, None]): 'Have not accessions',
        }

        df['msa'] = df.values.tolist()
        df['which'] = df['msa'].apply(lambda x: mat.get(str(x)))
        df = pd.DataFrame(df['which'].value_counts().to_dict(), index=[kind])


        def addsum(df, k):
            classes = [l for l in df.columns.tolist() if k in l]
            df['Total ' + k] = df[classes].sum().sum()
            return df

        df = addsum(df, 'Mammalia')
        df = addsum(df, 'Sauropsida')
        df = addsum(df, 'Actinopteri')
        df = df.T.sort_values(kind, ascending=False)

        return df
    
        # stats.whichdup('Convergent', rf, 0.15)
    
    # All duplication events
    def alldups(t, rf):
        dfc = stats.whichdup('Convergent', rf, t).reset_index()
        dft = stats.whichdup('Tandem', rf, t).reset_index()
        dfd = stats.whichdup('Divergent', rf, t).reset_index()

        dftc = pd.merge(dft, dfc, on='index')
        df = pd.merge(dftc, dfd, on='index')
        df = df.set_index('index')
        df.index.name = None

        return df
    
        #stats.alldups(0.15, rf)
        
    def outofmean(kind, threshold):

        db = json.load(open(cwd + '/database.json'))

        path = cwd + '/orthologues/' + kind + '.csv'
        df = pd.read_table(path, sep=';')
        df = df.T.iloc[3:]
        df = df.applymap(lambda x: db.get(x)['sequence lenght'] if not db.get(x) == None else None)
        df['len_m'] = df.mean(axis=1)

        def pairwise(iterable):
            "s -> (s0, s1), (s2, s3), (s4, s5), ..."
            a = iter(iterable)
            return zip(a, a)

        lp = list(pairwise(df.index))
        x = pd.DataFrame([df.loc[list(l)]['len_m'].tolist() for l in lp]).rename(columns={0: 'A', 1: 'B'})
        x['AB'] = x.max(axis=1) - x.min(axis=1)
        x['Outofmean'] = x['AB'].apply(lambda x: True if x > threshold else False)
        x = x.groupby('Outofmean').count()['AB']

        return x

class uni:

    def all_scores(kind, threshold):

        path = cwd + '/alignments/' + kind + '_alignments/'
        file = 'alignments.json'
        db = json.load(open(path + file))

        scores = []
        for l in db:
            pos = db[l]['Positions']
            for l1 in pos:
                s = pos[l1]['Score']
                if s > threshold:
                    scores.append(s)

        return scores
    
        #uni.all_scores('Convergent', 0)
    
    def filtdf(kind, threshold):

        # scrivo i range di posizioni che mi interessano (in questo caso -1 e +1 sia per begin che per end)
        df = pd.read_table(cwd + '/uniprot/' + kind + '_filtered_features.csv', sep=',').drop(columns='Unnamed: 0').reset_index()
        df['begin'] = df['begin'].apply(lambda x: None if '~' in str(x) else x)
        df['end'] = df['end'].apply(lambda x: None if '~' in str(x) else x)
        df['b+1'] = df['begin'].apply(lambda x: str(int(x)+1) if not x == None else x)
        df['b-1'] = df['begin'].apply(lambda x: str(int(x)-1) if not x == None else x)
        df['e+1'] = df['end'].apply(lambda x: str(int(x)+1) if not x == None else x)
        df['e-1'] = df['end'].apply(lambda x: str(int(x)-1) if not x == None else x)
        df['b+1'], df['b-1'] = df['b+1'].astype(float), df['b-1'].astype(float)
        df['e+1'], df['e-1'] = df['e+1'].astype(float), df['e-1'].astype(float)
        df['begin'], df['end'] = df['begin'].astype(float), df['end'].astype(float)

        uniprot_list = df[['ensembl_accession', 'index', 'begin', 'b+1', 'b-1', 'end', 'e+1', 'e-1']].values.tolist()

        positions = json.load(open(cwd + '/alignments/' + kind + '_alignments/alignments.json'))

        # confronto le posizioni trovate in uniprot con quelle trovate mediante allineamenti
        prova = []

        for group in positions:
            x = positions[group]
            acc_A, acc_B = x['Accessions'][0], x['Accessions'][1]

            for pos in x['Positions']:
                x2 = x['Positions'][pos]
                score = x['Positions'][pos]['Score']

                if score > threshold:
                    al_A, al_B = x2['alignment_A'], x2['alignment_B']
                    pos_A, pos_B = float(x2['Homo_position_A']), float(x2['Homo_position_B'])
                    res_A, res_B = x2['Homo_residue_A'], x2['homo_residue_B']

                    for line in uniprot_list:
                        accession = line[0]
                        index = line[1]
                        begin, begin_plus_one, begin_minus_one = line[2], line[3], line[4]
                        end, end_plus_one, end_minus_one = line[5], line[6], line[7]

                        if accession == acc_A:

                            if pos_A == begin or pos_A == begin_plus_one or pos_A == begin_minus_one or pos_A == end or pos_A == end_plus_one or pos_A == end_minus_one:
                                prova.append([acc_A, index, pos_A, res_A, res_B, al_A, al_B, score])
                            if pos_B == begin or pos_B == begin_plus_one or pos_B == begin_minus_one or pos_B == end or pos_B == end_plus_one or pos_B == end_minus_one:
                                prova.append([acc_A, index, pos_B, res_A, res_B, al_A, al_B, score])    

                        if accession == acc_B:
                            if pos_A == begin or pos_A == begin_plus_one or pos_A == begin_minus_one or pos_A == end or pos_A == end_plus_one or pos_A == end_minus_one:
                                prova.append([acc_B, index, pos_A, res_A, res_B, al_A, al_B, score])
                            if pos_B == begin or pos_B == begin_plus_one or pos_B == begin_minus_one or pos_B == end or pos_B == end_plus_one or pos_B == end_minus_one:
                                prova.append([acc_B, index, pos_B, res_A, res_B, al_A, al_B, score]) 

        positions_found = pd.DataFrame(prova, columns=['ensembl_accession', 'index', 'position_found', 'residue_A', 'residue_B', 'alignment_A', 'alignment_B', 'score'])

        final_df = pd.merge(df, positions_found, on=['index', 'ensembl_accession']).drop_duplicates()

        # recupero i full protein name di ensembl dal database (proteina in esame e partner)
        database = json.load(open('database.json'))
        orthos = pd.read_table(cwd + '/orthologues/' + kind + '.csv', sep=';')
        orthos = list(zip(orthos[orthos['Species'] == 'Homo_sapiens'].values.tolist()[0][3:], orthos[orthos['Species'] == 'Homo_sapiens'].columns.tolist()[3:]))
        orthos = pd.DataFrame(orthos, columns = ['ensembl_accession', 'orthogroup'])
        final_df = pd.merge(final_df, orthos, on=['ensembl_accession', 'orthogroup'])
        orthos['full_name'] = orthos['ensembl_accession'].apply(lambda x: database.get(x)[7] if not database.get(x) == None else x)

        final_df['full_name_ensembl_partner'] = final_df['orthogroup'].apply(lambda x: database.get(
            orthos[orthos['orthogroup'] == re.split('A|B', x)[0] + {'A': 'B', 'B': 'A'}.get(x[-1])].values.tolist()[0][0])[7] if not database.get(
            orthos[orthos['orthogroup'] == re.split('A|B', x)[0] + {'A': 'B', 'B': 'A'}.get(x[-1])].values.tolist()[0][0]) == None else None)

        final_df['full_name_ensembl'] = final_df['orthogroup'].apply(lambda x: database.get(
            orthos[orthos['orthogroup'] == x].values.tolist()[0][0])[7] if not database.get(
            orthos[orthos['orthogroup'] == x].values.tolist()[0][0]) == None else None)

        # ordino la tabella e la salvo in locale
        final_df = final_df.drop(columns=['index', 'b+1', 'b-1', 'e+1', 'e-1']).drop_duplicates()[['orthogroup', 'accession', 'ensembl_accession', 'gene', 'full_protein_name', 'full_name_ensembl', 'full_name_ensembl_partner', 'type', 'category', 'description', 'begin', 'end', 'position_found', 'score', 'residue_A', 'residue_B', 'alignment_A', 'alignment_B', 'EC_number']]
        final_df = final_df[~final_df['category'].str.contains('STRUCTURAL')]
        final_df = final_df[~final_df['category'].str.contains('TOPOLOGY')]

        return final_df

        #uni.filtdf('Convergent', 0)
        
    def types(kind, threshold):

        if threshold:
            df = uni.filtdf(kind, threshold)
        else:
            df = pd.read_table(cwd + '/uniprot/' + kind + '_features.csv', sep=',')
        types = df.groupby('type').count()['accession'].reset_index()
        types = types.rename(columns={'accession': kind})
        u_types = df.drop_duplicates(['type', 'accession'])
        u_types = u_types.groupby('type').count()['accession'].reset_index()
        u_types = u_types.rename(columns={'accession': 'Unique_' + kind.lower()})
        df = pd.merge(types, u_types, on='type')

        return df
    
        #uni.types('Convergent', 0)
        
    def allfeat(threshold):

        p = pd.merge(uni.types('Tandem', threshold), uni.types('Convergent', threshold), on='type', how='outer')
        df = pd.merge(p, uni.types('Divergent', threshold), on='type', how='outer').set_index('type')
        df.index.name = None

        return df

        #uni.allfeat(None)