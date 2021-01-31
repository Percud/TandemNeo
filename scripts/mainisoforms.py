from specieinfo import specieinfo as si
import pandas as pd
import urllib.request
import gzip
import os
import re
import natsort as ns
from Bio import SeqIO

cwd = os.path.dirname(os.getcwd())

class mainisoforms:

    def df(s):
        """retrieve a main isoforms containing txt 
        file from appris servers for a given specie 
        and return a dataframe"""

        a = si.assembly(s)
        cwd = 'pub/current_release'
        appris = 'http://apprisws.bioinfo.cnio.es/' + cwd
        file = '/appris_data.principal.txt'
        url = appris + '/datafiles/' + s.lower() + '/' + a + file

        df = pd.read_table(urllib.request.urlopen(url), header=None)
        df = df[df[4].str.contains('PRINCIPAL')]
        df = df.sort_values([1, 4]).drop_duplicates(1)

        return df

        #mainisoforms.df('Homo_sapiens')

    def writeappris(s):
        """store the main.isoforms dataframe 
        for a given specie as txt file"""

        df = mainisoforms.df(s)
        for l in df.values.tolist():
            f = cwd + '/appris/' + s + '_principal_appris.txt'
            print(*l, sep='\t', file=open(f, 'a'))

        #mainisoforms.writeappris('Homo_sapiens')

    def mlist(s):
        """return main isoforms for 
        a given specie as list"""

        df = mainisoforms.df(s)[2]
        return df.values.tolist()

        #mainisoforms.mlist('Homo_sapiens')
        
    def gtf(s):
        """filter the given specie gtf for the 
        given specie main isoforms list"""

        path = cwd + '/gtf/' + s + '.gtf.gz'
        gtf = pd.read_table(path, compression='gzip', comment='#', low_memory=False, header=None)
        gtf = gtf.join(gtf[8].str.split('"', expand=True).add_prefix('sec'))           
        gtf = gtf[gtf['sec5'].isin(mainisoforms.mlist(s))]

        return gtf

        #mainisoforms.gtf('Homo_sapiens')
        
    def gtflist(s):
        """return the filtered gtf accessions list"""

        return mainisoforms.gtf(s).values.tolist()
    
        #mainisoforms.gtflist('Homo_sapiens')
        
    def IDS(s):
        """filter the gtf lines to retrieve only CDS and their 
        chromosome, strand, gene ID, transcript ID, protein ID,
        gene name, start value, stop value and return a list """

        gtflist = mainisoforms.gtflist(s)
        
        principals = []
        for l in gtflist:
            if 'CDS' in l[2]:
                chrom = l[0]
                strand = l[6]
                gid = l[10]
                tid = l[14]
                gname = l[20]
                pattern = gid.split('0')[0][:-1] + 'P'
                pid = [l2 for l2 in l if str(l2).startswith(pattern)]
            elif 'start' in l[2]:
                start = l[3]
            elif 'stop' in l[2]:
                stop = l[4]
            principals.append([chrom, gid, strand, pid[0], start, stop, tid, gname])

        return principals
    
        #mainisoforms.IDS('Homo_sapiens')
        
    def msorted(s):
        """sorting based on appearence 
        order along the genome"""

        principals = mainisoforms.IDS(s)
        
        tsv = pd.DataFrame(principals, columns=[
            'Chromosome', 
            'Gene_id', 
            'Strand', 
            'Protein_id', 
            'Start', 
            'Stop', 
            'Transcript_id', 
            'Gene_name'])

        tsv['Chromosome'] = pd.Categorical(
            tsv['Chromosome'], 
            ordered=True, 
            categories=ns.natsorted(
                tsv['Chromosome'].unique()))

        tsv['Mean'] = tsv[['Start', 'Stop']].mean(axis=1)
        tsv = tsv.sort_values(['Chromosome', 'Mean'], ascending=(True, True))
        tsv = tsv.drop(columns='Mean')
        
        return tsv

        #mainisoforms.msorted('Homo_sapiens')
        
    def totsv(s):
        """store sorted main isoforms 
        in a tsv format file"""

        tsv = mainisoforms.msorted(s)
        path = cwd + '/main/' + s + '.tsv'
        tsv.to_csv((path), index=False, sep='\t')
        
        #mainisoforms.totsv(s)
        
    def fafile(s):
        """FASTA file for a given specie"""

        path = cwd + '/fa/' + s + '.fa.gz'
        file = SeqIO.parse(gzip.open(path, 'rt'), 'fasta')
        return file

        #mainisoforms.fafile('Homo_sapiens')
        
    def transcriptfa(t, s):
        
        f = mainisoforms.fafile(s)
        for l in f:
            splitted = re.split('\s', l.description)
            transcript = splitted[4].split('transcript:')[1].split('.')[0]
            if transcript == t:
                return '>' + l.id + '\n' + l.seq
                
        #mainisoforms.transcriptfa('ENST00000003084', 'Homo_sapiens')
        
    def tofa(s):
        """write a FASTA file containing all 
        the principal isoforms for a given specie"""
        
        path = cwd + '/main/' + s + '_main.fa'
        mlist = mainisoforms.mlist(s)
        for transcript in mlist:
            t = mainisoforms.transcriptfa(transcript, s)
            print(t, file=open(path, 'a'))
            
        #mainisoforms.tofa('Homo_sapiens')
        
    def blast(s, t, e, h, m):
        """perform a BLASTP intraspecie for a given specie and given parameters"""

        path_in = cwd + '/main/' + s + '_main.fa'
        path_out = cwd + '/blast_queries/' + s + '_blast.txt'
        makeblastdb(dbtype='prot', input_file=(path_in))()

        blastp(query=path_in, db=path_in, num_threads=t, evalue=e, 
               max_hsps=h, outfmt=7, max_target_seqs=m, out=path_out)()

        #mainisoforms.blast('Homo_sapiens', 20, '10e-6', 1, 5) #num_threads, evalue, max_hsps, max_target_seqs