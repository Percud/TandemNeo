from specieinfo import specieinfo as si
from ftplib import FTP
import os

cwd = os.path.dirname(os.getcwd())

class download:

    def downl(s, k):
        """download FASTA and GTF from ensembl ftp 
        servers for a given specie, 
        must be specified if FASTA or GTF"""

        ftp = FTP('ftp.ensembl.org')

        if k == 'fa':
            path = '/pub/current_fasta/' + s.lower() + '/pep/'
            retr = 'RETR ' + s + '.' + si.assembly(s) + '.pep.all.fa.gz'
            wr = open(cwd + '/fa/' + s + '.fa.gz', 'wb').write

        elif k == 'gtf':
            path = '/pub/current_gtf/' + s.lower() + '/'
            rs = si.assembly(s) + '.' + str(si.release())
            retr = "RETR " + s + '.' + rs + '.gtf.gz'
            wr = open(cwd + '/gtf/' + s + '.gtf.gz', 'wb').write

            ftp.login()
            ftp.cwd(path)
            ftp.retrbinary(retr, wr)

        ftp.quit()
        
    #download.downl('Homo_sapiens', 'fa')