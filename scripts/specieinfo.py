from Bio import Entrez
import pandas as pd
import requests
import os
import sys

"""retrieve the current assembly 
set from ensembl server"""
link = "http://rest.ensembl.org/info/species?"
assemblies = requests.get(link, 
                            headers={
                                "Content-Type": 
                                "application/json"
                            }).json()

class specieinfo:
    
    def __init__(self, name):

        self.name     = name
        self.genus    = specieinfo.genus_func(self)
        self.specie   = specieinfo.specie_func(self)
        self.taxid    = specieinfo.taxid_func(self)
        self.pubs     = specieinfo.pubs_func(self)
        self.taxlist  = specieinfo.taxlist_func(self)
        self.taxorder = specieinfo.taxorder_func(self)
        self.taxclass = specieinfo.taxclass_func(self)
        self.assembly = specieinfo.assembly_func(self)
        self.allinfo  = specieinfo.allinfo_func(self)
    
    def assembly_func(self):
        """return assembly"""

        return [
            l['assembly'] 
            for l in assemblies['species'] 
            if l['name'] == self.name.lower()
        ][0]
        
    def genus_func(self):
        """return genus"""

        return self.name.split('_')[0]
    
    def specie_func(self):
        """return specie"""

        return '_'.join(
            self.name.split('_')[1:])
    
    def taxid_func(self):
        """return taxid"""

        taxid = Entrez.read(
            Entrez.esearch(
                db="taxonomy", 
                term=self.name)
        )['IdList']
        return taxid[0] if not list(
            taxid) == [] else None

    def pubs_func(self):
        """return publications number"""

        return  Entrez.read(
            Entrez.esearch(
                db="pubmed", 
                term=self.name)
        )['Count']

    def taxlist_func(self):
        """return taxonomy list from Entrez"""

        if self.taxid:
            df = pd.DataFrame(
                Entrez.read(
                    Entrez.efetch(
                        db="taxonomy", 
                        id=self.taxid)
                )[0]['LineageEx'])
            return df
        else: 
            return None
    
    def taxorder_func(self):
        """return order"""

        if self.taxlist is not None:
            taxorder = self.taxlist[
                self.taxlist['Rank'] == 'order'][
                'ScientificName'].values
            return taxorder[0] if not list(
                taxorder) == [] else None
        else:
            return None

    def taxclass_func(self):
        """return class"""

        if self.taxlist is not None:
            if 'Sauropsida' in self.taxlist[
            'ScientificName'].values:
                taxclass = 'Sauropsida'
            else:
                taxclass = self.taxlist[
                    self.taxlist['Rank'] == 'class'][
                    'ScientificName'].values
                taxclass = taxclass[0] if not list(
                    taxclass) == [] else None
            return taxclass
        else:
            return None
            
    def allinfo_func(self):
        """return a list containing: class, order, 
        genus, specie, publications, taxid, assembly"""
        
        return [
            self.taxclass,
            self.taxorder,
            self.genus,
            self.specie,
            self.pubs,
            self.taxid,
            self.assembly
        ]