#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 10:06:08 2017

@author: giovanni
"""


from Bio import Entrez
Entrez.email = 'giovanni.pasquini2@unibo.it'


# Read file having gene names as a column
with open('/media/giovanni/500GbVolume/ANNA/Pan_Cancer/gene_list.csv') as input_file:
    for geneName in input_file:
        
        # Prepare query
        geneName=geneName.rstrip()
        
        query="(" + geneName + "[gene]) AND (Homo sapiens[orgn])"
    
        # Query NCBI Gene
        handle = Entrez.esearch(db="gene" ,term=query, tettype="gene", retmode="xml")
        
        record = Entrez.read(handle)
        
        record["IdList"][0]
        
        handle = Entrez.efetch(db="gene" ,id=record["IdList"][0], tettype="gene", retmode="xml")
        
        record = Entrez.read(handle)
        
        
        
        # Get Gene Location
        geneLoc=record[0]['Entrezgene_location'][0]['Maps_display-str']
    
        print(geneName, geneLoc)









