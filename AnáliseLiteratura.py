# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 10:15:17 2018

@author: sofia
"""

from Bio import Entrez
Entrez.email = "A.N.Other@example.com"
handle = Entrez.egquery(term="Legionella pneumophila Philadelphia-1 Pathway: L-arginine biosynthesis I")
record = Entrez.read(handle)


for row in record["eGQueryResult"]:
    if row["DbName"]=="pubmed":
        #conta o número de ficheiros associados ao genoma, no NCBI
        nf = row["Count"]
        
from Bio import Entrez
handle = Entrez.esearch(db="pubmed", term="Legionella pneumophila Philadelphia-1 Pathway: L-arginine biosynthesis I",
                        retmax=nf)
record = Entrez.read(handle)
handle.close()
idlist = record["IdList"]
#imprime os ids dos ficheiros associados ao genoma
#print(idlist)

from Bio import Medline
handle = Entrez.efetch(db="pubmed", id=idlist, rettype="medline", retmode="text")
#records é um iterador
records = Medline.parse(handle)

#guardar os records em forma de lista
records = list(records)

#guardar as informações dos artigos num ficheiro
file = open("analiseLiteratura.txt","w")
for record in records:
    file.write("Title: " + str(record.get("TI", "?")) + "\n")
    file.write("Authors: " + str(record.get("AU", "?")) + "\n")
    file.write("Source: " + str(record.get("SO", "?")) + "\n" + "\n")
file.close()

