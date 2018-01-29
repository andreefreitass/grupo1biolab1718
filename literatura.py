# -*- coding: utf-8 -*-
"""
Módulo responsável pela análise de literatura
"""

from Bio import Entrez
from Bio import Medline

#permite criar um ficheiro com as informações dos artigos para uma dada pesquisa na PubMed
def criaFicheiroLiteratura(query, file_name):
    Entrez.email = "sofiamgcarvalho@gmail.com"
    handle = Entrez.egquery(term=query)
    record = Entrez.read(handle)

    for row in record["eGQueryResult"]:
        if row["DbName"]=="pubmed":
            #conta o número de ficheiros associados à query
            nf = row["Count"]
        
        handle = Entrez.esearch(db="pubmed", term=query, retmax=nf)
        record = Entrez.read(handle)
        handle.close()
        idlist = record["IdList"]
        #imprime os ids dos ficheiros associados à query
        #print(idlist)

    handle = Entrez.efetch(db="pubmed", id=idlist, rettype="medline", retmode="text")
    #records é um iterador
    records = Medline.parse(handle)

    #guardar os records em forma de lista
    records = list(records)

    #guardar as informações dos artigos num ficheiro
    file = open(file_name,"w")
    for record in records:
        file.write("Title: " + str(record.get("TI", "?")) + "\n")
        file.write("Authors: " + str(record.get("AU", "?")) + "\n")
        file.write("Source: " + str(record.get("SO", "?")) + "\n" + "\n")
    file.close()
    
    return file_name

#literatura para o organismo
organismo = criaFicheiroLiteratura("Legionella pneumophila subsp. pneumophila str. Philadelphia 1", "literaturaOrganismo.txt")
print("Ficheiro " + organismo + " guardado")

#literatura para a via metabólica
via = criaFicheiroLiteratura("Legionella pneumophila Philadelphia-1 Pathway: L-arginine biosynthesis I", "literaturaVia.txt")
print("Ficheiro " + via + " guardado")