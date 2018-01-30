# -*- coding: utf-8 -*-
"""
Módulo responsável pela obtenção dos ficheiros dos genes de interesse e protéinas correspondentes
do NCBI e extração de algumas features e qualifiers úteis
"""
from Bio import Entrez, SeqIO

#cria os
def ficheirosGenes():
    Entrez.email = "sofiamgcarvalho@gmail.com"

    #posições dos genes de interesse no genoma
    lista_posicoes = [(3356827,3357996), (2964464, 2967667), (2263582, 2264706), (1283991, 1285145),
                     (2963973, 2964455), (1282693, 1283979), (2603700, 2604722), (534313, 535428), 
                     (531865, 533082), (532985, 534310), (531200, 531868), (530551, 531216)]
    
    lista_ficheiros = []
    
    #ficheiros com cada um dos genes de interesse
    for i in range(len(lista_posicoes)):
        #id do genoma no NCBI é 52840256
        handleGB = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id="52840256",
                           seq_start=lista_posicoes[i][0], seq_stop=lista_posicoes[i][1])
        seq_record = SeqIO.read(handleGB, "genbank")
        nome_ficheiro = 'sequenceGenbank'+ str(i) +'.gb'
        SeqIO.write(seq_record, nome_ficheiro, "genbank") #Guarda em formato genbank
        lista_ficheiros.append(nome_ficheiro)
        handleGB.close()
    
        #handleFASTA = Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text", id="52840256",
        #                            seq_start=listaPosicoes[i][0], seq_stop=listaPosicoes[i][1])
        #seq_record = SeqIO.read(handleFASTA, "fasta")
        #SeqIO.write(seq_record, 'sequenceFasta'+ str(i) +'.fasta', "fasta") #Guarda em formato fasta
        #handleFASTA.close()
    
        SeqIO.convert('sequenceGenbank'+ str(i) +'.gb', "genbank", 'sequenceF'+ str(i) +'.fasta', "fasta")

        #genes em ficheiros txt
        #txtFile= open('sequenceTxt.txt', 'w')
        #txtFile.write(str(seq_record.seq)) 
        #txtFile.close()
        
    return lista_ficheiros
        
def geneIdentification(lista):
    
    #lista = ["sequenceGenbank0.gb","sequenceGenbank1.gb","sequenceGenbank2.gb",
    #         "sequenceGenbank3.gb","sequenceGenbank4.gb","sequenceGenbank5.gb",
    #         "sequenceGenbank6.gb","sequenceGenbank7.gb","sequenceGenbank8.gb",
    #         "sequenceGenbank9.gb","sequenceGenbank10.gb","sequenceGenbank11.gb",]
    
    dic = {}
    
    for i in range(len(lista)):
        record = SeqIO.read(lista[i], "genbank")
        source_pos = str(record.features[0].location).replace("(+)","").replace("(-)","")
        
        for f in record.features:
            if f.type == "CDS":
                if source_pos == str(f.location).replace("(+)","").replace("(-)",""):
                    try:
                        gene_id = int(f.qualifiers["db_xref"][0].replace("GeneID:",""))
                    except:
                        gene_id = '-'
                    try:
                        locus_tag = f.qualifiers["locus_tag"][0]
                    except:
                        locus_tag = '-'
                    try:
                        name = f.qualifiers["gene"][0]
                    except:
                        name = '-'
                    try:
                        protein_id = f.qualifiers["protein_id"][0]
                    except:
                        protein_id = '-'
                    try:
                        ec_number = f.qualifiers["EC_number"][0]
                    except:
                        ec_number = '-'    
                    try:
                        product = f.qualifiers["product"][0] 
                    except:
                        product = '-'
                    try:
                        function = f.qualifiers["function"][0]
                    except:
                        function = '-'
                    try:
                        note = f.qualifiers["note"][0]
                    except:
                        note = '-'
                    strand = f.strand
                    dic[protein_id] = [gene_id, locus_tag, name, strand, ec_number, product, function, note]
    return dic

def ficheirosProteinas(dicionario):
    Entrez.email = "sofiamgcarvalho@gmail.com"
    
    i=0;
    lista_ficheiros = []
    
    for key in dicionario:
        handleGB = Entrez.efetch(db="protein", rettype="gb", retmode="text", id=key)
        seq_record = SeqIO.read(handleGB, "genbank")
        nome_ficheiro = 'sequenceProtGenbank'+ str(i) +'.gb'
        SeqIO.write(seq_record, nome_ficheiro, "genbank") #Guarda em formato genbank
        lista_ficheiros.append(nome_ficheiro)
        handleGB.close()
    
        SeqIO.convert('sequenceProtGenbank'+ str(i) +'.gb', "genbank", 'sequenceProtF'+ str(i) +'.fasta', "fasta")
        
        i+=1
        
    return lista_ficheiros
    

gene_files = ficheirosGenes()
print(gene_files)
print()
dic = geneIdentification(gene_files)
print(dic)

"""
import csv

with open('geneFeatures.csv', 'wb') as f:
    w = csv.DictWriter(f, str(dic.keys()))
    w.writeheader()
    w.writerow(dic)
"""

print()
protein_files = ficheirosProteinas(dic)
print(protein_files)
print()
                
                
            
    
    