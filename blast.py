# -*- coding: utf-8 -*-
"""
Módulo responsável pela execução do BLAST
"""

from Bio.Blast import NCBIXML 
from Bio.Blast import NCBIWWW 
from Bio import Entrez, SeqIO
import os
import re

#executar o blast para um ficheiro de uma proteína e guardar resultados num ficheiro xml 
def blastProt(database, file_name, file_format):
    record = SeqIO.read(open(file_name), format=file_format)
    
    print("BLAST runnning")
    result_handle = NCBIWWW.qblast("blastp", database, record.format(file_format))
    print("BLAST finnished")

    mo = re.search("\d", file_name)
    number = mo.group()
    
    xml_file = "blast-Prot" + number + ".xml"
    save_file = open(xml_file, "w")
    save_file.write(result_handle.read()) 
    save_file.close() 
    result_handle.close()
    
    return xml_file

#guarda informações do ficheiro xml do BLAST dos n melhores alinhamentos e retorna os ids das proteínas desses alinhamentos
def resultadosBlast(xml_file, n):
    mo = re.search("\d", xml_file)
    number = mo.group()
    
    file = open("blastBest" + str(n) + "-Prot" + number + ".txt", "w")
    
    ids_proteinas = []
    
    result_handle = open(xml_file)
    blast_records = NCBIXML.parse(result_handle)
    first_record = next(blast_records)
    
    file.write("Number of hits: " + str(len(first_record.alignments)) + "\n")
    file.write("Database: " + first_record.database + "\n")
    file.write("Matrix: " + first_record.matrix + "\n\n")
    
    i=0
    
    for alignment in first_record.alignments:
        ids_proteinas.append(alignment.accession)
        file.write("Accession: " + alignment.accession + "\n")
        file.write("Number of HSPs: " + str(len(alignment.hsps)) + "\n")
        file.write("Alignment length: " + str(alignment.length) + "\n")
        file.write("Definition: " + alignment.hit_def + "\n")
        for hsp in alignment.hsps:
            file.write("E-value: " + str(hsp.expect) + "\n")
            file.write("Score: " + str(hsp.score) + "\n")
            file.write("Identities: " + str(hsp.identities) + "\n")
            file.write("Query: " + hsp.query + "\n")           
            file.write("Match: " + hsp.match + "\n")
            file.write ("Subject: " + hsp.sbjct + "\n\n")
        i+=1
        if i == n:
            break
        
    return ids_proteinas

#cria ficheiro fasta com a proteína original e as proteínas homólogas do BLAST
def proteinasBlastToFasta(query_file, ids_proteinas):
    proteinas_files = []
    proteinas_files.append(query_file)
    
    Entrez.email = "a74619@alunos.uminho.pt"
    
    for i in range(len(ids_proteinas)):
        handle = Entrez.efetch(db="protein", rettype="fasta", retmode="text", id=ids_proteinas[i])
        seq_record = SeqIO.read(handle, "fasta")
        nome_ficheiro = str(i) + ".fasta"
        SeqIO.write(seq_record, nome_ficheiro, "fasta")
        proteinas_files.append(nome_ficheiro)
        handle.close()
    
    mo = re.search("\d", query_file)
    number = mo.group()
    
    file_homologas = 'homologasProt' + number + '.fasta' 
    with open(file_homologas,'w') as outfile:
        for fname in proteinas_files:
            with open(fname) as infile:
                outfile.write(infile.read())
                
    for i in range(1,len(proteinas_files)):
        try:
            os.remove(proteinas_files[i])
        except OSError:
            pass
        
    return file_homologas

#executar as funções anteriores para todas as proteínas
def blastAllProts(database, files, file_format, n):
    files_homologas = []
    
    for i in range(len(files)):
        xml_file = blastProt(database, files[i], file_format)
        ids = resultadosBlast(xml_file, n)
        file_h = proteinasBlastToFasta(files[i], ids)
        files_homologas.append(file_h)
    
    return files_homologas

#ficheiros correspondentes às sequências das proteínas dos nossos genes
files = ["sequenceProtF0.fasta","sequenceProtF1.fasta","sequenceProtF2.fasta",
         "sequenceProtF3.fasta","sequenceProtF4.fasta","sequenceProtF5.fasta",
         "sequenceProtF6.fasta","sequenceProtF7.fasta","sequenceProtF8.fasta",
         "sequenceProtF9.fasta","sequenceProtF10.fasta","sequenceProtF11.fasta",]

#execução das várias funcionalidades do BLAST implementadas             
op = input("Deseja executar BLAST para todas as sequências? (S ou N)")
if op == "S":
    db = input("Base de dados que deseja: (nr ou swissprot por exemplo)")
    n = input("Melhores n alinhamentos que deseja:")
    lista = blastAllProts(db, files, "fasta", n)
    print(lista)
elif op == "N":
    db = input("Base de dados que deseja: (nr ou swissprot por exemplo)")
    n = input("Melhores n alinhamentos que deseja:")
    f_name = input("Nome do ficheiro com a sequência da proteína:")
    xml_file = blastProt(db, f_name, "fasta")
    ids = resultadosBlast(xml_file,n)
    file_homologas = proteinasBlastToFasta(f_name, lista)
    print(file_homologas)
    
    
    
    

        
