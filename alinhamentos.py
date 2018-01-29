# -*- coding: utf-8 -*-
"""
Módulo responsável pelos alinhamentos e árvore filogenética
"""

from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio import Phylo

#lê o alinhamento efetuado no site Clustal Omega obtido para um ficheiro fasta com uma 
#proteína de interesse e as suas homólogas obtidas no BLAST  
def readAlinhamento(clustal_file):
    alClustal = AlignIO.read(clustal_file, "clustal")
    comprimento = alClustal.get_alignment_length()
    
    return alClustal,comprimento

#obtém a sequência consensus do alinhamento
def consensus(alClustal):
    summary_align = AlignInfo.SummaryInfo(alClustal)
    consensus_al = summary_align.dumb_consensus()
    
    return consensus_al

#lista das posições (colunas) onde o alinhamento é conservado e a sua percentagem
def conserved(alClustal):
    numcols = alClustal.get_alignment_length()
    numseqs = len(list(alClustal))
    conserved = []
    for c in range(numcols):
        col = alClustal[:,c]
        fc = col[0]
        if col.count(fc)==numseqs:
            conserved.append(c)      
    
    percentagem = len(conserved) / numcols * 100
    
    return conserved,percentagem

#determina a zona do alinhamento mais conservada, bem como as suas posições de início e de fim
def zonaMaisConservada(alClustal, conserved):
    inicio = 0
    end = 0
    best_inicio = 0
    best_end = 0
    maxs = 0
    pos = 1
    while pos < len(conserved):
        if conserved[pos] - conserved[pos-1] > 1:
            inicio = conserved[pos]
        end = conserved[pos]
        if end - inicio > maxs: 
            maxs = end - inicio
            best_inicio = inicio
            best_end = end
        pos = pos + 1

    best_conserved = alClustal[:,best_inicio:best_end+1]
    
    return best_conserved,best_inicio,best_end

#imprime a árvore filogenética obtida para o alinhamento
def printTree(tree_file):
    tree = Phylo.read(tree_file, "newick")
    Phylo.draw_ascii(tree)
    
    return tree

#execução das várias funcionalidades de alinhamentos implementadas
op = input("Deseja analisar o alinhamento (1) ou a árvore filogenética (2)?")
if op == "1":
    clustal_file = input("Indique o nome do ficheiro com o alinhamento:")
    #clustal_file = "alinhamentoHomologasProt0.clustal"
    
    alClustal, comprimento = readAlinhamento(clustal_file)
    consenso = consensus(alClustal)
    conservadas, percentagem = conserved(alClustal)
    mais_conservada, inicio, fim = zonaMaisConservada(alClustal, conservadas)
    consenso_conservada = consenso[inicio:fim+1]
    
    print("Alinhamento:", alClustal)
    print("Comprimento:", comprimento)
    print("Consensus:", consenso)

    print("\nPosições onde o alinhamento é conservado: ", conservadas)
    print("Percentagem de posições onde o alinhamento é conservado: ", percentagem)

    print("\nZona do alinhamento mais conservada:", mais_conservada)
    print("Posições de início e de fim: " , inicio, " ", fim)
    print ("Consensus da zona mais conservada:\n"  + consenso_conservada)
    
elif op == "2":
    tree_file = input("Indique o nome do ficheiro com a árvore:")    
    #tree_file = "arvoreHomologasProt0.txt"

    tree = printTree(tree_file)

