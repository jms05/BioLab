# -*- coding: utf-8 -*-

from Bio import Entrez
from Bio import Medline


def ProcurarArtigos():
    # Função que faz procura no PubMed de todos os artigos que contêm informações acerca da Neisseria gonorrhoeae
    Entrez.email = "joaosilva1993@live.com.pt"     
    handle = Entrez.egquery(term = "Neisseria gonorrhoeae")
    record = Entrez.read(handle)
    for row in record["eGQueryResult"]:
        if row["DbName"] == "pubmed":
            print "Numero de artigos no PubMed: ", (row["Count"])
            
    return record

def DownloadArtigos():
    handle = Entrez.esearch(db = "pubmed", term = "Neisseria gonorrhoeae", retmax = 11110) # 11110 existentes
    record = Entrez.read(handle)
    idlist = record["IdList"]
    print(idlist)
    
    return idlist

def ObterRegistos(idlist):
    # para obter registos Medline e extrair informação
    handle = Entrez.efetch(db = "pubmed", id = idlist, rettype = "medline", retmode = "text")
    records = Medline.parse(handle)
    
    # guarda os registo e estes são convertidos numa lista
    records = list(records)
    
    # Dá-nos algumas informações sobre os records
    for record in records:
        print ("Title: ", record.get("TI", "?"))
        print ("Authors: ", record.get("AU", "?"))
        print ("Source: ", record.get("SO", "?"))
        print ("")
        
    return records

def ProcuraTitulo(records):
    search_title = raw_input("Qual o artigo que procura: ") #Introduzir titulo que pretendemos procurar
        
    for record in records:
        if not "TI" in record: #TI-Title
            continue
        if search_title in record["TI"]:
            print ("%s encontrado: %s" % (search_title, record["SO"]))
            
def ProcuraAutor(records):
    search_author = raw_input("Qual o autor que procura: ") #Introduzir o nome do autor que procuramos
        
    for record in records:
        if not "AU" in record:  #AU-Authors
            continue
        if search_author in record["AU"]:
            print ("Autor %s encontrado: %s" % (search_author, record["SO"]))
            
            
def teste():
    lista = DownloadArtigos()
    registos = ObterRegistos(lista)
    ProcuraTitulo(registos)
    ProcuraAutor(registos)

        

if __name__ == '__main__':
    teste()
            
        
        