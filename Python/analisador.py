from Bio import SeqIO
import requests
from io import StringIO 
import re

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML


from Bio import ExPASy
import urllib
from Bio import SwissProt
from Bio.SwissProt import KeyWList
from Bio import ExPASy

import os
import random
from urllib.request import urlopen

DNA_BASES = ["A","T","C","G"]
PROTAIN_BASES= ['I', 'M', 'T', 'N', 'K', 'S', 'R',
				'P', 'H', 'Q', 'V', 'A', 'D', 'E',
				'G', 'F', 'L', 'Y', 'C', '_', 'W',]
START_BASE = 'M'
END_BASE = '_'

DNA_COMPLEMENT = {"A":"T","T":"A","C":"G","G":"C"}
CODONS = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }


def parseFile(filename):
	seq_record = SeqIO.read(filename, "genbank")
	return seq_record


def reverseComplement(dnaSeq):
	comp=""
	for base in dnaSeq:
		comp+=DNA_COMPLEMENT[base]
	return comp[::-1]


def parseXML(textxml):
	starf= False
	satrtG = False
	startFeat = False
	functions = []
	function = ""
	name  =""
	domains = []
	for line in textxml.splitlines():
		if "<comment type=\"function\">" in line:
			starf  =True
			function=""
		if(("</comment>" in line) and starf):
			starf = False
			functions.append(function)

		if(starf):
			function+=(line.strip())

		if "<feature" in line:
			startFeat=True
		if "</feature>" in line:
			startFeat=False
		if(startFeat):
			if "type=\"domain\"" in line:
				campos = line.split()
				for campo in campos:
					if "description" in campo:
						interesse = campo.split("=")[1].replace("\"","")
						domains.append(interesse)
		if "<gene>" in line:
			satrtG=True
		if "</gene>" in line:
			satrtG=False
		if(satrtG):
			if "primary" in line:
				name = line.split(">")[1].split("<")[0]

	ret = []
	for fun in functions:
		ret.append(fun.split("<")[-2].split(">")[1])	
	if(name ==""):
		name = "-"

	domianTex = ""
	for domain in domains:
		domianTex+=(domain+"|")

	if domianTex=="":
		domianTex="-"
	return(ret,name,domianTex)




def parseTXT(record):
	try:
		status = str(record.data_class)
	except Exception:
		status= "-"
	local = "-"
	name = "-"
	funcMolec = []
	bioPro = []
	for cr in record.cross_references:
		if(cr[0]== "GO"):
			(tipo,ids,cool,pis) =cr
			if(tipo=="GO"):
				cools = str(cool).split(":")
				if(cools[0]=='F'):
					funcMolec.append(cools[1])
				if (cools[0]=='P'):
					bioPro.append(cools[1])
				if (cools[0]=='C'):
					local=cools[1]
	return (status,local,funcMolec,bioPro)

def downloadSwiss(idfasta,ext):
	target_url = "http://www.uniprot.org/uniprot/" + idfasta+ "."+ext
	ret = urlopen(target_url).read()
	return ret

def parseSwissProt(idswiss):
	txt  =downloadSwiss(idswiss,"txt")
	xml = downloadSwiss(idswiss,"xml")
	#ficheirs ja aqui
	f =open("tmp.dat","w")
	f.write(txt.decode('utf-8'))
	f.close()
	handle = open("tmp.dat")
	(status,local,funcMol,bioPro) = parseTXT(SwissProt.read(handle))
	#parse ao txt feito
	(functions,name,domianTex)  = parseXML(xml.decode('utf-8'))
	return(status,local,funcMol,bioPro,functions,name,domianTex)
	#print(str(record))


def getinfosfromgem(genbank):
	ACC = "NC_002942"
	ret =[]
	i=0
	dna = genbank.seq
	for feat in rec.features:
		#pri(str(feat))
		if feat.type == 'CDS':
			strand  = str(feat.location).split("]")[1]
			inter= str(feat.location).split("]")[0].replace("[","")
			start = int(inter.split(":")[0])
			end = int(inter.split(":")[1])
			
			if(strand=="(-)"):
				strand="(-)"
				seqdnaprot= reverseComplement(dna[start:end])
			else:
				strand="(+)"
				seqdnaprot= str((dna[start:end]))

			geneID = feat.qualifiers["db_xref"][0]
			ID_prot = feat.qualifiers["protein_id"][0]
			try:
				function= feat.qualifiers["function"][0]
			except Exception:
				function= "Unknown"
			try:
				genName= feat.qualifiers["gene"][0]
			except Exception:
				genName= "-"

			tradu = feat.qualifiers["translation"][0]
			locus = feat.qualifiers["locus_tag"][0]
			prorainNAme = feat.qualifiers["product"][0]
			try:
				ecNumber = feat.qualifiers["EC_number"][0]
			except Exception:
				ecNumber= "-"
			geninfo = (geneID,ACC,locus,genName,strand,seqdnaprot)
			protinfo = ("uniorotID_rev",ID_prot,prorainNAme,len(tradu),"local",function,tradu)
			ret.append((geninfo,protinfo,ecNumber))
	return ret



def getInfouniprot(protainID):
	params = {"query": protainID, "format": "fasta"}
	r = requests.get("http://www.uniprot.org/uniprot/", params)
	i=0
	for record in SeqIO.parse(StringIO(r.text), "fasta"):
		idfasta = record.id
		idfasta = str(idfasta).split("|")[1]
		#break
	parst = parseSwissProt(idfasta)
	return(idfasta,parst)



def shortSeq(sequence,tam):
	return sequence # tirar isto para fazer short
	start = sequence[:tam]
	end = sequence[(-1*tam):]
	return(start+"..."+end)

def juntaFuncoes(listaF,fungb,sep):
	ret = fungb
	for f in listaF:
		ret += (sep+f)
	return ret

def juntaLista(lista,sep):
	ret = ""
	for elem in lista:
		ret += elem+sep
	return ret[:-1]

def createCVSRecord(gbData,swissData,sep):
	grauRev = "---"
	(genInfo,protInfo,EC) = gbData
	(geneID,NCIgual,locusTag,geneName,strand,dnaSeq) =genInfo
	(lixo,assNCBI,protName,protLen,lixo2,protFungb,protSeq) = protInfo
	(idSwiss,parse) = swissData
	(protStatus,protLocal,funcaoMolec,processBiol,funcoes,geneNameSwiss,domianTex) = parse
	funcaoMolec = juntaLista (funcaoMolec, "_")
	processBiol = juntaLista( processBiol,"_")
	if(funcaoMolec==""):
		funcaoMolec="-"
	if(processBiol==""):
		processBiol="-"
	funcoes = juntaFuncoes(funcoes,protFungb,"_")
	geneID=geneID.split(":")[1]
	##feito na
	if(geneName=="-"):
		geneName = geneNameSwiss
	data = geneID+sep+geneName+sep+NCIgual+sep+locusTag+sep+strand+sep+shortSeq(dnaSeq,7)
	data = data + sep + assNCBI+ sep+idSwiss+sep+protName+sep+shortSeq(protSeq,7)
	data = data+ sep+ str(protLen)+ sep+ protStatus+ sep + grauRev + sep + protLocal
	data = data +sep +EC+ sep + funcaoMolec+ sep+ processBiol+ sep + funcoes+ sep + domianTex

	return data


sep =";"
rec = parseFile("grupo6.txt")
datas = getinfosfromgem(rec)
filecsv = open("tabela_com_seq_Completa.csv","w")
filecsv.write("sep="+sep+"\n")

cabeca="geneID"+sep+"GeneName"+sep+"GeneAccessNumber"+sep+"locusTag"+sep+"strand"+sep+"DNA_SEQ"
cabeca+=sep+"AccessNumberNCBI"+sep+"idSwiss"+sep+"protName"+sep+"PROT_SEQ"
cabeca+=sep+"PROT_Tamanho"+sep+"protStatus"+sep+"grauRev"+sep+"protLocal"
cabeca+=sep+"EC"+sep+"(GO)funcaoMolec"+sep+"(GO)processBiol"+sep+"funcoes"+sep+"Domain"

filecsv.write(cabeca+"\n")

i=1
for data in datas:
	(gene,prot,ec)=data
	(rev,ID_prot,prorainNAme,tam,local,function,tradu)=prot
	swissinfo = getInfouniprot(ID_prot)
	dataCS  = createCVSRecord(data,swissinfo,sep)
	filecsv.write(dataCS+"\n");
	print("mais UMA: "  + str(i)) 
	i=i+1

filecsv.close()