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
	functions = []
	function = ""
	for line in textxml.splitlines():
		if "<comment type=\"function\">" in line:
			starf  =True
			function=""
		if(("</comment>" in line) and starf):
			starf = False
			functions.append(function)

		if(starf):
			function+=(line.strip())

	ret = []
	for fun in functions:
		ret.append(fun.split("<")[-2].split(">")[1])	
	return(ret)




def parseTXT(record):
	try:
		status = str(record.data_class)
	except Exception:
		status= "-"
	local = "--"
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
	functions  = parseXML(xml.decode('utf-8'))
	return(status,local,funcMol,bioPro,functions)
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
				strand="reverse"
				seqdnaprot= reverseComplement(dna[start:end])
			else:
				strand="forward"
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
		break
	parst = parseSwissProt(idfasta)
	return(idfasta,parst)









def blast(seq):
	print("TO Blast :" + seq )
	#fasta_string = open("m_cold.fasta").read()
	result=NCBIWWW.qblast("blastp","swissprot ",seq)
	
	save_file=open("blastteste.xml","w")
	save_file.write(result.read())
	save_file.close()
	result.close()
	
	
	resultxml = open("blastteste.xml") 
	blast= NCBIXML.parse(resultxml)
	
	E_VALUE_THRESH = 0.05
	for blast_record in blast:
	    for alignment in blast_record.alignments:
	        for hsp in alignment.hsps:
	            if hsp.expect < E_VALUE_THRESH:
	                print ('\n')
	                print ('****Alignment****')
	                print ('sequence:'+ str(alignment.title))
	                print ('e value:'+  str(hsp.expect))
	                print (hsp.query[0:75])
	                print (hsp.match[0:75])
	                print (hsp.sbjct[0:75])
	raw_input()


rec = parseFile("grupo6.txt")
datas = getinfosfromgem(rec)

for data in datas:
	(gene,prot,ec)=data
	(rev,ID_prot,prorainNAme,tam,local,function,tradu)=prot
	print ("-"*20)
	print ("EC: " + str(ec))
	print ("GEN: " + str(gene))
	print ("Prot: " + str(prot))
	swissinfo = getInfouniprot(ID_prot)
	print("UNIPROT id :" + str(swissinfo[0]))
	print("UNIPROT uniprot :" + str(swissinfo[1:]))
	print ("-"*20)
	a  =input()
'''
seq=""

i=0
for feat in rec.features:
	findStrand(str(feat))
	if feat.type == 'CDS':
		i+=1
		print ("xref -> " + feat.qualifiers["db_xref"][0])
		print ("ID_prot -> " + feat.qualifiers["protein_id"][0])
		try:
			print ("Funcao-> " + feat.qualifiers["function"][0])
		except Exception:
			print ("Funcao-> Unknown")
		
		print ("Traducao -> " + feat.qualifiers["translation"][0])
		procuraUni(feat.qualifiers["protein_id"][0])
		blastuniprot(feat.qualifiers["translation"][0])
		print("-"*20)
		#if(feat.qualifiers["function"][0]!= "Unknown"):
			#blast(feat.qualifiers["translation"][0])
		

print ("Total: " + str(i))
'''





