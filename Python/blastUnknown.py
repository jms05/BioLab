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

'''


'''
def blast(seq,database="swissprot"):
	print("blasting....")
	#fasta_string = open("m_cold.fasta").read()
	result=NCBIWWW.qblast("blastp",database,seq)
	print("blast done")
	save_file=open("blastteste.xml","w")
	save_file.write(result.read())
	save_file.close()
	result.close()
	
	
	resultxml = open("blastteste.xml") 
	blast= NCBIXML.parse(resultxml)
	
	E_VALUE_THRESH = 0.5
	for blast_record in blast:
		for alignment in blast_record.alignments:
			for hsp in alignment.hsps:
				print(str(hsp)) 
				print("E_VALUE_THRESH: " + str(hsp.expect))
				if hsp.expect < E_VALUE_THRESH:
					print ('\n')
					print ('****Alignment****')
					print ('sequence:'+ str(alignment.title))
					print ('e value:'+  str(hsp.expect))
					print (hsp.query[0:75])
					print (hsp.match[0:75])
					print (hsp.sbjct[0:75])
	input()


def parseCSV(filename):
	filecsvlines = open(filename,"r").read().splitlines()
	sep = filecsvlines[0].split("=")[1]
	filecsvlines= filecsvlines[1:]
	dic = {}
	for line in filecsvlines:
		campos = line.split(sep)
		status  = campos[17].upper()
		idprot = campos[0]
		seqAA = campos[9]
		if "unknown".upper() in status: 
			dic[idprot] = seqAA

	return dic
lpeq =9999
peq =""
d = parseCSV("tabela_com_seq_Completa.csv")
for k in d.keys():
	print(k +"="+d[k])
	if(len(d[k])<lpeq):
		lpeq = len(d[k])
		peq=d[k]

print(str(len(d)))
print("PEQ: " + peq)
blast(peq)

