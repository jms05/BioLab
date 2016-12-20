from Bio import SeqIO
import requests
from io import StringIO 
import re

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

def parseFile(filename):
	seq_record = SeqIO.read(filename, "genbank")
	return seq_record


def findStrand(r):
	strand=""
	aux = str(re.findall("location: \[.+\:.+\]\\([+-]\)",r))
	i=2
	for i in xrange(1,10):
		pass
   	for i in range(len(aux)):
   		if aux[i]=='(':
   			strand=aux[i+1]
   	print(strand)
   	return -1

def blastuniprot(seq):
	#"http://www.uniprot.org/blast/uniprot/B201612208A530B6CA0138AFAA6D2B97CE8C2A924B8F0ECU"
	r = requests.get("http://www.uniprot.org/blast/uniprot/"+seq+"/")
	print(str(r))

def procuraUni(idp):
	params = {"query": idp, "format": "fasta"}
	r = requests.get("http://www.uniprot.org/uniprot/", params)
	for record in SeqIO.parse(StringIO(r.text), "fasta"):
		print(record)
	a=input("lixo")

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






