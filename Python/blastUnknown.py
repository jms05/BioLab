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


def juntaLista(lista,sep):
	ret = ""
	for elem in lista:
		ret += elem+sep
	return ret[:-1]


def shortSeq(sequence,tam):
	return sequence # tirar isto para fazer short
	start = sequence[:tam]
	end = sequence[(-1*tam):]
	return(start+"..."+end)
#<fullName>Prolyl endopeptidase</fullName>
def parseXML(textxml):
	starf= False
	satrtG = False
	startFeat = False
	functions = []
	function = ""
	geneNAme = ""
	organismSientif = ""
	organisComun = ""
	namesProt =[]
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
		if("</fullName>" in line):
				nameProt= line.split(">")[1].split("<")[0]
				namesProt.append(nameProt)
		if("<name type=\"scientific\">" in line):
			#<name type="scientific">Sus scrofa</name>
			organismSientif= line.split(">")[1].split("<")[0]

		if("<name type=\"common\">" in line):
			#<name type="scientific">Sus scrofa</name>
			organisComun= line.split(">")[1].split("<")[0]

		if "<gene>" in line:
			satrtG=True
		if "</gene>" in line:
			satrtG=False
		if(satrtG):
			if "primary" in line:
				geneNAme = line.split(">")[1].split("<")[0]




	ret = []
	for fun in functions:
		ret.append(fun.split("<")[-2].split(">")[1])	
	if(geneNAme ==""):
		geneNAme = "-"

	domianTex = ""
	for domain in domains:
		domianTex+=(domain+"|")

	if domianTex=="":
		domianTex="-"

	print(str(namesProt))
	if(len(namesProt))==0:
		protName="-"
	else:
		protName=juntaLista(namesProt,"|")
	print("Protaina: " +protName)
	hostName = organismSientif+"("+organisComun+")"
	return(ret,geneNAme,domianTex,protName,hostName)


def parseTXT(record):
	try:
		status = str(record.data_class)
	except Exception:
		status= "-"
	local = "-"
	name = "-"
	try:
		sequence = shortSeq(str(record.sequence),5)
	except Exception:
		sequence= "-"
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
	funcoes =""
	for fun in funcMolec:
		funcoes+= fun+"|"
	if(funcoes==""):
		funcoes="-"
	else:
		funcoes=funcoes[:-1]
	return (status,local,funcoes,bioPro,sequence)

def downloadSwiss(idfasta,ext):
	target_url = "http://www.uniprot.org/uniprot/" + idfasta+ "."+ext
	ret = urlopen(target_url).read()
	return ret

def parseSwissProt(idswiss):
	txt  =downloadSwiss(idswiss,"txt")
	xml = downloadSwiss(idswiss,"xml")
	#ficheirs ja aqui
	f =open("tmp_Blast.dat","w")
	f.write(txt.decode('utf-8'))
	f.close()
	handle = open("tmp_Blast.dat")
	(status,local,funcMol,bioPro,sequence) = parseTXT(SwissProt.read(handle))
	#parse ao txt feito
	(functions,name,domianTex,protName,hostName)  = parseXML(xml.decode('utf-8'))
	return(status,local,funcMol,bioPro,functions,name,domianTex,sequence,protName,hostName)
	#print(str(record))


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



def createCVSBlast(idprot,idBlastMatch,evalu,score,dataMatch,sep):
	grauRev = "---"
	(protStatus,protLocal,funcaoMolec,processBiol,funcoes,geneNameSwiss,domianTex,sequence,protName,hostName) = dataMatch
	funcaoMolec = juntaLista (funcaoMolec, "|")
	processBiol = juntaLista( processBiol,"|")
	#protName = juntaLista(protName,"|")
	funcoes = juntaLista(funcoes,"|")

	if(funcaoMolec==""):
		funcaoMolec="-"
	if(processBiol==""):
		processBiol="-"
	if(funcoes==""):
		funcoes="-"
	if(protName==""):
		protName="-"
	line  = idprot + sep 
	line += idBlastMatch+ sep
	line += str(evalu)+ sep
	line += str(score)+ sep 
	line +=protName+sep
	line +=protStatus +sep 
	line +=grauRev + sep
	line +=hostName+sep 
	line += protLocal + sep 
	line += funcaoMolec + sep 
	line += processBiol + sep 
	line += funcoes + sep 
	line += geneNameSwiss+ sep 
	line += domianTex + sep +sequence
	print("LINHA: " + line)
	return line

def getInfosFromBlastResult(idsBlasts):
	ret ={}
	print ("results to get info: " + str(idsBlasts))
	for (idB,evalu,score) in idsBlasts:
		print ("Getting info for blast id: " + idB) 
		(idfasta,details) = getInfouniprot(idB)
		print("Got this details: "+ str(details) )
		ret[idfasta]=(details,evalu,score)

	return ret


def blast(seq,database="swissprot",maxMatch=5):
	matches=0
	idsGoodBlast= []
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
	SCORE =95
	for blast_record in blast:
		for alignment in blast_record.alignments:
			for hsp in alignment.hsps:
				#print(str(hsp)) 
				if hsp.expect < E_VALUE_THRESH and hsp.score> SCORE:
					print("NICE\n\tE_VALUE_THRESH: " + str(hsp.expect) + " SCORE: " + str(hsp.score))
					data = (alignment.title.split("|")[3],hsp.expect,hsp.score)
					idsGoodBlast.append(data)
					matches +=1
					if (matches>=maxMatch):
						return idsGoodBlast
					#print ('\n')
					#print ('****Alignment****')
					#print ('sequence:'+ str(alignment.title))
					#print ('e value:'+  str(hsp.expect))
					#print (hsp.query[0:75])
					#print (hsp.match[0:75])
					#print (hsp.sbjct[0:75])
	
	return idsGoodBlast




def doRandomDelay(minTime,maxTime):
	#TODO
	return

def tryGetInfoFrom(idprot,ProtSeq,sep):
	csvFileLines  =[]
	blastResults = blast(ProtSeq) #list of id,eval,score
	dicInfosFromBlast  =getInfosFromBlastResult(blastResults)
	i=0
	for idProtainMatch in dicInfosFromBlast.keys():
		i+=1
		(detailsMatch,evalu,score) = dicInfosFromBlast[idProtainMatch]
		stringCSV = createCVSBlast(idprot,idProtainMatch,evalu,score,detailsMatch,sep)
		csvFileLines.append(stringCSV)
	if(i==0):
		csvFileLines.append(idprot+sep+"NoResult")
	print("CSVLINES: " +str(csvFileLines) )
	return csvFileLines



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

	return (dic,sep)


def readinfileBlast(filename):
	try:
		filecsvlines = open(filename,"r").read().splitlines()
		sep = filecsvlines[0].split("=")[1]
		filecsvlines= filecsvlines[1:]
		idsProcess=[]
		for line in filecsvlines:
			campos = line.split(sep)
			idprot = campos[0]
			if idprot not in idsProcess and "198" in  idprot:
				idsProcess.append(idprot)
		return idsProcess

	except Exception:
		return None

	
def startCSVFile(filename,sep,head):
	newCSVfile = open (filename,"w")
	newCSVfile.write("sep=" + sep+"\n")
	newCSVfile.write(head+"\n")
	newCSVfile.close()

def saveResultsonCSV(filename,lines):
	newCSVfile = open (filename,"a")
	for line in lines:
		print("Savig line " + line)
		newCSVfile.write(line+"\n")
		print("line Saved")
	newCSVfile.close()

	


minTime = 30
maxTime = 60
blastCSVName = "BlastCSV.csv"
(d,sep) = parseCSV("tabela_com_seq_Completa.csv")

todoBlast =  len(d)

newCSVHead = "Idprot" + sep + "idBlastMatch"+  sep+"evalu"+ sep+ "Score"+ sep +"protName"+sep+ "protStatus" +sep + "grauRev" + sep +"host"+sep
newCSVHead+= "protLocal" + sep + "funcaoMolec" + sep +"processBiol" + sep +"funcoes "+ sep
newCSVHead+= "geneNameSwiss"+ sep +"domianTex" + sep + "sequence"


blastDone = readinfileBlast(blastCSVName)
if(blastDone==None):
	startCSVFile(blastCSVName,sep,newCSVHead)
	blastDone=[]

todoBlast =  len(d) -len(blastDone)
print(str(todoBlast)+ " BLASTS TO DO")


for idProt in d.keys():
	if(idProt not in blastDone):
		print(idProt +"="+d[idProt])
		print(str(todoBlast)+ " BLASTS TO DO")
		lines = tryGetInfoFrom(idProt,d[idProt],sep)
		saveResultsonCSV(blastCSVName,lines)
		doRandomDelay(minTime,maxTime)
		todoBlast=todoBlast-1
	else:
		print("Blast Done for " + idProt)
	
	




