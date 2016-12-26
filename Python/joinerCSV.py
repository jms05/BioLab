import os
import random

def getIndex(dic,key):
	return dic[key]

def getinfoFromBlastCSV(filename):
	dicheadpos ={}
	i=0
	f = open(filename,"r")
	filecsvlines = f.read().splitlines()
	f.close()
	sep = filecsvlines[0].split("=")[1]
	heads = filecsvlines[1].split(sep)
	for head in heads:
		dicheadpos[head]=i
		i+=1
	csvRecors= filecsvlines[2:]
	parseddata = {}
	for record in csvRecors:
		campos = record.split(sep)
		idprot = campos[dicheadpos["Idprot"]]
		idmatch = campos[dicheadpos["idBlastMatch"]]
		if(idmatch!="NoResult"):
			evalu  = campos[dicheadpos["evalu"]]
			score =  campos[dicheadpos["Score"]]
			protrNames  = campos[dicheadpos["protName"]].split("|")
			protStatus = campos[dicheadpos["protStatus"]]
			grauRev = campos[dicheadpos["grauRev"]]
			host = campos[dicheadpos["host"]].replace("()","")
			local  =campos[dicheadpos["protLocal"]].replace("-","Unknown")
			functions = campos[dicheadpos["funcaoMolec"]]
			if functions == "-":
				functions=[]
				functions.append("Unknown")
			else:
				functions = functions.split("|")

			processBio = campos[dicheadpos["processBiol"]]
			if processBio == "-":
				processBio=[]
				processBio.append("Unknown")
			else:
				processBio = processBio.split("|")
			genName  =campos[dicheadpos["geneNameSwiss"]].replace("-","Unknown")

			domain = campos[dicheadpos["domianTex"]]
			if domain == "-":
				domain=[]
				domain.append("Unknown")
			else:
				domain = domain.split("|")
			seq = campos[dicheadpos["sequence"]]
			parseddata[idprot] = (idmatch,evalu,score,protrNames,protStatus,grauRev,host,local,functions,processBio,genName,domain,seq)
		else:
			parseddata[idprot] = (idmatch)

	for head in dicheadpos:
		dicheadpos[head]-=1

	return (dicheadpos,parseddata)

def getinfoFromOriginal(filename):
	dicheadpos ={}
	i=0
	f = open(filename,"r")
	filecsvlines = f.read().splitlines()
	f.close()
	sep = filecsvlines[0].split("=")[1]
	heads = filecsvlines[1].split(sep)
	for head in heads:
		dicheadpos[head]=i
		i+=1
	csvRecors= filecsvlines[2:]
	parseddataNeedBlast = {}
	parseddataNotNeedBlast = {}

	for record in csvRecors:
		campos = record.split(sep)
		geneID = campos[dicheadpos["geneID"]]
		geneName = campos[dicheadpos["GeneName"]].replace("-","Unknown")
		geneAccNUM = campos[dicheadpos["GeneAccessNumber"]].replace("-","Unknown")
		locusTag = campos[dicheadpos["locusTag"]]
		strand = campos[dicheadpos["strand"]]
		dnaSeq = campos[dicheadpos["DNA_SEQ"]]
		AccNumNCBI  =campos[dicheadpos["AccessNumberNCBI"]]
		idSwiss = campos[dicheadpos["idSwiss"]]
		protName = campos[dicheadpos["protName"]]
		aaSeq = campos[dicheadpos["PROT_SEQ"]]
		protLen = campos[dicheadpos["PROT_Tamanho"]]
		status  =campos[dicheadpos["protStatus"]]
		grauRev  = campos[dicheadpos["grauRev"]]
		local  =campos[dicheadpos["protLocal"]].replace("-","Unknown")
		ec  =campos[dicheadpos["EC"]].replace("-","Unknown")

		functionsGO = campos[dicheadpos["(GO)funcaoMolec"]]
		if functionsGO == "-":
			functionsGO=[]
			functionsGO.append("Unknown")
		else:
			functionsGO = functionsGO.split("_")

		processBioGO = campos[dicheadpos["(GO)processBiol"]]
		if processBioGO == "-":
			processBioGO=[]
			processBioGO.append("Unknown")
		else:
			processBioGO = processBioGO.split("_")

		functions = campos[dicheadpos["funcoes"]]
		functions = functions.split("_")

		domain = campos[dicheadpos["Domain"]]
		if domain == "-":
			domain=[]
			domain.append("Unknown")
		else:
			domain = domain.split("|")
		##Alldata reed
		data  =(geneName,geneAccNUM,locusTag,strand,dnaSeq,AccNumNCBI,idSwiss,protName,aaSeq,protLen,status,grauRev,local,ec,functionsGO,processBioGO,functions,domain)
		if(functions[0]=="Unknown"): #TODOBLAST
			parseddataNeedBlast[geneID]=data
		else: ## not need blast
			parseddataNotNeedBlast[geneID]=data

	for head in dicheadpos:
		dicheadpos[head]-=1
	return (dicheadpos,parseddataNeedBlast,parseddataNotNeedBlast)


def displaySEQ(seq,maxonline=65,tablen=2):
	ret=""
	i=0
	for base in seq:
		if(i>=maxonline):
			ret+="\n"+"\t"*tablen
			i=0
		ret+=base
		i+=1
	return ret

def stringfromLis(lista,tabs=3,maxlin=5,sep=","):
	funs = ""
	i=0
	for fun in lista:
		if(i>=maxlin):
			i=0
			funs+="\n"+"\t"*tabs
		funs+=(" "+fun+",")
		i+=1
	funs=funs[:-1]
	return funs

def createFilenoBlast(geneID,geneName,geneAccNUM,locusTag,
		strand,dnaSeq,AccNumNCBI,idSwiss,protName,aaSeq,protLen,status,
		grauRev,local,ec,functionsGO,processBioGO,functions,domain):
	filename = "gene_NoBlast_" + geneID+".txt"
	f = open(filename,"w")
	f.write("Protain Name: " + protName +"\t"*3 + "Uniprot ID: " + idSwiss +"\n\n")
	funs = stringfromLis(functions,maxlin=1)
	f.write("Functions:" + funs+ "\n\n")
	f.write("Gene Ontology:\n\n")
	funs = stringfromLis(functionsGO,maxlin=1,tabs=7)
	f.write("\t"*2+"Nuclear Function: " + funs + "\n\n")
	funs = stringfromLis(processBioGO,maxlin=1,tabs=7)
	f.write("\t"*2+"Biological Process: " + funs + "\n\n")
	f.write("Location: " + local+"\n\n")
	f.write("Revision: "+ status+ " ("+grauRev+"/5)\n\n" )
	prot = displaySEQ( aaSeq)
	f.write("Seq AA ("+str(protLen)+"): " + prot+"\n\n")
	f.close()
	return

def createBlastHead(message="Blast Results",carct="-",times=32):
	return carct*times+message+carct*times

def createFileBlastNotFound(geneID,geneName,geneAccNUM,locusTag,
		strand,dnaSeq,AccNumNCBI,idSwiss,protName,aaSeq,protLen,status,
		grauRev,local,ec,functionsGO,processBioGO,functions,domain):
	filename = "gene_Blast_Not_Found_" + geneID+".txt"
	f = open(filename,"w")
	f.write("Protain Name: " + protName +"\t"*3 + "Uniprot ID: " + idSwiss +"\n\n")
	funs = stringfromLis(functions,maxlin=1)
	f.write("Functions:" + funs+ "\n\n")
	f.write("Gene Ontology:\n\n")
	funs = stringfromLis(functionsGO,maxlin=1,tabs=7)
	f.write("\t"*2+"Nuclear Function: " + funs + "\n\n")
	funs = stringfromLis(processBioGO,maxlin=1,tabs=7)
	f.write("\t"*2+"Biological Process: " + funs + "\n\n")
	f.write("Location: " + local+"\n\n")
	f.write("Revision: "+ status+ " ("+grauRev+"/5)\n\n" )
	prot = displaySEQ( aaSeq)
	f.write("Seq AA ("+str(protLen)+"): " + prot+"\n")
	f.write("\n"+createBlastHead()+"\n\n")
	f.write("Result: NOT FOUND\n")
	f.close()
	return

def createFileBlastFound(geneID,geneName,geneAccNUM,locusTag,
		strand,dnaSeq,AccNumNCBI,idSwiss,protName,aaSeq,protLen,status,
		grauRev,local,ec,functionsGO,processBioGO,functions,domain,
		idBlastMatch,evalu,score,protrNamesBlast,protStatusBlast,grauRevBlast,
		hostBlast,localBlast,functionsBlast,processBioBlast,genNameBlast,domainBlast,seqBlast):
	filename = "gene_Blast_Found_" + geneID+".txt"
	f = open(filename,"w")
	f.write("Protain Name: " + protName +"\t"*3 + "Uniprot ID: " + idSwiss +"\n\n")
	funs = stringfromLis(functions,maxlin=1)
	f.write("Functions:" + funs+ "\n")
	f.write("Gene Ontology:\n\n")
	funs = stringfromLis(functionsGO,maxlin=1,tabs=7)
	f.write("\t"*2+"Nuclear Function: " + funs + "\n\n")
	funs = stringfromLis(processBioGO,maxlin=1,tabs=7)
	f.write("\t"*2+"Biological Process: " + funs + "\n\n")
	f.write("Location: " + local+"\n\n")
	f.write("Revision: "+ status+ " ("+grauRev+"/5)\n\n" )
	prot = displaySEQ( aaSeq)
	f.write("Seq AA ("+str(protLen)+"): " + prot+"\n")
	f.write("\n"+createBlastHead()+"\n")
	f.write("Result: "+ idBlastMatch + "\t"+ "E-Value: "+ evalu + "\t"+ "Score: "+ score+"\n\n")
	f.write("Host: " + hostBlast + "\n\n")
	protrNamesBlast1 = stringfromLis(protrNamesBlast,maxlin=1,tabs=7)
	f.write("Protain Name: " + protrNamesBlast1 +"\n\n")
	funs = stringfromLis(functionsBlast,maxlin=1)
	f.write("Functions:" + funs+ "\n\n")
	f.write("Gene Ontology:\n\n")
	funs = stringfromLis(processBioBlast,maxlin=1,tabs=7)
	f.write("\t"*2+"Biological Process: " + funs + "\n\n")
	f.write("Location: " + localBlast+"\n\n")
	f.write("Revision: "+ protStatusBlast+ " ("+grauRevBlast+"/5)\n\n" )
	prot = displaySEQ( seqBlast)
	f.write("Seq AA ("+str(protLen)+"): " + prot+"\n")
	f.close()
	return


def processBlast(dataOriginal,blastData):
	for geneID in dataOriginal.keys():
		print("Crating file for: " + geneID)
		(geneName,geneAccNUM,locusTag,strand,dnaSeq,AccNumNCBI,idSwiss,protName,aaSeq,protLen,status,grauRev,local,ec,functionsGO,processBioGO,functions,domain) = dataOriginal[geneID]
		blastDatagene = blastData[geneID]
		#print("--"+str(blastDatagene)+"--")
		if(blastDatagene=="NoResult"):
			print("NO BLAST FOUND")
			createFileBlastNotFound(geneID,geneName,geneAccNUM,locusTag,strand,dnaSeq,AccNumNCBI,idSwiss,protName,aaSeq,
				protLen,status,grauRev,local,ec,functionsGO,processBioGO,functions,domain)
		else:
			print("BLAST FOUND")
			(idmatch,evalu,score,protrNames,protStatus,grauRev,host,local,functions,processBio,genName,domain,seq) = blastDatagene
			createFileBlastFound(geneID,geneName,geneAccNUM,locusTag,strand,dnaSeq,AccNumNCBI,idSwiss,protName,aaSeq,
				protLen,status,grauRev,local,ec,functionsGO,processBioGO,functions,domain,
				idmatch,evalu,score,protrNames,protStatus,grauRev,host,local,functions,processBio,genName,domain,seq)
		print("File created for "+ geneID)

def processNoblast(data):
	for geneID in data.keys():
		print("Crating file for: " + geneID)
		(geneName,geneAccNUM,locusTag,strand,dnaSeq,AccNumNCBI,idSwiss,protName,aaSeq,protLen,status,grauRev,local,ec,functionsGO,processBioGO,functions,domain) = data[geneID]
		createFilenoBlast(geneID,geneName,geneAccNUM,locusTag,strand,dnaSeq,AccNumNCBI,idSwiss,protName,aaSeq,protLen,
			status,grauRev,local,ec,functionsGO,processBioGO,functions,domain)
		print("File created for "+ geneID)


def joinData(originalData,blastD):
	(dicOrig,blasted,noblast) = originalData
	(dicBlast,blastData)=blastD
	processNoblast(noblast)
	processBlast(blasted,blastData)

originalData  =getinfoFromOriginal("tabela_com_seq_Completa.csv")
blastData = getinfoFromBlastCSV("BlastCSV.csv")

joinData(originalData,blastData)
