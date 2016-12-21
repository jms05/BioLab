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



def getInfouniprot(protainID):
	params = {"query": protainID, "format": "fasta"}
	r = requests.get("http://www.uniprot.org/uniprot/", params)
	i=0
	for record in SeqIO.parse(StringIO(r.text), "fasta"):
		idfasta = record.id
		idfasta = str(idfasta).split("|")[1]
		#break

	return idfasta
	parst = parseSwissProt(idfasta)
	return(idfasta,parst)

a = "P23687.1"

print (str(getInfouniprot(a)))