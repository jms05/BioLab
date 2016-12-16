from Bio import SeqIO


def parseFile(filename):
	seq_record = SeqIO.read(filename, "genbank")
	return seq_record


rec = parseFile("grupo6.txt")

i=0
for feat in rec.features:
	if feat.type == 'CDS':
		i+=1
		print ("xref -> " + feat.qualifiers["db_xref"][0])
		print ("ID_prot -> " + feat.qualifiers["protein_id"][0])
		try:
			print ("Funcao-> " + feat.qualifiers["function"][0])
		except Exception:
			print ("Funcao-> Unknown")
		
		print ("Traducao -> " + feat.qualifiers["translation"][0])
		print("-"*20)

print ("Total: " + str(i))