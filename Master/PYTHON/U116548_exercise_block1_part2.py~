def get_sequence_from_FASTA_file(filename,residue,threshold=0.3):
	fd= open(filename)
	a=0
	residue=residue.upper()
	for line in fd: 
		if not line.startswith(">"):
			if float(line.count(residue))/len(line) >= threshold:
				a+=1 
	print a
	return a
	

def print_sequence_tails(filename,first_n=10,last_m=10):
	fd= open(filename)
	for line in fd:
		if line.startswith(">"):
			id=line.strip("\n")
			
		else:
			print (id+"\t"+line[:first_n]+"\t"+line[(-last_m-1):-1])





get_sequence_from_FASTA_file("sample_fasta1.fa","k",0.000005) 

print_sequence_tails("sample_fasta1.fa",5,5)
