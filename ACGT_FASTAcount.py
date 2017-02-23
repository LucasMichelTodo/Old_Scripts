from tqdm import tqdm
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

def get_sequence_from_FASTA_file(filename):
	fd= open(filename)
	seq=""
	for line in tqdm(fd): 
		if not line.startswith(">"):
				seq = seq + Seq(line.strip("\n"))
	print seq
	print len(seq)
	print seq.count("A")/float(len(seq))*100
	print seq.count("T")/float(len(seq))*100
	print seq.count("G")/float(len(seq))*100
	print seq.count("C")/float(len(seq))*100
	return seq

get_sequence_from_FASTA_file("/home/lucas/ISGlobal/Scripts/Master/PYTHON/sample_fasta1.fa")
	

# def print_sequence_tails(filename,first_n=10,last_m=10):
# 	fd= open(filename)
# 	for line in fd:
# 		if line.startswith(">"):
# 			id=line.strip("\n")
			
# 		else:
# 			print (id+"\t"+line[:first_n]+"\t"+line[(-last_m-1):-1])


