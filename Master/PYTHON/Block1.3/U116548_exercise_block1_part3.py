def calculate_aminoacid_frequencies(fasta_filename, subsequences_filename):
    fasta = open(fasta_filename)
    fragments = open(subsequences_filename)
    i=0
    for line in fd:
        if line[0]==">":
            i +=1
    print(i)
    print(string(fragments.countline))
    fasta.close()
    fragments.close()

calculate_aminoacid_frequencies("uniprot_sprot_short.fasta")
