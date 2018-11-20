#!/usr/bin/env python
from Bio import pairwise2
from Bio import Align
from Bio import SeqIO

# ind_system = SeqIO.to_dict(SeqIO.parse("/home/lucas/ISGlobal/Chip_Seq_Oriol/Fastas_insercions/insercions_3E_F1.fa", "fasta"))
# reference_genome = SeqIO.to_dict(SeqIO.parse("/home/lucas/ISGlobal/Chip_Seq_Oriol/Reference_Genomes/PlasmoDB-39_Pfalciparum3D7_Genome_3E_F1.fasta", "fasta"))

ind_system = SeqIO.to_dict(SeqIO.parse("/home/lucas/ISGlobal/Chip_Seq_Oriol/Fastas_insercions/insercions_induits.fa", "fasta"))
reference_genome = SeqIO.to_dict(SeqIO.parse("/home/lucas/ISGlobal/Chip_Seq_Oriol/Reference_Genomes/PlasmoDB-39_Pfalciparum3D7_Genome_Ind_System_In.fasta", "fasta"))

# ind_system = SeqIO.to_dict(SeqIO.parse("/home/lucas/ISGlobal/Chip_Seq_Oriol/Fastas_insercions/insercions_NOinduits.fa", "fasta"))
# reference_genome = SeqIO.to_dict(SeqIO.parse("/home/lucas/ISGlobal/Chip_Seq_Oriol/Reference_Genomes/PlasmoDB-39_Pfalciparum3D7_Genome_Ind_System_NI.fasta", "fasta"))

pf3D7 = SeqIO.to_dict(SeqIO.parse("/home/lucas/ISGlobal/Gen_Referencies/PlasmoDB-31_Pfalciparum3D7_Genome.fasta", "fasta"))

## Bowtie command
## "bowtie2 -f -x Ind_System_In -U ../Fastas_insercions/insercions_induits.fa  > ../bowtie2_aligned_inducedsystem.sam"
#
# with open("/home/lucas/ISGlobal/Chip_Seq_Oriol/bowtie2_aligned_3E_F1.sam") as samfile:
#     for line in samfile:
#         if line.startswith("@"):
#             pass
#         else:
#             if line.split()[2] in reference_genome.keys():
#                 print line.split()[0]
#                 print line.split()[2]
#                 print ind_system[line.split()[0]].seq in reference_genome[line.split()[2]].seq
#                 print "From: ", line.split()[3], " to ", int(line.split()[3])+len(ind_system[line.split()[0]].seq)
#                 print ""
#                 print ">", line.split()[0]
#                 print ind_system[line.split()[0]].seq
#                 print ">", line.split()[2]
#                 print reference_genome[line.split()[2]][int(line.split()[3])-1:int(line.split()[3])+len(ind_system[line.split()[0]].seq)-1].seq
#                 print "---------------------------\n"
#             else:
#                 print "No hit for {}" .format(line.split()[0])

# alignment = pairwise2.align.globalxx(ind_system["pL7ap2-g_induit"], reference_genome["Pf3D7_12_v3"][907042:907042+len(ind_system["pL7ap2-g_induit"])].seq)
#
# print(pairwise2.format_alignment(*alignment[0]))
# print alignment[2:]

# print ">Modified"
# print reference_genome["Pf3D7_12_v3"][897203:924501].seq
# print ">Reference_3D7"
# print pf3D7["Pf3D7_12_v3"][897203:924501].seq


print "Modified Genome\n"
for key, value in reference_genome.iteritems():
    print key, len(value.seq)
print ""

print "Inserts\n"
for key, value in ind_system.iteritems():
    print key, len(value.seq)
print ""

print "Reference 3D7 Genome\n"
for key, value in pf3D7.iteritems():
    print key, len(value.seq)
