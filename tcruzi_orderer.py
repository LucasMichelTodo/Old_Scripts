#!/usr/bin/env python

import subprocess as sp

cmd = "cd /media/lucas/Disc4T/Projects/tcruzi_Actual/Reference_fastas"
sp.Popen(cmd, shell=True).wait()

cmd = "fasta_tagger.py TriTrypDB-35_TcruziSylvioX10-1_ORFs_AA.fasta 01 > TcruziSylvioX10-1_tagged.fa"
sp.Popen(cmd, shell=True).wait()
cmd = "fasta_tagger.py TriTrypDB-38_TcruziSylvioX10-1-2012_ORFs_AA.fasta 02 > TcruziSylvioX10-1-2012_tagged.fa"
sp.Popen(cmd, shell=True).wait()
cmd = "fasta_tagger.py TriTrypDB-35_TcruziJRcl4_ORFs_AA.fasta 03 > TcruziJRcl4.fa"
sp.Popen(cmd, shell=True).wait()
cmd = "fasta_tagger.py TriTrypDB-35_TcruziTulacl2_ORFs_AA.fasta 04 > Tulacl2_tagged.fa"
sp.Popen(cmd, shell=True).wait()
cmd = "fasta_tagger.py TriTrypDB-35_TcruziEsmeraldo_ORFs_AA.fasta 05 > TcruziEsmeraldo_tagged.fa"
sp.Popen(cmd, shell=True).wait()
cmd = "fasta_tagger.py TriTrypDB-35_TcruziDm28c_ORFs_AA.fasta 06 > TcruziDm28c_tagged.fa"
sp.Popen(cmd, shell=True).wait()
cmd = "fasta_tagger.py TriTrypDB-38_TcruzicruziDm28c_ORFs_AA.fasta 07 > TcruzicruziDm28c_tagged.fa"
sp.Popen(cmd, shell=True).wait()
cmd = "fasta_tagger.py Tcruzi_TCC.fasta 08 > Tcruzi_TCC_tagged.fa"
sp.Popen(cmd, shell=True).wait()
cmd = "fasta_tagger.py TriTrypDB-35_TcruziCLBrener_ORFs_AA.fasta 09 > TcruziCLBrener_tagged.fa"
sp.Popen(cmd, shell=True).wait()
cmd = "fasta_tagger.py TriTrypDB-35_TcruziCLBrenerEsmeraldo-like_ORFs_AA.fasta 10 > TcruziCLBrenerEsmeraldo-like_tagged.fa"
sp.Popen(cmd, shell=True).wait()
cmd = "fasta_tagger.py TriTrypDB-35_TcruziCLBrenerNon-Esmeraldo-like_ORFs_AA.fasta 11 > TcruziCLBrenerNon-Esmeraldo-liketagged.fa"
sp.Popen(cmd, shell=True).wait()
cmd = "fasta_tagger.py Tcruzi_TcIII_strain231_predicted_proteins.fasta 12 > Tcruzi_TcIII_strain231_predicted_tagged.fa"
sp.Popen(cmd, shell=True).wait()
cmd = "fasta_tagger.py ../exposed_taged.fa 10 > exposed_taged_tagged.fa"
sp.Popen(cmd, shell=True).wait()

cmd = "mkdir Tagged_fastas"
sp.Popen(cmd, shell=True).wait()

cmd = "mv ./*\_tagged.fa Tagged_fastas"
sp.Popen(cmd, shell=True).wait()
