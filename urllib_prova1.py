#!/usr/bin/env python

import urllib
import urllib2
 
data = urllib.urlencode({'seqfile': '/home/lucas/ISGlobal/Cruzi/tcruzi_epitopes_vaccine/Gen_referencies/PDB_fastas/fasta.fasta', 'outform':'-short'})
url = 'http://www.cbs.dtu.dk/services/TMHMM/cgi-bin/webface2.fcgi'
full_url = url + '?' + data
response = urllib2.urlopen(full_url)

print response.read.decode() 