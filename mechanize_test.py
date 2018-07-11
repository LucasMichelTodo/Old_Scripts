#!/usr/bin/env python

import mechanize 

# Browser
br = mechanize.Browser()

# Ignore robots.txt
br.set_handle_robots( False )
# Google demands a user-agent that isn't a robot
br.addheaders = [('User-agent', 'Firefox')]

# Retrieve the Google home page, saving the response
br.open( "http://www.cbs.dtu.dk/services/SignalP/" )

# Select the search box and search for 'foo'
for form in br.forms():
	print "Form name:", form.name
	print form

br.select_form(nr = 0)
br.form[ 'SEQPASTE' ] = 'AETEVSRDPTDEHSD'


response = br.submit()
print response.read()
