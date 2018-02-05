#!/usr/bin/env python

while True:
    inp = raw_input('Enter Fahrenheit Temperature or \"quit\":')
    if inp == "quit":
        break
    try:
        fahr = float(inp)
        cel = (fahr - 32.0) * 5.0 / 9.0
        print cel
    except:
        print 'Please enter a number'
        continue
