#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 28 16:07:41 2018

@author: lucas
"""

def split(a, n):
	k, m = divmod(len(a), n)
	splits = [a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in xrange(n)]
	return [[i[0], i[-1]] for i in splits]


splits = split(range(20,31), 3)

print splits

   
coverages = "\t".join(str(i) for i in sum(splits, []))

print coverages