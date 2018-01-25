#!/usr/bin/env python

import pandas as pd
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
import re
from tqdm import tqdm

epitopes_df = pd.read_csv("/home/lucas/ISGlobal/Cruzi/tcruzi_epitopes_vaccine/epitopes.txt", sep='\t', header = None, skiprows = 373)
epitopes_df.columns = ["EpSeq", "length", "ID", "Ag", "AgID", "Start", "End"]

epitopes_df["SP"] = ""
epitopes_df["TM"] = ""
epitopes_df["TP"] = ""

i = 0

for epitope in tqdm(epitopes_df.iloc[:,0]):

	driver = webdriver.Firefox()
	driver2 = webdriver.Firefox()
	driver3 = webdriver.Firefox()

	driver.get("http://www.cbs.dtu.dk/services/SignalP/")
	driver2.get("http://www.cbs.dtu.dk/services/TMHMM/")
	driver3.get("http://www.cbs.dtu.dk/services/TargetP/")

	element = driver.find_element_by_name("SEQPASTE")
	element2 = driver2.find_element_by_name("SEQ")
	element3 = driver3.find_element_by_name("SEQPASTE")

	element.send_keys(epitope)
	element2.send_keys(epitope)
	element3.send_keys(epitope)

	element.submit()
	element2.submit()
	element3.submit()

	WebDriverWait(driver, 3600).until(EC.title_contains("P 4.1 Server - prediction results"))  
	WebDriverWait(driver2, 3600).until(EC.title_contains("TMHMM result"))
	WebDriverWait(driver3, 3600).until(EC.title_contains("TargetP 1.1 Server - prediction results"))

	predictionSP = driver.page_source
	predictionTM = driver2.page_source
	predictionTP = driver3.page_source

	regexSP = re.compile("SP='\w*'")
	regexTM = re.compile("predicted TMHs:\s*\d+")
	regexTP = re.compile("Sequence.*")

	epitopes_df.loc[i,"SP"] = re.findall(regexSP, predictionSP)[0].replace("SP=", "").replace("'", "")
	epitopes_df.loc[i,"TM"] = re.findall(regexTM, predictionTM)[0].replace("predicted TMHs:", "")
	epitopes_df.loc[i,"TP"] = re.findall(regexTP, predictionTP)[0].split()[5]


	driver.quit()
	driver2.quit()
	driver3.quit()

	epitopes_df.iloc[0:i+1,:].to_csv(path_or_buf="/home/lucas/ISGlobal/Cruzi/tcruzi_epitopes_vaccine/epitopes_table_progress.csv", index = False)

	i += 1




