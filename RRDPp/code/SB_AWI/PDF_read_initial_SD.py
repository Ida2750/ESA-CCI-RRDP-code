# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 16:01:15 2023

@author: Ida Olsen

Read PDF file 
"""

# importing required modules
import PyPDF2
import os
import re
import numpy as np

# creating a pdf file object
directory = "C:/Users/Ida Olsen/Documents/work/RRDPp/RawData/SDB_AWI/Arctic/DATA_Past_Buoys/2019S93_data"

def pdf_read_initial(file):
    pdfFileObj = open(file, 'rb')
      
    # creating a pdf reader object
    pdfReader = PyPDF2.PdfReader(pdfFileObj)
      
    # creating a page object - get first page
    pageObj = pdfReader.pages[0]
    
    # extract text from page
    text = pageObj.extract_text()
    
    # substitue special characters
    text = re.sub('[^A-Za-z0-9.]+', ' ', text)
    
    # search for initial sit
    sit = re.search('Ice thickness ', text)
    SIT_init = text[sit.span()[1]:sit.span()[1]+4].replace('m','')
    try:    
        SIT_init=float(SIT_init.strip())
    except:  # no initial sea ice thicknes recorded
        SIT_init = np.nan
    # search for initial snowdepths on PDF
    try:
        search1 = re.search('Initial sh1 ', text)
        sh1 = text[search1.span()[1]:search1.span()[1]+4].replace('m','')
    except:
        sh1=np.nan
    try:
        search2 = re.search('Initial sh2 ', text)
        sh2 = text[search2.span()[1]:search2.span()[1]+4].replace('m','')
    except:
        sh2=np.nan
    try:
        search3 = re.search('Initial sh3 ', text)
        sh3 = text[search3.span()[1]:search3.span()[1]+4].replace('m','')
    except:
        sh3=np.nan
    try:
        search4 = re.search('Initial sh4 ', text)
        sh4 = text[search4.span()[1]:search4.span()[1]+4].replace('m','')
    except:
        sh4=np.nan
        
    SD_init = np.nanmean([float(sh1), float(sh2), float(sh3), float(sh4)])  
    # closing the pdf file object
    pdfFileObj.close()
    
    return [float(sh1), float(sh2), float(sh3), float(sh4), SD_init, SIT_init]