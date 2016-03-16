#/usr/bin/python

import os.path
import re
import subprocess
from optparse import OptionParser
import glob
import sys
import csv
import numpy as np

directory = os.getcwd()
print(directory)
parser = OptionParser()
parser.add_option("-i", "--input") # Requited
parser.add_option("-c", "--comparison", default='n') # Required
parser.add_option("-t", "--timeseries", default='n') # Required
(options, args) = parser.parse_args()

# 
base_dir = options.input
comparison = options.comparison
timeseries = options.timeseries

print "We are working on traits in " + base_dir
os.chdir(base_dir)

 
def loadphe(infile):
    with open(infile) as f:
        lines = f.read().splitlines()
        return lines
    

# Get a list of raw phenotypes
phe_r = loadphe("phe.r.txt")
nphe_r = phe_r.pop(0)
print str(phe_r) + "\n"

if re.match('y', comparison):
    phe_d = loadphe("phe.d.txt")
    nphe_d = phe_d.pop(0)
    print str(phe_d) + "\n"

final_result = np.empty((0,20), str)

if re.match('y', comparison):
  for p in phe_r:
    #print p
    (treat, trait) = re.split('\.', str(p))
    name_of_summary_table = "summary.table.mqm." + treat + "." + trait + ".csv" 
    path_to_summary_table =  trait + "/" + treat + "/" + "mqm.out" + "/" + name_of_summary_table
    if os.path.exists(path_to_summary_table):
        print "Yes"
        print path_to_summary_table
        reader=csv.reader(open(path_to_summary_table,"r"),delimiter=',')
        x=list(reader)
        result=np.array(x).astype('str')
        d = result.shape
        rownum = d[0]
        print rownum
        exp = trait[-4:-2]
        yr = trait[-2:]
        t1 = [trait] * rownum
        t2 = [treat] * rownum
        t3 = [exp] * rownum
        t4 = [yr] * rownum
        t5 = ["raw"] * rownum
        t1[0] = "trait"
        t2[0] = "treatment"
        t3[0] = "exp"
        t4[0] = "year"
        t5[0] = "type"
        
        print(t1)
        print(type(t1))
        #print(np.shape(np.asarray(t1)))
        
        t1 = np.asarray(t1)
        t1 = np.transpose(t1)
        print(type(t1))
        print(np.shape(t1))
        print(t1)
        
        #print(result[:,5])
        
        print(np.shape(result))
        metadata = np.array((t1, t2, t3, t4, t5))
        metadata = np.transpose(metadata)
        print(np.shape(result))
        print(np.shape(metadata))
        #metadata = np.delete(metadata, 0,0)
        #result = np.delete(result, 0,0)
        result = np.hstack((result, metadata))
        header = result[0,:]
        print(header)
        result = np.delete(result, 0,0)

        
        print(np.shape(result))
        print(type(result))
        #print(result)
        final_result = np.vstack((final_result, result))
        print(np.shape(final_result))
        #print(final_result)

    else:
        print "No"
  for p in phe_d:
    print(p)
    (trait, treat1, treat2) = re.split('\.', str(p))
    name_of_summary_table = "summary.table.mqm." + p + ".csv" 
    path_to_summary_table =  directory + "/" + trait + "/comparison/mqm.out/" + name_of_summary_table
    
    # Get comparison traits
    pwd = os.getcwd()
    print(pwd)
    print(name_of_summary_table)
    print(path_to_summary_table)
    if os.path.exists(path_to_summary_table):
        print "Yes"
        print path_to_summary_table
        reader=csv.reader(open(path_to_summary_table,"r"),delimiter=',')
        x=list(reader)
        result=np.array(x).astype('str')
        d = result.shape
        rownum = d[0]
        print rownum
        exp = trait[-4:-2]
        yr = trait[-2:]
        t1 = [trait] * rownum
        t2 = ["diff"] * rownum
        t3 = [exp] * rownum
        t4 = [yr] * rownum
        t5 = ["comp"] * rownum
        t1[0] = "trait"
        t2[0] = "treatment"
        t3[0] = "exp"
        t4[0] = "year"
        t5[0] = "type"
        
        print(t1)
        print(type(t1))
        #print(np.shape(np.asarray(t1)))
        
        t1 = np.asarray(t1)
        t1 = np.transpose(t1)
        print(type(t1))
        print(np.shape(t1))
        print(t1)
        
        #print(result[:,5])
        
        print(np.shape(result))
        metadata = np.array((t1, t2, t3, t4, t5))
        metadata = np.transpose(metadata)
        print(np.shape(result))
        print(np.shape(metadata))
        #metadata = np.delete(metadata, 0,0)
        #result = np.delete(result, 0,0)
        result = np.hstack((result, metadata))
        header = result[0,:]
        print(header)
        result = np.delete(result, 0,0)

        
        print(np.shape(result))
        print(type(result))
        #print(result)
        final_result = np.vstack((final_result, result))
        print(np.shape(final_result))
        #print(final_result)
    else:
        print "No"
        
if re.match('n', comparison):
    for p in phe_r:
      print(p)
    name_of_summary_table = "summary.table.mqm." + p + ".csv" 
    path_to_summary_table =  trait + "/" + "comparison" + "/" + "mqm.out" + "/" + name_of_summary_table
    print(name_of_summary_table)
    if os.path.exists(path_to_summary_table):
        print "Yes"
        print path_to_summary_table
        reader=csv.reader(open(path_to_summary_table,"r"),delimiter=',')
        x=list(reader)
        result=np.array(x).astype('str')
        d = result.shape
        rownum = d[0]
        print rownum
        exp = trait[-4:-2]
        yr = trait[-2:]
        t1 = [trait] * rownum
        t2 = ["none"] * rownum
        t3 = [exp] * rownum
        t4 = [yr] * rownum
        t5 = ["raw"] * rownum
        t1[0] = "trait"
        t2[0] = "treatment"
        t3[0] = "exp"
        t4[0] = "year"
        t5[0] = "type"
        
        print(t1)
        print(type(t1))
        #print(np.shape(np.asarray(t1)))
        
        t1 = np.asarray(t1)
        t1 = np.transpose(t1)
        print(type(t1))
        print(np.shape(t1))
        print(t1)
        
        #print(result[:,5])
        
        print(np.shape(result))
        metadata = np.array((t1, t2, t3, t4, t5))
        metadata = np.transpose(metadata)
        print(np.shape(result))
        print(np.shape(metadata))
        #metadata = np.delete(metadata, 0,0)
        #result = np.delete(result, 0,0)
        result = np.hstack((result, metadata))
        header = result[0,:]
        print(header)
        result = np.delete(result, 0,0)

        
        print(np.shape(result))
        print(type(result))
        #print(result)
        final_result = np.vstack((final_result, result))
        print(np.shape(final_result))
        #print(final_result)
    else:
        print "No"
   
concatenated_summary_table= base_dir + "_" + "concatenated_summary_table.csv"
np.savetxt(concatenated_summary_table, final_result, delimiter=",", fmt="%s")

os.chdir(directory)
banan_banner_command = "Rscript mk.banan.banner.R " + base_dir + " " + comparison  
subprocess.call(banan_banner_command, shell=True)


if re.match('y', timeseries):
    if re.match('y', comparison):
        print("Working on the timeseries qtl")
        print("In base directory:")
        print(base_dir)
        print("Treatments are:")
        print(treat1)
        print(treat2)
        time_series_qtl_command = "Rscript time_series_qtl_comp.R " + base_dir + " " + treat1 + " " + treat2
        subprocess.call(time_series_qtl_command, shell=True)
    if re.match('n', comparison):
        print("Timeseries QTL analysis in a single treatment is not yet supported.")

