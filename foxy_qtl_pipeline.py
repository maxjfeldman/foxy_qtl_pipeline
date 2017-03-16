#/usr/bin/python

# import librariers
import sys
from optparse import OptionParser
import os.path
import subprocess
import re
import multiprocessing as mp
import glob

# Define variables used in the program
version = "1.10"  # Version of the pipeline
pwd = os.getcwd()

# Get arguments
parser = OptionParser()
parser.add_option("-i", "--input") # Requited
parser.add_option("-o", "--output") # Required
parser.add_option("-c", "--comparison") # Required
parser.add_option("-m", "--map") # Optional
parser.add_option("-s", "--server") # Optional (Are you on a server?)
parser.add_option("-d", "--model", default = "normal") # Optional(what type of distribution describes the trait?)
parser.add_option("-q", "--qtl", default= "hk") # Optional (what is the method you want to use for qtl analysis?)
parser.add_option("-t", "--type", default= "riself") # Optional(what type of cross is this?)

(options, args) = parser.parse_args()


# Print a statement that tells user the program is running and version of pipeline
print "\nYou are running foxy_qtl_pipeline.py pipeline program version " + version + "\n"

# Re-assign values from command line into more intuitive variable names
infile = options.input
outfile = options.output
comparison = options.comparison
g_map = options.map
server = options.server
trait_model = options.model
qtl_method = options.qtl
cross_type = options.type

# Can't do a comparison with a binary trait
if trait_model == 'binary':
    comparison = 'n'


# If no directory to store data is present, go ahead and create it
if not os.path.exists(pwd+"/"+outfile):
    os.makedirs(pwd+"/"+outfile)

# If the infile exists run the exploratory data analysis script in R
if os.path.exists(pwd+"/"+infile):
    print "\nYour input file " + infile + " was found. Proceeding with EDA...\n\n"
    run_EDA_command = "Rscript QTL_format_and_EDA.R " + infile + " " + outfile
    #run_EDA_command = "echo Hello World"
    print run_EDA_command
    subprocess.call(run_EDA_command, shell=True)

print "Checking for phenotype input file congruent with the genetic map...\nIf there is one present it is named: 'qtl.phenofile.csv'...\nIf none are present one will be created...\n\n";

# Prepare genetic map, only keeping individuals with marker data, build a set of R/qtl cross objects (individual treatments, difference between treatments or plots and BLUP models)
print "The name of the map is: " + g_map + "\n"

run_reformat_command_raw = "Rscript format_and_mk.cross.obj.raw.R " + g_map + " " + outfile + " " + cross_type 
#run_reformat_command_raw = "echo reformat_raw command here";
subprocess.call(run_reformat_command_raw, shell=True)

if re.match('y', comparison):
    run_reformat_command_diff = "Rscript format_and_mk.cross.obj.diff.R " + g_map + " " + outfile + " " + cross_type
    #run_reformat_command_diff = "echo reformat_diff command here";
    subprocess.call(run_reformat_command_diff, shell=True)
else:
    print "Command-line arguments indicate no comparison was specified.\n"

def loadphe(infile):
    with open(infile) as f:
        lines = f.read().splitlines()
        return lines

# Get a list of raw phenotypes
phe_r = loadphe(outfile + "/phe.r.txt")
nphe_r = phe_r.pop(0)
print str(phe_r) + "\n"

# If comparison is desired get list of phenotypes that are the difference between treatments
if re.match('y', comparison):
    phe_d = loadphe(outfile + "/phe.d.txt")
    nphe_d = phe_d.pop(0)


cond=dict()
if re.match('y', comparison):
    for p in phe_r:
        (treat, trait) = re.split('\.', str(p))
        if treat in cond:    
            cond[treat].append(trait)
        else:
            cond[treat] = [trait]
    # Get a list of the keys in 
    comp = cond.keys()
    for k in comp:
        print k + "\t" + str(cond.get(k)) + "\n"
        print "This is a key: " + k + "\n"
    print "You have indicated you'd like to do a comparison...\n"
    print "Checking for directories to store the analysis results in...\nIf not found they will be created...\n";
    
    # Get a list of traits (it is important that all traits be present in all treatments)
    traits = cond.values()[0]
    
    for t in traits:
        s = "/"
        seq = (pwd, outfile, t)
        trait_dir = s.join(seq)
        if not os.path.exists(s.join(seq)):
            print s.join(seq)
            os.makedirs(s.join(seq))
            
            s = "/"
            seq = (pwd, outfile, t, 'comparison')
            print s.join(seq)
            os.makedirs(s.join(seq))
        
            s = "/"
            seq = (pwd, outfile, t, 'comparison', 'scanone.out')
            print s.join(seq)
            os.makedirs(s.join(seq))
            
            s = "/"
            seq = (pwd, outfile, t, 'comparison', 'mqm.out')
            print s.join(seq)
            os.makedirs(s.join(seq))

            for c in comp:
                cond_trait_dir = trait_dir + "/" + c 
                print cond_trait_dir
                if not os.path.exists(cond_trait_dir):
                    os.makedirs(cond_trait_dir)
                cond_trait_dir_so = trait_dir + "/" + c + "/" + "scanone.out"
                print cond_trait_dir_so
                if not os.path.exists(cond_trait_dir_so):
                    os.makedirs(cond_trait_dir_so)
                 
                cond_trait_dir_mqm = trait_dir + "/" + c + "/" + "mqm.out"
                print cond_trait_dir_mqm
                if not os.path.exists(cond_trait_dir_mqm):
                    os.makedirs(cond_trait_dir_mqm)
else:
    for p in phe_r:
        print p + "\n"
        s = "/"
        seq = (pwd, outfile, p)
        trait_dir = s.join(seq)
        print trait_dir
        if not os.path.exists(trait_dir):
            os.makedirs(trait_dir)
            
        s = "/"
        seq = (pwd, outfile, p, "scanone.out")
        trait_dir_so = s.join(seq)
        print trait_dir_so
        if not os.path.exists(trait_dir_so):
            os.makedirs(trait_dir_so)
        
        s = "/"
        seq = (pwd, outfile, p, "mqm.out")
        trait_dir_mqm = s.join(seq)
        print trait_dir_mqm
        if not os.path.exists(trait_dir_mqm):
            os.makedirs(trait_dir_mqm)    

### Perform scanone, scantwo and stepwise MQM QTL analysis

# Check to see if analysis is being done on server. If so do as below if not see lower else statement
if re.match('y', server):
    # Check to see if you are doing a treatment comparison.
    # If so run each trait individually as well as the difference between treatments and BLUP traits
    if re.match('y', comparison):
        def parallel_comp(x):
            phe_r_col_number = x + 2
            #for p in phe_r:
            (treat, trait) = re.split('\.', str(phe_r[x]))
            print "Working on treatment: " + treat
            print "Working on trait: " + trait
            print "This is phe_r column number: " + str(phe_r_col_number)
            qtlcommand_raw_phenotype = "Rscript do.qtl.comp.R " + str(phe_r_col_number) + " " + treat + " " + trait + " " + outfile + " " +  trait_model + " " + qtl_method
            print qtlcommand_raw_phenotype
            subprocess.call(qtlcommand_raw_phenotype, shell=True)
            #phe_r_col_number += 1
            # Need to add parallelization routine here
        # Call analysis of each trait in parallel    
        #processes_comp = [mp.Process(target=parallel_comp, args=(x,)) for x in range(0,int(nphe_r)-1,1)]
        # Begin looping through processes
        #for p in processes_comp:
            #p.start()
        p = mp.Pool(processes=20)
        p.map_async(parallel_comp, range(0,int(nphe_r)-1,1)).get(9999999)
        #for p in processes_comp:
            #p.join()
        def parallel_diff(x):    
            phe_d_col_number = x + 2
            #for d in phe_d:
            #    print d
            (trait, type, treat1, treat2) = re.split('\.', str(phe_d[x]))
            print trait
            print type
            treat_diff = treat1 + "." + treat2
            print treat_diff
            print "Working on trait: " + phe_d[x] + " name of physical trait is " + trait
            print "This is phe_d column number: " + str(phe_d_col_number)
            qtlcommand_diff_phenotype = "Rscript do.qtl.R " + str(phe_d_col_number) + " " + outfile + " " + trait + " " + " cross.object.diff.Rdata"  + " " +  trait_model + " " + qtl_method
            print qtlcommand_diff_phenotype
            #phe_d_col_number += 1
            subprocess.call(qtlcommand_diff_phenotype, shell=True)
        # Call analysis of each trait in parallel    
        #processes_diff = [mp.Process(target=parallel_diff, args=(x,)) for x in range(0,int(nphe_d)-1,1)]
        # Begin looping through processes
        #for p in processes_diff:
            #p.start()
            
        #for p in processes_diff:
            #p.join()
        p = mp.Pool(processes=20)
        p.map_async(parallel_diff, range(0,int(nphe_d)-1,1)).get(9999999)
        
    else:
        def parallel_raw(x):
            print "Processing individual traits with no comparison..."
            # What phenotype column in R/qtl object are you acting on?
            phe_r_col_number = x + 2
            # Get trait name
            trait = phe_r[x]
            # Print status of what trait is being worked on
            print "Working on trait: " + trait
            print "This is phe_r column number: " + str(phe_r_col_number)
            # Build a text command to invoke in the shell
            qtlcommand_raw_phenotype = "Rscript do.qtl.R " + str(phe_r_col_number) + " " + outfile + " " + trait + " cross.object.raw.Rdata"  + " " +  trait_model + " " + qtl_method
            print qtlcommand_raw_phenotype
            subprocess.call(qtlcommand_raw_phenotype, shell=True)
        # Call analysis of each trait in parallel    
        #processes = [mp.Process(target=parallel_raw, args=(x,)) for x in range(0,int(nphe_r)-1,1)]
        # Begin looping through processes
        #for p in processes:
            #p.start()
            
        #for p in processes:
            #p.join()
        p = mp.Pool(processes=20)
        p.map_async(parallel_raw, range(0,int(nphe_r)-1,1)).get(9999999)
else:
    print "Not on a server system. Will process each trait one-by-one..."
    if re.match('y', comparison):
        phe_r_col_number = 2
        for p in phe_r:
            (treat, trait) = re.split('\.', str(p))
            print "Treatment is: " + treat
            print "Trait is: " + trait
            print "This is phe_r column number: " + str(phe_r_col_number)
            # Build command
            qtlcommand_raw_phenotype = "Rscript do.qtl.comp.R " + str(phe_r_col_number) + " " + treat + " " + trait + " " + outfile + " " +  trait_model + " " + qtl_method
            print qtlcommand_raw_phenotype
            subprocess.call(qtlcommand_raw_phenotype, shell=True)
            phe_r_col_number += 1
        
        phe_d_col_number = 2
        for d in phe_d:
            #    print d
            (trait, type, treat1, treat2) = re.split('\.', str(phe_d[x]))
            print "Trait is: " + trait
            treat_diff = treat1 + "." + treat2
            print "Treatment comparison is: " + treat_diff
            print "Working on trait: " + d + " name of physical trait is " + trait
            print "This is phe_d column number: " + str(phe_d_col_number)
            qtlcommand_diff_phenotype = "Rscript do.qtl.R " + str(phe_d_col_number) + " " + outfile + " " + trait + " " + " cross.object.diff.Rdata"  + " " +  trait_model + " " + qtl_method
            print qtlcommand_diff_phenotype
            subprocess.call(qtlcommand_diff_phenotype, shell=True)
            phe_d_col_number += 1
    else:
        print "Processing individual traits with no comparison..."
        phe_r_col_number = 2
        for p in phe_r:
            print "Working on trait: " + p
            print "This is phe_r column number: " + str(phe_r_col_number)
            qtlcommand_raw_phenotype = "Rscript do.qtl.R " + str(phe_r_col_number) + " " + outfile + " " + p + " " + " cross.object.diff.Rdata"  + " " +  trait_model + " " + qtl_method 
            print qtlcommand_raw_phenotype
            subprocess.call(qtlcommand_raw_phenotype, shell=True)
            phe_r_col_number += 1
            
 
 # If a comparison is specified, make plots of LOD trace for each trait + treatment combination
if re.match('y', comparison):
    # Object 'traits' is defined above 
    for t in traits:
        os.chdir('/' + pwd + '/' + outfile + '/' + t)
        print t
        qtl12comp = glob.glob("[a-z]*")
        files_kept = list()
        for q in qtl12comp:
            if re.match('comparison*', q):
                next
            else:
                files_kept.append(q)
        print(files_kept[0] + "\t" + files_kept[1])
        os.chdir('/' + pwd)
        cond1 = files_kept[0]
        cond2 = files_kept[1]
        compCommand = "Rscript mk.comp.plots.R" + " " + t + " " + cond1 + " " + cond2 + " " + outfile;
        print "Executing " + compCommand
        subprocess.call(compCommand, shell=True)


### Get all genes within the confidence interval for each QTL from scanone and MQM

# Check to see if analysis is being done on server. If so do as below if not see lower else statement

if re.match('y', server):
    
    if re.match('y', comparison):
        def parallel_get_gene_comp(x):
            #(treat, trait) = re.split('\.', str(phe_r[x]))
            (treat, trait) = re.split('\.', phe_r[x])
            print "Getting genes from treatment: " + treat
            print "Getting genes from trait: " + trait
            # Get genes from  output
            # This is directory where QTLs info is stored
            s = "/"
            seq = (pwd, outfile, trait)
            trait_dir = s.join(seq)
            trait_dir_treat = trait_dir + "/" + treat 
            # Build command
            qtlcommand_get_gene_comp = "Rscript get.gene.comp.R " + trait_dir_treat + " " + trait + " " + treat
            print qtlcommand_get_gene_comp
            subprocess.call(qtlcommand_get_gene_comp, shell=True)

        # Call analysis of each trait in parallel    
        #processes_comp = [mp.Process(target=parallel_get_gene_comp, args=(x,)) for x in range(0,int(nphe_r)-1,1)]
        # Begin looping through processes
        #for p in processes_comp:
            #p.start()
            
        #for p in processes_comp:
            #p.join()
        p = mp.Pool(processes=20)
        p.map_async(parallel_get_gene_comp, range(0,int(nphe_r)-1,1)).get(9999999)
        
        
        
        def parallel_get_gene_diff(x):
            (trait, type, treat1, treat2) = re.split('\.', str(phe_d[x]))
            diff_trait = str(phe_d[x])
            print "Getting genes from trait: " + trait
            print "Getting genes from difference trait: " + diff_trait
            # Get genes from  output
            # This is directory where QTLs info is stored
            s = "/"
            seq = (pwd, outfile, trait)
            trait_dir = s.join(seq)
            trait_dir_comp = trait_dir + "/" + 'comparison' 
            # Build command
            qtlcommand_get_gene_diff = "Rscript get.gene.R " + trait_dir_comp + " " + diff_trait
            print qtlcommand_get_gene_diff
            subprocess.call(qtlcommand_get_gene_diff, shell=True)

        # Call analysis of each trait in parallel    
        #processes_comp = [mp.Process(target=parallel_get_gene_diff, args=(x,)) for x in range(0,int(nphe_d)-1,1)]
        # Begin looping through processes
        #for p in processes_comp:
            #p.start()
            
        #for p in processes_comp:
            #p.join()
        p = mp.Pool(processes=20)
        p.map_async(parallel_get_gene_diff, range(0,int(nphe_r)-1,1)).get(9999999)

    else:
        def parallel_get_gene(x):
            # Get trait name
            trait = phe_r[x]
            # Print status of what trait is being worked on
            print "Getting genes from trait: " + trait
            s = "/"
            seq = (pwd, outfile, trait)
            trait_dir = s.join(seq)
            #trait_dir = trait_dir 
            #trait_dir_treat = trait_dir_mqm + trait 
            # Build a text command to invoke in the shell
            qtlcommand_get_gene = "Rscript get.gene.R " + trait_dir + " " + trait
            print qtlcommand_get_gene
            subprocess.call(qtlcommand_get_gene, shell=True)
        # Call analysis of each trait in parallel    
        #processes = [mp.Process(target=parallel_get_gene, args=(x,)) for x in range(0,int(nphe_r)-1,1)]
        # Begin looping through processes
        #for p in processes:
            #p.start()
            
        #for p in processes:
            #p.join()
        p = mp.Pool(processes=20)
        p.map_async(parallel_get_gene, range(0,int(nphe_r)-1,1)).get(9999999)
else:
    
    if re.match('y', comparison):
        for t in traits:
            s = "/"
            seq = (pwd, outfile, t)
            trait_dir = s.join(seq)    
            for c in comp:
                cond_trait_dir = trait_dir + "/" + c 
                print cond_trait_dir
                if os.path.exists(cond_trait_dir):
                    cond_trait_dir_so_qtl = cond_trait_dir + "/" + "scanone.out" + "/" + "summary.table.so." + c + "." + t + ".csv"
                    print cond_trait_dir_so_qtl
                    if os.path.exists(cond_trait_dir_so_qtl):
                        #os.makedirs(cond_trait_dir_so)
                        print "Found it: " + cond_trait_dir_so_qtl
                 
                    cond_trait_dir_mqm_qtl = cond_trait_dir + "/" + "mqm.out" + "/" + "summary.table.mqm." + c + "." + t + ".csv"
                    print cond_trait_dir_mqm_qtl
                    if os.path.exists(cond_trait_dir_mqm_qtl):
                        #os.makedirs(cond_trait_dir_mqm)
                        print "Found it: " + cond_trait_dir_mqm_qtl
    else:
        for p in phe_r:
            print p + "\n"
            s = "/"
            seq = (pwd, outfile, p)
            trait_dir = s.join(seq)
            
            s = "/"
            seq = (pwd, outfile, p, "scanone.out", "summary.table.so." + p + ".csv")
            trait_dir_so_qtl = s.join(seq)
            print trait_dir_so_qtl
            if os.path.exists(trait_dir_so_qtl):
                #os.makedirs(trait_dir_so_qtl)
                print "Found it: " + trait_dir_so_qtl
        
            s = "/"
            seq = (pwd, outfile, p, "mqm.out", "summary.table.mqm." + p + ".csv")
            trait_dir_mqm = s.join(seq)
            print trait_dir_mqm
            if os.path.exists(trait_dir_mqm):
                #os.makedirs(trait_dir_mqm)
                print "Found it: " + trait_dir_mqm_qtl


print "Thank you for using our QTL pipeline =)"