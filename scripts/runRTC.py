import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt

import scipy
from scipy import stats
import sys
import numpy as np
import pandas as pd
from pandas import Series,DataFrame
import sys
import codecs
 
import random
random.seed(0)
from random import sample

import random
random.seed(0)
from random import sample

import os
import load_params as LP
import prepsteps as PS
import allstepsH as AS
import eqtl_funcs as EF
import gtex_funcs as GF
import deprecated_funcs as DF
import allstepslist as ASL


def run_RTC(paramDict, geneName, chrNum, pos, tissue):
    ### RUn this for the real eQTL
    import load_params as LP
    import eqtl_funcs as EF
    import allstepsH as AS
    import numpy as np
    
    #geneName = "ENSG00000172404.4"
    #chrNum = "22"
    #pos = "40407887"
    #tissue = "Breast_Mammary_Tissue"

    print("Start") #SLOW
    slope, intercept,  residuals, indivVector = EF.get_eqtl_stats(chrNum, pos, geneName, tissue, paramDict)  
    print("RTC1")

    
    def calculate_pp(slope_r, intercept, residuals, chrNum, pos, indivVector, paramDict):
        import numpy as np
        permutedResiduals = random.sample(list(residuals), len(residuals))
        genoVector, indivGenoVector = GF.get_dos_vector(chrNum, pos, paramDict["vcfDir"])
        filteredGenoVector, indivVector = GF.filter_and_sort_genotype_vector(genoVector, indivGenoVector, indivVector)
        pp = float(slope_r)*np.array(filteredGenoVector) + permutedResiduals + intercept   
        return pp



   # random.seed(0)

    param_file = "param_files/binkley_local.csv" 

    paramDict = LP.load_param_file(param_file)
    print("RTC2")
    #geneName = "ENSG00000172404.4"
    #chrNum = "22"
    #pos = "40407887"
    #tissue = "Breast_Mammary_Tissue"

    #print(pos)
    slope_r, intercept,  residuals, indivVector = EF.get_eqtl_stats(chrNum, pos, geneName, tissue, paramDict)  
    print("RTC3")

    chrNumLinked = "22"
    ### THIS is the position of the randomly selected
    #posLinked = output_pos_rand_variant(paramDict["gtexDir"] + "linkedeQTLpos.txt", paramDict)
    posLinked = AS.output_pos_rand_variant2(paramDict["gtexDir"] + "linkedeQTLpos.txt", paramDict)
    print("RTC4") #SLOW
    #print("position", posLinked, posLinked[0])
    posLinked = posLinked[0]

    pp = calculate_pp(slope_r, intercept, residuals, chrNumLinked, posLinked, indivVector, paramDict)
    print("RTC5")

    def output_bed(indivVector, expVector, filepath, geneName, paramDict):
        selected = output_gene_info(geneName, paramDict)
        header = output_header2(indivVector)
        with open(filepath, 'w') as outfile:
            outfile.write(header + "\n")
            geneLocation = selected
            outfile.write("\t".join([str(x) for x in geneLocation]) + "\t" + "\t".join([str(x) for x in expVector]) + "\n")
    #print(geneName)        
    output_bed(indivVector, pp, paramDict["tmpSmall"] + "/GTExExpfile.bed", geneName, paramDict)
    print("RTC6") #SLOW

    #This prepares the index files for the pseudophenotype
    print("RTC6A")
    os.system("bgzip /users/michaelbinkley/desktop/RTCstuffs/tmpsmall/GTExExpfile.bed")
    print("RTC6B")
    os.system("tabix -p bed /users/michaelbinkley/desktop/RTCstuffs/tmpsmall/GTExExpfile.bed.gz")
    print("RTC6C")
    #This calculates the RTC score
    os.system("QTLtools rtc --vcf /users/michaelbinkley/desktop/GTEx/chr22_subset_gtex.vcf.gz --bed /users/michaelbinkley/desktop/RTCstuffs/tmpsmall/GTExExpfile.bed.gz --hotspot /users/michaelbinkley/desktop/qtltools_test/hotspots.bed --gwas-cis /users/michaelbinkley/desktop/GTEx/linkedGWAS.txt /users/michaelbinkley/desktop/RTCstuffs/tmpSmall/permutations.txt --normal --out rtc_results2.txt")
    print("RTC7")
#run_RTC(paramDict)

def output_gene_info(geneName, paramDict):
    
    fileIN = open(paramDict["tmpDir"] + '/chrname.txt')
    selected = list()
    for i in fileIN: 
        i = i.rstrip().split('\t')
        
        if i[3]== geneName:
            selected = i
            break
            
    return selected

  
def output_header2(indVector): 

            
            
    header = "\t".join(['#chr',  'start', 'end' , 'gene' , 'length' , 'strand'] + [str(x) for x in indVector] )
    return header

    fileIN.close()

