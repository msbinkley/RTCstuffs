import os
import numpy as np
import eqtl_funcs as EF
import gtex_funcs as GF
import deprecated_funcs as DF
import load_params as LP
#param_file = "param_files/ryo_local.csv"
param_file = "param_files/binkley_local.csv" 
paramDict = LP.load_param_file(param_file)


import random
random.seed(0)
from random import sample

geneName = "ENSG00000172404.4"
    
######This step randomly picks an eQTL in the coldspot and the position

def select_rand_variant(outfilename, paramDict): 
    fileIN = open('filteredgenotype.txt')
    fileOUT = open(outfilename, "w")
    rnd = list()
    for i in fileIN:

        i = i.rstrip().split('\t')
        rnd.append(i[2])

   # seed = random.randrange(sys.maxsize)
   # rng = random.Random(seed)
    rand = sample(rnd[0:],1)


 #   print("Seed was:", seed)
   # print(rnd[0:])
    print("eqtl", rand)
    for i in rand:
        fileOUT.write("\t".join([str(x) for x in rand] ) +  "\n")

    fileOUT.close()
    fileIN.close()
    return rand
#select_rand_variant(paramDict["gtexDir"] + "eQTLcausal.txt", paramDict)

def output_pos_rand_variant(outfilename, paramDict): 
    fileIN = open('filteredgenotype.txt')
    fileOUT = open(outfilename, "w")
    sand = open(paramDict["gtexDir"] +  'eQTLcausal.txt')
    for i in sand: 
        i = i.rstrip().split('\t') 
        rand = str(i[0])

    for j in fileIN: 
        j = j.rstrip().split('\t')
        SNP = str(j[2])
        random = str(rand)

        if SNP == random: 

            fileOUT.write(str(j[1]) + "\n")
            break
    fileOUT.close()
    fileIN.close()
#output_pos_rand_variant(paramDict["gtexDir"] + "eQTLcausalpos.txt", paramDict)    


# In[10]:

#####This step randomly picks a GWAScausal in the coldspot

def select_rand_variant_gwas(outfilename, paramDict): 
    fileIN = open('filteredgenotype.txt')
    fileOUT = open(outfilename, "w")
    rnd = list()
    for i in fileIN:

        i = i.rstrip().split('\t')
        rnd.append(i[2])
    rand = sample(rnd[0:],1)
    print("gwas", rand)
    for i in rand:
        fileOUT.write("\t".join([str(x) for x in rand] ) +  "\n")
    fileOUT.close()
    fileIN.close()
#select_rand_variant(paramDict["gtexDir"] + "GWAScausal.txt")

def output_pos_rand_variant_gwas(outfilename, paramDict): 
    fileIN = open('filteredgenotype.txt')
    fileOUT = open(outfilename, "w")
    sand = open(paramDict["gtexDir"] +  'GWAScausal.txt')
    for i in sand: 
        i = i.rstrip().split('\t') 
        rand = str(i[0])

    for j in fileIN: 
        j = j.rstrip().split('\t')
        SNP = str(j[2])
        random = str(rand)

        if SNP == random: 

            fileOUT.write(str(j[1]) + "\n")
            break
    fileOUT.close()
    fileIN.close()
#output_pos_rand_variant(paramDict["gtexDir"] + "GWAScausalpos.txt", paramDict)    



# In[12]:

### This will be the GWAS to run the RTC score 
###This will be the eQTL to run the RTC score 
#### THIS IS NOW FIXED, it will output randomly linked variant that is within the cold spot

def filter_genov_output_gwas(outfilename, paramDict):
    fileIN = open('filteredgenotype.txt')
    #fileOUT = open(outfilename, "w")
    pos = list()
    for j in fileIN: 
        j = j.rstrip().split('\t')
        SNP = str(j[2])  
        pos.append(SNP)
        #fileOUT.write(str(SNP) + "\n")

    #fileOUT.close()
    fileIN.close()
    return pos
#filter_genov_output_gwas(paramDict["gtexDir"] + "filteredSNPs.txt", paramDict)

        
def filter_plink_output_gwas(outfilename, paramDict):
    fileIN = open('/users/michaelbinkley/desktop/RTCstuffs/plink.tags')
    fileOUT = open(outfilename, "w") 
    pos = filter_genov_output_gwas(paramDict["gtexDir"] + "filteredSNPs.txt", paramDict)
    for i in fileIN: 
        i = i.rstrip().split('\t') 
        var = str(i[0])

            
        if var in pos: 
            fileOUT.write(str(var) + "\n")

        else: 
            continue

    fileOUT.close()
    fileIN.close()
filter_plink_output_gwas(paramDict["gtexDir"] + "filteredplink.txt", paramDict) 
    
def select_rand_variant_gwas2(outfilename, paramDict): 
    fileIN = open(paramDict["gtexDir"] + "filteredplink.txt")
    fileOUT = open(outfilename, "w")
    rnd = list()
    for i in fileIN:

        i = i.rstrip().split('\t')
        rnd.append(i[0])

        
    rand = sample(rnd[0:],1)

    for i in rand:
        fileOUT.write("\t".join([str(x) for x in rand] ) +"\t" + 'Cancer' + "\n")
    fileOUT.close()
    fileIN.close()
#select_rand_variant_gwas2(paramDict["gtexDir"] + "linkedGWAS.txt", paramDict)

def output_pos_rand_variant_gwas2(outfilename, paramDict): 
    fileIN = open('filteredgenotype.txt')
    fileOUT = open(outfilename, "w")
    sand = open(paramDict["gtexDir"] +  'linkedGWAS.txt')
    posI = list()
    for i in sand: 
        i = i.rstrip().split('\t') 
        rand = str(i[0])

    for j in fileIN: 
        j = j.rstrip().split('\t')
        SNP = str(j[2])
        random = str(rand)

        if SNP == random: 

            fileOUT.write(str(j[1]) + "\n")
            posI.append(j[1])
            break
    return posI
    fileOUT.close()
    fileIN.close()
#output_pos_rand_variant2(paramDict["gtexDir"] + "linkedGWASpos.txt", paramDict)  


# In[13]:

###This will be the eQTL to run the RTC score 
#### THIS IS NOW FIXED, it will output randomly linked variant that is within the cold spot

def filter_genov_output(outfilename, paramDict):
    fileIN = open('filteredgenotype.txt')
    #fileOUT = open(outfilename, "w")
    pos = list()
    for j in fileIN: 
        j = j.rstrip().split('\t')
        SNP = str(j[2])  
        pos.append(SNP)
        #fileOUT.write(str(SNP) + "\n")

    #fileOUT.close()
    fileIN.close()
    return pos
#filter_genov_output(paramDict["gtexDir"] + "filteredSNPs.txt", paramDict)

        
def filter_plink_output(outfilename, paramDict):
    fileIN = open('/users/michaelbinkley/desktop/RTCstuffs/plink.tags')
    fileOUT = open(outfilename, "w") 
    pos = filter_genov_output(paramDict["gtexDir"] + "filteredSNPs.txt", paramDict)
    for i in fileIN: 
        i = i.rstrip().split('\t') 
        var = str(i[0])

            
        if var in pos: 
            fileOUT.write(str(var) + "\n")

        else: 
            continue

    fileOUT.close()
    fileIN.close()
#filter_plink_output(paramDict["gtexDir"] + "filteredplink.txt", paramDict) 
    
def select_rand_variant2(outfilename, paramDict): 
    fileIN = open(paramDict["gtexDir"] + "filteredplink.txt")
    fileOUT = open(outfilename, "w")
    rnd = list()
    for i in fileIN:

        i = i.rstrip().split('\t')
        rnd.append(i[0])

        
    rand = sample(rnd[0:],1)

    for i in rand:
        fileOUT.write("\t".join([str(x) for x in rand] ) +"\t" + 'Cancer' + "\n")
    fileOUT.close()
    fileIN.close()
#select_rand_variant2(paramDict["gtexDir"] + "linkedeQTL.txt", paramDict)

def output_pos_rand_variant2(outfilename, paramDict): 
    fileIN = open('filteredgenotype.txt')
    fileOUT = open(outfilename, "w")
    sand = open(paramDict["gtexDir"] +  'linkedeQTL.txt')
    posI = list()
    for i in sand: 
        i = i.rstrip().split('\t') 
        rand = str(i[0])

    for j in fileIN: 
        j = j.rstrip().split('\t')
        SNP = str(j[2])
        random = str(rand)

        if SNP == random: 

            fileOUT.write(str(j[1]) + "\n")
            posI.append(j[1])
            break
    return posI
    fileOUT.close()
    fileIN.close()
#output_pos_rand_variant2(paramDict["gtexDir"] + "linkedeQTLpos.txt", paramDict)   



def prepare_perm_file(outfilename): 
    fileIN = open('/users/michaelbinkley/desktop/RTCstuffs/tmpsmall/GTExExpfile.bed')
    causalvar = open('/users/michaelbinkley/desktop/GTEx/linkedeQTL.txt')
    fileIN.readline()
    fileOUT = open(outfilename, "w")
    for j in causalvar: 
        j = j.rstrip().split('\t')
        cvar = j[0]
    for i in fileIN: 
        i = i.rstrip().split('\t')
        gene = i[3]
        chrn = i[0]
        start = i[1]
        end = i[2]
        dist = i[4]
        sense = i[5]
        
        fileOUT.write("\t".join([ str(gene), str(chrn),  str(start) , str(end), sense, dist, dist, str(cvar), chrn , 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA' ] ) + "\n")

    fileIN.close()
    fileOUT.close()

