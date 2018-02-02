
# coding: utf-8

# In[1]:

import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
get_ipython().magic('matplotlib inline')
import scipy
from scipy import stats
import sys
import numpy as np
import pandas as pd
from pandas import Series,DataFrame
import sys
import codecs
 

from random import sample

import seaborn as sns
sns.set_style('whitegrid')

# guarantee unicode string
_u = lambda t: t.decode('UTF-8', 'replace') if isinstance(t, str) else t
_uu = lambda *tt: tuple(_u(t) for t in tt) 
# guarantee byte string in UTF8 encoding
_u8 = lambda t: t.encode('UTF-8', 'replace') if isinstance(t, unicode) else t
_uu8 = lambda *tt: tuple(_u8(t) for t in tt) 




# In[19]:

def get_rpkm_column_indices_vector(outFileName):
    TPKM = open('/Users/michaelbinkley/Desktop/GTEx/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct')
    fileOUT = open(outFileName, "w")
    TPKM.readline() #Read first line
    TPKM.readline() #Read second line
    #header = TPKM.readline() #Read the third line and store it in header
    #headerSplit = header.rstrip().split("\t")
    for i in TPKM:
        i = i.rstrip().split("\t")
        fileOUT.write("\t".join([str(x) for x in i])+ "\n")
    fileOUT.close()
    TPKM.close()
get_rpkm_column_indices_vector('v7tpkm.txt')


# In[29]:

#extract gene name, chromosome, start position, end position
def get_gene_chr_position(outFileName):
    encode = open('/Users/michaelbinkley/Desktop/RyoProj/gencode.v19.annotation.gtf')
    fileOUT=open(outFileName, "w")
    chromlist3=list()
    encode.readline()
    encode.readline()
    encode.readline()
    encode.readline()
    encode.readline() #five lines to skip at top of file
    for i in encode:
        i = i.rstrip().split("\t")
        chrnum=i[0]

        
        startpos=int(i[3])
        endpos=int(i[4])
        dist=abs(int(i[4]) - int(i[3]))
        sense=i[6]
        gene=i[8]

        chromlist3.append(i[0])
        chromlist3.append(i[3])
        chromlist3.append(i[4])
        chromlist3.append(i[6])
        chromlist3.append(i[8])
        #, float(i[3]), float(i[4]), 
        #fileOUT.write(  str(gene) + "\n")
        fileOUT.write( str(chrnum) + "\t" + str(startpos) + "\t" + str(endpos) +"\t" + str(dist)+ "\t" + str(sense) + "\n")
    print("Number of chr", len(chromlist3))
    fileOUT.close()
    encode.close()
    return chromlist3
#get_gene_chr_position('chrom3.csv')
get_gene_chr_position('chrom2.txt')


# In[18]:

def output_geneID(outFileName):
    fileOUT = open(outFileName, "w")
    IN = open('chrom3.csv')
    chromlist2 = list()
    for i in IN:
        i = i.rstrip().split(";")
        genename = i[0]
        chromlist2.append(i[0])
        fileOUT.write(str(i[0]) + "\n")
    fileOUT.close()
    IN.close()
    return chromlist2
output_geneID('chrom4.txt')
    


# In[30]:

def get_gene_name(outFileName):
    encode = open('chrom4.txt')
    #encode = open('/Users/michaelbinkley/Desktop/RyoProj/gencode.v19.annotation.gtf')
    chromlist3 = open('chrom2.txt')
    fileOUT=open(outFileName, "w")
    genenamelist=list()
    #encode.readline()
    #encode.readline()
    #encode.readline()
    #encode.readline()    

    #encode.readline() #five lines to skip at top of file    
    for i in encode:
        i = i.rstrip().split(";")
        genename = i[0]
        for num in i[0]: #Get rid of "gene_name"
            if num in "gene_id ":
                genename = genename.replace(num, '')
                
            if num in '" "': #get rid of the quotes
                genename = genename.replace(num, '')

        genenamelist.append(genename)

    for  counter, i in enumerate(chromlist3):
        i = i.rstrip().split("\t")
        chrnum=i[0]
        startpos=i[1]
        endpos=i[2]
        dist=i[3]
        sense=i[4]

        #genenamelist.append(chrnum)
        #genenamelist.append(startpos)
        #genenamelist.append(endpos)  
        #genenamelist.append(dist)
        #genenamelist.append(sense)
    
        fileOUT.write(str(chrnum) + "\t" + str(startpos) + "\t" + str(endpos)+  "\t" + genenamelist[counter] + "\t" + str(dist) + "\t" + str(sense)  + "\n")
    #
    #print("Number of genes", len(genename))
    fileOUT.close()
    encode.close()
    chromlist3.close
    return genenamelist

get_gene_name('chrname2.txt')


# In[62]:

def extract_longest_gene():
    fileIN = open("chrname2.txt")    
    longestGeneCoordinates = list()
    currentLine = fileIN.readline().rstrip().split('\t') 
    currentGeneCoordinates = [currentLine[0], int(currentLine[1]), int(currentLine[2]), str(currentLine[3]), int(currentLine[4]), currentLine[5]]#This is a variable that holds the information for the "current" gene.
    #"Current" gene is the gene that we are currently working to figuring out the longest coordinates.
    for counter, i in enumerate(fileIN): #Go through each line
        i = i.rstrip().split("\t")
        i =[i[0], int(i[1]), int(i[2]), str(i[3]), int(i[4]), i[5]]
        
        if i[3] == currentGeneCoordinates[3]: #If this line has the same gene as "current" gene, then continue to compare the coordiantes.
            currentDistance = abs(currentGeneCoordinates[2] - currentGeneCoordinates[1])
            thisLinesDistance = abs(i[2] - i[1])
            if thisLinesDistance > currentDistance: #If this line is longer than our "current" gene, then update the "current" gene.
                currentGeneCoordinates = i
                
        if i[3] != currentGeneCoordinates[3]: 
            #If this line's gene does not match the currentGene's gene, then we must have moved on to the next gene.
            #Since we're now at a new gene, let's store the currentGene
            longestGeneCoordinates.append(currentGeneCoordinates)
            #Also, let's update the currentGene to match this gene. 
            currentGeneCoordinates = i
        if counter %100000==0:
            print(counter)
        # if counter==1000: 
       #     break
   
    #Now, that we've gone thorough all of totalData, le'ts print the info out
    fileOUT = open("chrname2.longest.txt", "w")
    for i in longestGeneCoordinates:
        print(i)
        fileOUT.write("\t".join([str(x) for x in i]) + "\n")  #Join will take a list, and create ONE string. where each element of the list is separated by the character in the quotation marks (here it is a "\t")
        #Also,we're converting the items into a string - because we can't output integers. 
    fileOUT.close() 
    
extract_longest_gene()


# In[20]:

def get_tpkm_vector(outFileName):
    TPKM = open('v7tpkm.txt')
    TPKM.readline()
    fileOUT = open(outFileName, "w")
    geneVector = list()
    #header = TPKM.readline() #Read the third line and store it in header
    #headerSplit = header.rstrip().split("\t")
    for i in TPKM:
        i = i.rstrip().split("\t")
        gene = i[0]
        geneVector.append(i[0])
        fileOUT.write(i[0] + "\n")
    fileOUT.close()
    TPKM.close()
    return geneVector
get_tpkm_vector('v7tpkmV.txt')

def get_tpkm_vector(outFileName):
    TPKM = open('chrname2.longest.txt')
    TPKM.readline()
    fileOUT = open(outFileName, "w")
    geneVector2 = list()
    #header = TPKM.readline() #Read the third line and store it in header
    #headerSplit = header.rstrip().split("\t")
    for i in TPKM:
        i = i.rstrip().split("\t")
        gene = i[3]
        geneVector2.append(i[3])
        fileOUT.write(i[3] + "\n")
    fileOUT.close()
    TPKM.close()
    return geneVector2
get_tpkm_vector('chrname2V.txt')


# In[22]:

fileIN = open('GTExchrname2.txt')

def file_len(fileIN):
    with open(fileIN) as f:
        for i, l in enumerate(f):
            pass
    return i + 1
file_len("GTExchrname2.txt")


# In[21]:

def output_expression_for_TS(outfilename):
    fileIN = open("chrname2.longest.txt")
    gv = get_tpkm_vector('v7tpkmV.txt')


    fileOUT = open(outfilename, "w")
    TopExp = list()
    SortedGenes = list()
    for i in fileIN:
        i = i.rstrip().split("\t")
        gene = i[3]
        if gene in gv:
            fileOUT.write("\t".join([str(x) for x in i]) + "\n")            

        else:
            continue

                    
    fileIN.close()
    fileOUT.close()
print('done')
output_expression_for_TS("GTExchrname2.txt")



# In[ ]:




# In[6]:

def mean_pos_top_genes(outfilename):
    tgp = open("GTExchrname2.txt")
    fileOUT = open(outfilename, "w")
    meanpos=list()
    for i in tgp:
        i = i.rstrip().split("\t")

        chrnum=i[0]
        startpos=float(i[1])
        endpos=float(i[2])
        gene=str(i[3])
        dist=float(i[4])
        sense=i[5]
        
        
        meanpos.append([chrnum, startpos, endpos, gene, dist, sense])
        fileOUT.write("\t".join([str(x) for x in i]) + "\n")

    fileOUT.close()
    tgp.close()
    return meanpos
#mean_pos_top_genes("mpos.txt")


def GTEx_expression_matrix(outfilename):
    print("started")
    gtexFile = open('v7tpkm.txt')
    gtexFile.readline()
    fileOUT = open(outfilename, "w")
   # mposData = mean_pos_top_genes("mpos.txt")
    with open("mpos.txt") as fileIN:
        mposData = [x.rstrip().split("\t") for x in fileIN.readlines()]
    print("Finished")
    counter = 0
    for i in gtexFile:
        counter +=1
        
        i=i.rstrip().split("\t")

        genex = i[0]
        tpkm = i[2:]
        
        if counter%1000==0: print(counter)
        for chrnum, startpos, endpos, gene, dist, sense in mposData:
        #for mposGene, mposChr, mposPos in mposData:
    
            if gene != genex:

                continue
                
            else:

                fileOUT.write("\t".join([chrnum, startpos, endpos, gene, dist, sense]) + '\t' + "\t".join(tpkm) + "\n")
    fileOUT.close()
    gtexFile.close()   
#GTEx_expression_matrix("gtexQTLexp.txt")


# In[5]:

def output_header(outfilename):

    fileIN = open('v7tpkm.txt')
    fileOUT = open(outfilename, "w")
#fileIN.readline()
#fileIN.readline()
    fileOUT.write(fileIN.readline() + "\n")
    print(fileIN.readline())
#output_header("GTEXheader.txt")


# In[60]:

def sort_combined_gwas(outfilename):
    fileIN = open("gtexQTLexp.txt")
    sortedData = list()
    fileOUT = open(outfilename, "w")
    for i in fileIN:
        i=i.rstrip().split("\t")
        chrnum=str(i[0])
        startpos=float(i[1])
        endpos=i[2]
        gene=i[3]
        dist=i[4]
        sense=i[5]
        tpkm=i[6:]
        sortedData.append( [str(i[0]), float(i[1]), float(i[2]), i[3], float(i[4]), i[5]])
  

    
    sortedScoreData = sorted(sortedData, key = lambda x : (x[0], x[1]), reverse=False) #Sort by the element in the column index 1
    
    #print (sortedScoreData)


    for i in sortedScoreData:
        #print(i)
        fileOUT.write("\t".join([str(x) for x in i]) + '\t' + "\t".join(tpkm) + "\n")
    #return [x[0] for x in sortedScoreData]
    fileOUT.close()
    fileIN.close()
sort_combined_gwas('sortedgtexQTLexp.txt')


# In[64]:

def add_header():
    filenames = ['']
filenames = ['GTEXheader2.txt', 'sortedgtexQTLexp.txt']
with open('GTExExpfile.bed', 'w') as outfile:
    for fname in filenames:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)
add_header()


# In[38]:

fileIN = open('/Users/michaelbinkley/Desktop/GTEx_Analysis_v7_eQTL/Breast_Mammary_Tissue.v7.egenes.txt')


#fileIN.readline()

print(fileIN.readline())


# In[86]:

import fileinput
with open('try2.txt', 'w') as fout, fileinput.input('GTEXheader2.txt', 'sortedgtexQTLexp.txt') as fin:
    for line in fin:
        fout.write(line)


# In[3]:

Y = "Thyroid"
def get_tissue_list():
    #Output the list of all tissues.
    Ann = open('/Users/michaelbinkley/Desktop/GTEx/GTEx_v7_Annotations_SampleAttributesDS.txt')
    Ann.readline() #Reads first line
    sampleIdList = list() #create a new list with sample name and tissue
    tissueList = list()
    for i in Ann: #for each row
        i = i.rstrip().split('\t')
        if i[6].rstrip()=="": continue  #Remove any tissues that have a blank name
        tissueList.append(i[6])
    Ann.close()

    uniqueTissueList = sorted(list(set(tissueList)))
    
    #print("Number of unique tissues:", len(uniqueTissueList))
    return uniqueTissueList
#get_tissue_list()

def get_sample_names(Y):
    Ann = open('/Users/michaelbinkley/Desktop/GTEx/GTEx_v7_Annotations_SampleAttributesDS.txt')
    Ann.readline() #Reads first line
    sampleIdList = list() #create a new list with sample name and tissue
    for i in Ann: #for each row
        i = i.rstrip().split('\t')
        if i[6]==Y :
            sampleId = i[0]
            sampleIdList.append(sampleId)
        #print(len(sampleIdList))
    Ann.close()
    print("Number of ", Y, "samples found:", len(sampleIdList))
    return sampleIdList


####Identifying the specific tissue in the expression data
def get_rpkm_column_indices_vector(Y):
    RPKM = open('/Users/michaelbinkley/Desktop/GTEx/GTExExpfile.bed')

    header = RPKM.readline() #Read the third line and store it in header
    headerSplit = header.rstrip().split("\t")
    
    sampleNames = get_sample_names(Y)
    rpkmColumnIndexVector = list()
    for i in range(6, len(headerSplit)): #Start at 2 to skip the first two columns (which is the name of the gene and the "description" of the gene)
        if headerSplit[i] in sampleNames:
            rpkmColumnIndexVector.append(1)
        else:
            rpkmColumnIndexVector.append(0)
    RPKM.close()
    return rpkmColumnIndexVector

### This one - will output the actual list of indices that have that vector (instead of 0's and 1's)
def get_rpkm_column_indices_list(Y):
    RPKM = open('/Users/michaelbinkley/Desktop/GTEx/GTExExpfile.bed')

    header = RPKM.readline() #Read the first line and store it in header
    headerSplit = header.rstrip().split("\t")
    
    sampleNames = get_sample_names(Y)
    rpkmColumnIndexList = list()
    for i in range(6, len(headerSplit)): #Start at 6 to skip the first two columns (which is the name of the gene and the "description" of the gene)
        if headerSplit[i] in sampleNames:
            rpkmColumnIndexList.append(i)
    #print("Here's the list", i)
    print(rpkmColumnIndexList)
    return rpkmColumnIndexList


    RPKM.close()
#get_sample_names(Y)
#get_rpkm_column_indices_vector(Y)
#get_rpkm_column_indices_list(Y)


# In[4]:


def create_average_expression_matrix():
    
    #tissueList = get_tissue_list()
    tissue = "Thyroid"
    RPKM = open('/Users/michaelbinkley/Desktop/GTEx/GTExExpfile.bed')

    header = RPKM.readline() #Read the first line and store it in header
    RPKM.readline()
    #rpkmColumnIndexVector
    #Change added here to speed up the code
    #Added a dictionary ---- where tissueDict["tissueName"] will output the list of indices for the columns of the rpkm file

    indexList = get_rpkm_column_indices_list(tissue)
    sampleList = get_sample_names(tissue)
   
    #Filter the tissue list to only include tissues with at least one index matching to it. 
    #filteredTissueList = [Y  for Y in tissueList if len(tissueDict[Y])>0]
    #Y = "Thyroid"    
    fileOUT = open(tissue + "averageExpressionv7.txt", "w")
    #Add a header line
    fileOUT.write("\t".join(["#chr"] + ["start"] + ["end"] + ["gene"] + ["length"] + ["strand"] + sampleList) + "\n")
    
     # Go through each gene.
    counter = 0
    for rpkmLine in RPKM:
        rpkmLine = rpkmLine.rstrip().split("\t")

        chrnum=rpkmLine[0]
        startpos=rpkmLine[1]
        endpos=rpkmLine[2]
        gene=rpkmLine[3]
        dist=rpkmLine[4]
        sense=rpkmLine[5]
        
        geneList = list()
        geneList.append(rpkmLine[0])
        geneList.append(rpkmLine[1])
        geneList.append(rpkmLine[2])
        geneList.append(rpkmLine[3])
        geneList.append(rpkmLine[4])
        geneList.append(rpkmLine[5])
        #Go through each tissue.
        
            #print("tissue name: ", tissue)
        
        tissueSpecificExpression = [str(rpkmLine[idx])   for idx in indexList]
        geneList = geneList + tissueSpecificExpression
            #average = np.mean(tissueSpecificExpression)
            #geneList.append(str(indexList)) 
        fileOUT.write("\t".join(geneList) + "\n")
        counter = counter +1
        
        #if counter > 100: break  #We're not using this line now, but this line is good for testing. "break" means leave the for loop.
        if counter%1000==0: #This section is also for testing. "%" means remainder (in the basic arithmetic sense). 0%1000==0 , 1%1000=1
            print(counter)
            break
    fileOUT.close()
    RPKM.close()
#create_average_expression_matrix()

    
    


# In[2]:

def ID_hotspots(outfilename): 
    fileIN = open('/Users/michaelbinkley/Desktop/hotspots_b37_hg19.bed')
    hotspots = list()
    fileOUT = open(outfilename, "w")
    for i in fileIN: 
        i = i.rstrip().split('\t')

        
        Chr = i[0]
        for num in i[0]: #Get rid of "chr in the file"
            if num in "chr":
                Chr = Chr.replace(num, '')       
        start = int(i[1])
        end = int(i[2])
        middle = ((start+end)/2)
        hotspots.append(Chr)
        hotspots.append(start)
        hotspots.append(end)
        hotspots.append(middle)
        fileOUT.write(Chr + "\t" + str(start) + "\t" + str(end) + "\t" + str(middle)  + "\n")
    fileIN.close()
    fileOUT.close()
    return hotspots
ID_hotspots('hotspots.txt')


# In[3]:

def output_coldspot(outfilename):

    egene = open('/Users/michaelbinkley/Desktop/GTEx_Analysis_v7_eQTL/Breast_Mammary_Tissue.v7.egenes.txt')
    egene.readline()

    fileOUT = open(outfilename, "w")
    hotspots = open('hotspots.txt')
    for j in hotspots: 
        j = j.rstrip().split('\t')
        Chr = str(j[0])
        start = float(j[1])
        end = float(j[2])
        
    
    for i in egene:
        
        i = i.rstrip().split('\t')
        chrnum = str(i[2])
        gene = i[0]
        position = float(i[14])
        if chrnum == Chr and position > start and position <end :
            continue
        else: 
            fileOUT.write("\t".join([str(x) for x in i])+ "\n")
            
    
    egene.close()
    fileOUT.close()
    hotspots.close()
output_coldspot('coldspotlist.txt')


# In[39]:

#Next steps: 
# Need to make a list of GWAS, eQTLs within a cold spot, because you need to randomly pick from them, 
#BUT, that's a lot of adata, and a nonefficient way to ask

# MAYBE, you first need to start with the GWAS variant, then identify the nearest hot spot, you have the hotspot file
# then, identify eQTLs that are nearby (nearest neighbor)





# In[12]:

def ID_hotspots(outfilename): 
    fileIN = open('coldspots.txt')
    hotspots = list()
    fileOUT = open(outfilename, "w")
    for i in fileIN: 
        i = i.rstrip().split('\t')

        
        Chr = i[0]
        for num in i[0]: #Get rid of "chr in the file"
            if num in "chr":
                Chr = Chr.replace(num, '')       
        start = int(i[1])
        end = int(i[2])
        middle = str(i[3])
        bound = str (i[4])
        hotspots.append([Chr, start, end, middle, bound])
        #hotspots.append(start)
        #hotspots.append(end)
        #hotspots.append(middle)
        fileOUT.write(Chr + "\t" + str(start) + "\t" + str(end) + "\t" + str(middle)  +  "\t" + str(bound) + "\n")
    fileIN.close()
    fileOUT.close()
    return hotspots
#ID_hotspots('coldspots2.txt')

def find_variant_coldspot(outfilename):
    fileIN = open('/Users/michaelbinkley/Desktop/exSNP3.txt')
    
    fileOUT = open(outfilename, "w")
    hspots = ID_hotspots('coldspots2.txt')



    for i in fileIN:
        
        i = i.rstrip().split('\t')
        chrnum = str(i[0])
        snp = (i[2])
        position = float(i[1])
        
        for  Chr, start, (end), middle, (bound) in hspots:
            if Chr  != chrnum :
                continue
            else: 
                if position < float(bound) and position > float(end): 
                    #distance = abs(position - middle)
                    fileOUT.write("\t".join([ str(Chr), str(start), str(end), str(middle), str(bound) , snp, str(position)] ) + "\n")
                else: 
                    continue
    fileIN.close()
    fileOUT.close()

find_variant_coldspot('excoldspot2.txt')


# In[105]:

def sort_gwas(outfilename):
    fileIN = open('excoldspot.txt')
    sortedData = list()
    fileOUT = open(outfilename, "w")
    for i in fileIN:
        i=i.rstrip().split("\t")
        if i[2]=="X": continue
        sortedData.append( [i[0], (i[1]), (i[2]), (i[3]), (i[4]), (i[5]), float(i[6])] )
  

    
    sortedScoreData = sorted(sortedData, key = lambda x : (x[6]), reverse=False) #Sort by the element in the column index 1
    
    #print (sortedScoreData)


        for i in sortedScoreData:
        #print(i)
            fileOUT.write("\t".join([str(x) for x in i]) + "\n")
        return [x[0] for x in sortedScoreData]
    fileOUT.close()
    fileIN.close()
#sort_gwas('sortedSNPex.txt')


# In[104]:

#### need a function to then identify the shortest distance
def extract_longest_gene():
    fileIN = open("sortedSNPex.txt")    
    longestGeneCoordinates = list()
    currentLine = fileIN.readline().rstrip().split('\t') 
    currentGeneCoordinates = [float(currentLine[0]),currentLine[1], currentLine[2], float(currentLine[3]), str(currentLine[4]), float(currentLine[5]), float(currentLine[6])]#This is a variable that holds the information for the "current" gene.
    #"Current" gene is the gene that we are currently working to figuring out the longest coordinates.
    for counter, i in enumerate(fileIN): #Go through each line
        i = i.rstrip().split("\t")
        i =[float(i[0]), i[1], i[2], float(i[3]), str(i[4]), float(i[5]), float(i[6])]

          
        if i[4] == currentGeneCoordinates[4]: #If this line has the same gene as "current" SNP, then continue to compare the coordiantes.
            currentDistance = abs(currentGeneCoordinates[6])
            thisLinesDistance = abs(i[6])
            if thisLinesDistance < currentDistance: #If this line is shorter than our "current" gene, then update the "current" gene.
                currentGeneCoordinates = i
                
        if i[4] != currentGeneCoordinates[4]: 
            #If this line's gene does not match the currentGene's gene, then we must have moved on to the next gene.
            #Since we're now at a new gene, let's store the currentGene
            longestGeneCoordinates.append(currentGeneCoordinates )
            #Also, let's update the currentGene to match this gene. 
            currentGeneCoordinates = i
        if counter %100000==0:
            print(counter)

       
    fileOUT = open("SNPhotspot.txt", "w")
    for i in longestGeneCoordinates:
        print(i)
        fileOUT.write("\t".join([str(x) for x in i]) + "\n")  #Join will take a list, and create ONE string. where each element of the list is separated by the character in the quotation marks (here it is a "\t")
        #Also,we're converting the items into a string - because we can't output integers. 
    fileOUT.close() 
    
extract_longest_gene()


# In[13]:

def create_list():
    fileIN = open("excoldspot2.txt")
    gwas = list()
    for i in fileIN:
        i = i.rstrip().split('\t')
        Chr = str(i[0] )
        start = i[1] 
        end = float(i[2]) 
        middle = float(i[3] ) 
        bound = float(i[4])
        snp = i[5] 
        position = float(i[6]) 

        gwas.append([Chr, start, end, middle, bound, snp, position])    
    return gwas
create_list()

def find_variant_coldspot(outfilename):
    fileIN = open('/Users/michaelbinkley/Desktop/GTEx_Analysis_v7_eQTL/Breast_Mammary_Tissue.v7.egenes.txt')
    fileIN.readline()
    fileOUT = open(outfilename, "w")
    gwas = create_list()



    for i in fileIN:
        
        i = i.rstrip().split('\t')
        chrnum = str(i[13])
        eqtl = (i[18])
        posQ = float(i[14])
        
        for  Chr, start, end, middle, bound, snp, position in gwas:
            if Chr  != chrnum :
                continue
            else: 
                
                
                if posQ < bound and posQ > end:
                    fileOUT.write("\t".join([ str(Chr), str(start), str(end), str(middle), str(bound), snp, str(position), eqtl, str(posQ)] ) + "\n")
                else: 
                    continue
    fileIN.close()
    fileOUT.close()

find_variant_coldspot('eQTLcoldspot.txt')


# In[ ]:

### calc pearson r2, for a particular SNP, each individual will have 0, 1, or 2. For this, we want to find
## a linked SNP w/in the cold spot
### Use Chr 22 file
### Possibly use Plink tool. , plink 1.9, choose best eQTL (based on p value for each gene)

