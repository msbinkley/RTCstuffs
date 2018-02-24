#!/usr/bin/env python3

import gzip, os, subprocess 


def get_exp_vector(geneName, tissue, gtexDir):
    '''
    Output a vector of expression levels for a given gene and tissue 

    Note:  Warning, the order of this expression vector may not be the same as the genotype vector.
    '''

    filePath = gtexDir + "/" + tissue + ".v7.normalized_expression.bed.gz"
    fileIN = gzip.open(filePath)
    counter = 0
    matchingLine = ""
    header = fileIN.readline().decode().rstrip().split('\t')
    individualNames = header[4:]
    for i in fileIN:
        i =i.decode().rstrip().split("\t")
        if i[3]==geneName:
            matchingLine = i
            break
    fileIN.close()
    expVector = matchingLine[4:]
    expVector = [float(x) for x in expVector]
    return expVector, individualNames


def get_dos_vector(chrNum, pos, vcfDir):
    '''
    Outputs the dosage vector for a particular chromoosome and position. 
    
    #Make sure file has a .tbi. 
    If it doesn't, then do the following: 
    DN527o9v:vcfDir ryosukekita$ gunzip chr22_subset_gtex.vcf.gz
    DN527o9v:vcfDir ryosukekita$ bgzip chr22_subset_gtex.vcf
    DN527o9v:vcfDir ryosukekita$ tabix -p vcf chr22_subset_gtex.vcf.gz
    DN527o9v:vcfDir ryosukekita$ tabix chr22_subset_gtex.vcf.gz

    TEST SNP:  40051275        22_40051275_G_GC_b37  
    '''
    print("\tGetting dosage vector for :", chrNum, pos, vcfDir)
    vcfFilePath= vcfDir + "/chr" + chrNum + "_subset_gtex.vcf.gz"
    command = ["tabix", vcfFilePath ,  chrNum + ":" + pos + "-" + pos]
    proc = subprocess.Popen(command, stdout = subprocess.PIPE)
    output, err = proc.communicate()
    if len(output)==0:
        print("Error: There is a problem with the tabix output. Check to make sure the vcf file is tabix-indexed. May need to reindex. Exiting. ")
        sys.exit()
    genotypes = output.decode().split("\t")[9:]
    genotypes = [x.split(":")[0] for x in genotypes]
    dosageSplit = [[int(y) for y in x.split("/")]    for x in genotypes]
    dosage =  [sum(x)    for x in dosageSplit]
   

    print("\tGetting header")
    command = ["tabix", vcfFilePath ,  "-H", chrNum + ":"  + pos + "-" + pos]
    proc = subprocess.Popen(command, stdout = subprocess.PIPE)
    output, err = proc.communicate()
    indivs = output.decode().split("\n")[-2].split("\t")[9:]
    return dosage, indivs


def filter_and_sort_genotype_vector(genoV, indivGenoV, indivExpV):
    '''
    Filters and sorts the genotype Vector so that the indiivduals matching the indivExpVector

    '''
    indivGenoVIdx = [indivGenoV.index(x) for x in indivExpV]
    filteredAndSortedGenoV = [genoV[x]  for x in indivGenoVIdx]
    return filteredAndSortedGenoV, indivExpV





