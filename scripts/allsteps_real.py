import os
import load_params as LP
import prepsteps as PS
import allstepsH as AS
import eqtl_funcs as EF
import gtex_funcs as GF
import deprecated_funcs as DF
import allstepslist as ASL
import runRTC as rrtc
import numpy as np

import load_params as LP
import eqtl_funcs as EF



def prep_real_GWAS(outfilename):
    fileIN = open("bestGWAS.txt")

    fileOUT = open(outfilename, "w")

    for i in fileIN:
        #i = i[0]
        
        i=i.rstrip().split("\t")
        if len(i)<6: continue
        print(i)
        varG = str(i[3]) 

        fileOUT.write(str(varG) + "_b37" + "\t" + "Cancer"  )

    fileOUT.close()
    fileIN.close()
#prep_real_GWAS("bestGWAS_var.txt")


expV2, indivVector2 = EF.get_eqtl_stats_real(chrNum, pos, geneName, tissue, paramDict)  

#print(slope)
#print(intercept)
#print(residuals)
print(expV2)

def output_gene_info_real(geneName, paramDict, chrNum, pos, tissue):
    
    fileIN = open(paramDict["tmpDir"] + '/chrname.txt')
    selected = list()
    for i in fileIN: 
        i = i.rstrip().split('\t')
        
        if i[3]== geneName:
            selected = i
            break
            
    return selected

  
def output_header2_real(indivVector2, geneName, chrNum, pos, tissue): 

            
            
    header = "\t".join(['#chr',  'start', 'end' , 'gene' , 'length' , 'strand'] + [str(x) for x in indivVector2] )
    return header

    fileIN.close()

def output_bed_real(indivVector2, expV2, filepath, geneName, paramDict, chrNum, pos, tissue):
    selected = output_gene_info_real(geneName, paramDict, chrNum, pos, tissue)
    header = output_header2_real(indivVector2, geneName, chrNum, pos, tissue)
    with open(filepath, 'w') as outfile:
        outfile.write(header + "\n")
        geneLocation = selected
        outfile.write("\t".join([str(x) for x in geneLocation]) + "\t" + "\t".join([str(x) for x in expV2]) + "\n")
print(geneName)        
#output_bed_real(indivVector2, expV2, paramDict["tmpSmall"] + "/GTExExpfile_realeQTL.bed", geneName, paramDict)







def run_perm_real(paramDict, geneName, chrNum, pos, tissue):
    numperm = 1
    outlist = list()

    for i in range(numperm):
      
        prep_real_GWAS("bestGWAS_var.txt")
        output_bed_real(indivVector2, expV2, paramDict["tmpSmall"] + "/GTExExpfile_realeQTL.bed", geneName, paramDict, chrNum, pos, tissue)
        prepare_perm_file_real( "/users/michaelbinkley/desktop/RTCstuffs/tmpSmall/permutations_real.txt")
        os.system("bgzip /users/michaelbinkley/desktop/RTCstuffs/tmpsmall/GTExExpfile_realeQTL.bed")
        os.system("tabix -p bed /users/michaelbinkley/desktop/RTCstuffs/tmpsmall/GTExExpfile_realeQTL.bed.gz")



        run_real_eQTL()
        
        fileIN = open('/users/michaelbinkley/desktop/RTCstuffs/rtc_results_realeQTL.txt')
        fileIN.readline()
        genes = list()
        for j in fileIN:
            g = j.rstrip().split(" ")
            name = g[19]

            #genes.append(name)
            outlist.append([geneName, chrNum, pos, tissue, name])
            #first = (j[0]) 
        
        fileIN.close()
    return outlist
    #print((genes))                         
    #fileOUT = open("RTCperm_H1" + str(numperm) + ".txt", "w")
    #for k in outlist:
        #fileOUT.write(k + "\n")       
    #fileOUT.close()


def prepare_perm_file_real(outfilename):  
    fileIN = open('/users/michaelbinkley/desktop/RTCstuffs/tmpsmall/GTExExpfile.bed')
    causalvar = open("besteQTL_var.txt")
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

    
#prepare_perm_file_real( "/users/michaelbinkley/desktop/RTCstuffs/tmpSmall/permutations_real.txt")
  


