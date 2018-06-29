

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

### RUn this for the real eQTL, this does not use residuals or calculate residuals.
import load_params as LP
import eqtl_funcs as EF

#geneName = "ENSG00000172404.4"
#chrNum = "22"
#pos = "40407887"
#tissue = "Breast_Mammary_Tissue"
#22_40407887_G_A_b37

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

  
def output_header2_real(indivVector2): 

            
            
    header = "\t".join(['#chr',  'start', 'end' , 'gene' , 'length' , 'strand'] + [str(x) for x in indivVector2] )
    return header

    fileIN.close()

def output_bed_real(indivVector2, expV2, filepath, geneName, paramDict, chrNum, pos, tissue):
    selected = output_gene_info_real(geneName, paramDict, chrNum, pos, tissue)
    header = output_header2_real(indivVector2)
    with open(filepath, 'w') as outfile:
        outfile.write(header + "\n")
        geneLocation = selected
        outfile.write("\t".join([str(x) for x in geneLocation]) + "\t" + "\t".join([str(x) for x in expV2]) + "\n")
print(geneName)        
#output_bed_real(indivVector2, expV2, paramDict["tmpSmall"] + "/GTExExpfile_realeQTL.bed", geneName, paramDict)



