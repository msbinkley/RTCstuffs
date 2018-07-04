import deprecated_funcs as DF
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

def prep_real_eQTL(paramDict, outfilename):
    DF.extract_best_eqtl_real(paramDict)
    fileIN = open("besteQTLs.txt")
    fileOUT = open(outfilename, "w")
    for i in fileIN:
        i=i.rstrip().split("\t")
        varG = str(i[3])
        fileOUT.write(str(varG) + "\n")

    fileOUT.close()
    fileIN.close()
    
    
    
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
    
### RUn this for the real eQTL, this does not use residuals or calculate residuals.
import load_params as LP
import eqtl_funcs as EF


def get_bed_real(chrNum, pos, geneName, tissue, paramDict):
    import deprecated_funcs as DF
    import eqtl_funcs as EF  
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
        selected = output_gene_info_real(geneName, paramDict)
        header = output_header2_real(indivVector2)
        with open(filepath, 'w') as outfile:
            outfile.write(header + "\n")
            geneLocation = selected
            outfile.write("\t".join([str(x) for x in geneLocation]) + "\t" + "\t".join([str(x) for x in expV2]) + "\n")
    print(geneName)   
    
    #output_bed_real(indivVector2, expV2, paramDict["tmpSmall"] + "/GTExExpfile_realeQTL.bed", geneName, paramDict)

#get_bed_real(chrNum, pos, geneName, tissue, paramDict):
    


