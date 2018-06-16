
import os
import load_params as LP
import prepsteps as PS
import allstepsH as AS
import eqtl_funcs as EF
import gtex_funcs as GF
import deprecated_funcs as DF
#param_file = "param_files/ryo_local.csv"
param_file = "param_files/binkley_local.csv" 

paramDict = LP.load_param_file(param_file)



#This step is incredibly fast. No need to optimize.
def run_all_steps_H0(paramDict):  
    #This selects a 'random' eQTL and GWAS within the selected coldspot
    AS.select_rand_variant(paramDict["gtexDir"] + "eQTLcausal.txt", paramDict)
    AS.output_pos_rand_variant(paramDict["gtexDir"] + "eQTLcausalpos.txt", paramDict)  
    AS.select_rand_variant_gwas(paramDict["gtexDir"] + "GWAScausal.txt", paramDict)
    AS.output_pos_rand_variant_gwas(paramDict["gtexDir"] + "GWAScausalpos.txt", paramDict)
    
    #This identifies the linked eQTL
    os.system("plink --vcf /users/michaelbinkley/desktop/GTEx/chr22_subset_gtex_copy2.vcf --show-tags /users/michaelbinkley/desktop/GTEx/eQTLcausal.txt --tag-r2 0.5")
    AS.filter_genov_output(paramDict["gtexDir"] + "filteredSNPs.txt", paramDict)
    AS.filter_plink_output(paramDict["gtexDir"] + "filteredplink.txt", paramDict) 
    AS.select_rand_variant2(paramDict["gtexDir"] + "linkedeQTL.txt", paramDict)
    AS.output_pos_rand_variant2(paramDict["gtexDir"] + "linkedeQTLpos.txt", paramDict)   
    AS.prepare_perm_file( "/users/michaelbinkley/desktop/RTCstuffs/tmpSmall/permutations.txt")
    
    #This now identifies the linked GWAS variant
    os.system("plink --vcf /users/michaelbinkley/desktop/GTEx/chr22_subset_gtex_copy2.vcf --show-tags /users/michaelbinkley/desktop/GTEx/GWAScausal.txt --tag-r2 0.5")
    AS.filter_plink_output_gwas(paramDict["gtexDir"] + "filteredplink.txt", paramDict)
    AS.select_rand_variant_gwas2(paramDict["gtexDir"] + "linkedGWAS.txt", paramDict)
    AS.output_pos_rand_variant_gwas2(paramDict["gtexDir"] + "linkedGWASpos.txt", paramDict)
#run_all_steps_H0(paramDict)

def run_all_steps_H1(paramDict):  
    #This selects a 'random' eQTL and GWAS within the selected coldspot
    AS.select_rand_variant(paramDict["gtexDir"] + "eQTLcausal.txt", paramDict)
    AS.output_pos_rand_variant(paramDict["gtexDir"] + "eQTLcausalpos.txt", paramDict)  
    AS.select_rand_variant_gwas(paramDict["gtexDir"] + "GWAScausal.txt", paramDict)
    AS.output_pos_rand_variant_gwas(paramDict["gtexDir"] + "GWAScausalpos.txt", paramDict)
    #This identifies the linked eQTL
    os.system("plink --vcf /users/michaelbinkley/desktop/GTEx/chr22_subset_gtex_copy2.vcf --show-tags /users/michaelbinkley/desktop/GTEx/eQTLcausal.txt --tag-r2 0.5")
    AS.filter_genov_output(paramDict["gtexDir"] + "filteredSNPs.txt", paramDict)
    AS.filter_plink_output(paramDict["gtexDir"] + "filteredplink.txt", paramDict) 
    AS.select_rand_variant2(paramDict["gtexDir"] + "linkedeQTL.txt", paramDict)
    AS.output_pos_rand_variant2(paramDict["gtexDir"] + "linkedeQTLpos.txt", paramDict)   
    AS.prepare_perm_file( "/users/michaelbinkley/desktop/RTCstuffs/tmpSmall/permutations.txt") 
    #This now identifies the linked GWAS variant
    #os.system("plink --vcf /users/michaelbinkley/desktop/GTEx/chr22_subset_gtex_copy.vcf --show-tags /users/michaelbinkley/desktop/GTEx/GWAScausal.txt --tag-r2 0.5")
    AS.filter_plink_output_gwas(paramDict["gtexDir"] + "filteredplink.txt", paramDict)
    AS.select_rand_variant_gwas2(paramDict["gtexDir"] + "linkedGWAS.txt", paramDict)
    AS.output_pos_rand_variant_gwas2(paramDict["gtexDir"] + "linkedGWASpos.txt", paramDict)
#run_all_steps_H1(paramDict)


def run_RTC(paramDict):
    ### RUn this for the real eQTL
    import load_params as LP
    import eqtl_funcs as EF
    import allstepsH as AS
    
    geneName = "ENSG00000172404.4"
    chrNum = "22"
    pos = "40407887"
    tissue = "Breast_Mammary_Tissue"

    slope, intercept,  residuals, indivVector = EF.get_eqtl_stats(chrNum, pos, geneName, tissue, paramDict)  


    
    def calculate_pp(slope_r, intercept, residuals, chrNum, pos, indivVector, paramDict):
 
        permutedResiduals = random.sample(list(residuals), len(residuals))
        genoVector, indivGenoVector = GF.get_dos_vector(chrNum, pos, paramDict["vcfDir"])
        filteredGenoVector, indivVector = GF.filter_and_sort_genotype_vector(genoVector, indivGenoVector, indivVector)
        pp = float(slope_r)*np.array(filteredGenoVector) + permutedResiduals + intercept   
        return pp



   # random.seed(0)

    param_file = "param_files/binkley_local.csv" 

    paramDict = LP.load_param_file(param_file)

    geneName = "ENSG00000172404.4"
    chrNum = "22"
    pos = "40407887"
    tissue = "Breast_Mammary_Tissue"

    #print(pos)
    slope_r, intercept,  residuals, indivVector = EF.get_eqtl_stats(chrNum, pos, geneName, tissue, paramDict)  


    chrNumLinked = "22"
    ### THIS is the position of the randomly selected
    #posLinked = output_pos_rand_variant(paramDict["gtexDir"] + "linkedeQTLpos.txt", paramDict)
    posLinked = AS.output_pos_rand_variant2(paramDict["gtexDir"] + "linkedeQTLpos.txt", paramDict)

    #print("position", posLinked, posLinked[0])
    posLinked = posLinked[0]

    pp = calculate_pp(slope_r, intercept, residuals, chrNumLinked, posLinked, indivVector, paramDict)


    def output_bed(indivVector, expVector, filepath, geneName, paramDict):
        selected = AS.output_gene_info(geneName, paramDict)
        header = AS.output_header2(indivVector)
        with open(filepath, 'w') as outfile:
            outfile.write(header + "\n")
            geneLocation = selected
            outfile.write("\t".join([str(x) for x in geneLocation]) + "\t" + "\t".join([str(x) for x in expVector]) + "\n")
    #print(geneName)        
    AS.output_bed(indivVector, pp, paramDict["tmpSmall"] + "/GTExExpfile.bed", geneName, paramDict)


    #This prepares the index files for the pseudophenotype
    os.system("bgzip /users/michaelbinkley/desktop/RTCstuffs/tmpsmall/GTExExpfile.bed")
    os.system("tabix -p bed /users/michaelbinkley/desktop/RTCstuffs/tmpsmall/GTExExpfile.bed.gz")
    #This calculates the RTC score
    os.system("QTLtools rtc --vcf /users/michaelbinkley/desktop/GTEx/chr22_subset_gtex.vcf.gz --bed /users/michaelbinkley/desktop/RTCstuffs/tmpsmall/GTExExpfile.bed.gz --hotspot /users/michaelbinkley/desktop/qtltools_test/hotspots.bed --gwas-cis /users/michaelbinkley/desktop/GTEx/linkedGWAS.txt /users/michaelbinkley/desktop/RTCstuffs/tmpSmall/permutations.txt --normal --out rtc_results2.txt")

#run_RTC(paramDict)
