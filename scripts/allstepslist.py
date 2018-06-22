
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

