#!/usr/bin/env python3


import gtex_funcs as GF
import ols_funcs as OF
import scipy.stats
import load_params as LP

param_file = "param_files/binkley_local.csv" 

paramDict = LP.load_param_file(param_file)

def get_eqtl_stats(chrNum, pos, geneName, tissue, paramDict):
    '''

    '''
    #expVector, indivExpVector = GF.get_exp_vector(geneName, tissue, paramDict["gtexExpDir"])
    
    #genoVector, indivGenoVector = GF.get_dos_vector(chrNum, pos, paramDict["vcfDir"])
    expVector, indivExpVector = GF.get_exp_vector(geneName, tissue, paramDict["gtexDir"])
    
    genoVector, indivGenoVector = GF.get_dos_vector(chrNum, pos, paramDict["gtexDir"])
    genoVector = GF.filter_and_sort_genotype_vector(genoVector, indivGenoVector, indivExpVector)

    results = scipy.stats.linregress(genoVector,expVector)
    residuals = OF.calculate_residuals(genoVector, expVector, results.slope, results.intercept)

    return results.slope, results.intercept, residuals
