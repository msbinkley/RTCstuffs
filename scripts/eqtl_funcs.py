#!/usr/bin/env python3


import gtex_funcs as GF
import ols_funcs as OF
import scipy.stats

def get_eqtl_stats(chrNum, pos, geneName, tissue, paramDict):
    '''

    '''
    expVector, indivExpVector = GF.get_exp_vector(geneName, tissue, paramDict["gtexExpDir"])
    genoVector, indivGenoVector = GF.get_dos_vector(chrNum, pos, paramDict["vcfDir"])
    genoVector = GF.filter_and_sort_genotype_vector(genoVector, indivGenoVector, indivExpVector)

    results = scipy.stats.linregress(genoVector,expVector)
    residuals = OF.calculate_residuals(genoVector, expVector, results.slope, results.intercept)

    return results.slope, results.intercept, residuals
