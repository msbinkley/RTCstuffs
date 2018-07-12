
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

#geneName = "ENSG00000172404.4"
#chrNum = "22"
#pos = "40407887"
#tissue = "Breast_Mammary_Tissue"

def run_perm_H1(paramDict, geneName, chrNum, pos, tissue, numperm):
    numperm = numperm
    outlist = list()

    for i in range(numperm):
      
        ASL.run_all_steps_H1(paramDict)
        rrtc.run_RTC(paramDict, geneName, chrNum, pos, tissue)
        
        fileIN = open('/users/michaelbinkley/desktop/RTCstuffs/rtc_results2.txt')
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
#run_perm_H1(paramDict, geneName, chrNum, pos, tissue)



def run_perm_H0(paramDict, geneName, chrNum, pos, tissue, numperm):
    numperm = numperm
    outlist2 = list()

    for i in range(numperm):
      
        ASL.run_all_steps_H0(paramDict)
        rrtc.run_RTC(paramDict, geneName, chrNum, pos, tissue)

        
        fileIN = open('/users/michaelbinkley/desktop/RTCstuffs/rtc_results2.txt')
        fileIN.readline()
        genes = list()
        for j in fileIN:
            g = j.rstrip().split(" ")
            name = g[19]

            #genes.append(name)
            outlist2.append([geneName, chrNum, pos, tissue, name])
            #first = (j[0]) 
        
        fileIN.close()
    return outlist2



