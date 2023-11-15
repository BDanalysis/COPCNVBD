#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 17:09:24 2022

@author: yyj
"""
import copod_preprocess

if __name__ == "__main__":
    
    binSize = 1000
    col = 50
    k = 10

    for i in range(1,23):
    #bam = 'chr'+str(i)+'.bam'
    #name = 'chr'+str(i)
        chrnum = i
        bam = 'chr' + str(chrnum) + '.bam'
        name = 'chr' + str(chrnum)
        genename = 'EGAR00001004802_2053_1'
        bam_path = '/EGA/'+genename+'/'
        refpath = '/hg37'
        outpath = 'EGA/'+genename+'/result/copod_step1'
        
        
        segpath = outpath +str("/seg_")+name
        p_value_file = outpath + '/' + bam + '.score.txt'
        outfile = outpath + '/'+bam +".result.txt"
        params = (bam_path, bam, refpath, outpath, segpath, binSize, col, k, chrnum)
        copod_preprocess.main(params)
