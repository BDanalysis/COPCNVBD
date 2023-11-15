#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  3 06:40:22 2023

@author: yyj
"""

import numpy as np
import pandas as pd
import pysam
import cvxpy as cp

def readFasta(filename):
    seq = ''
    fread = open(filename)
    #delete row one
    fread.readline()

    line = fread.readline().strip()
    while line:
        seq += line
        line = fread.readline().strip()  
    return seq


#function to read readCount file , generate readCount array
def readRd(filename, seqlen):
    print(seqlen)
    readCount = np.full(seqlen, 0.0)
    samfile = pysam.AlignmentFile(filename, "rb")
    for line in samfile:
        if line.reference_name:
            chr = line.reference_name.strip('chr')
            if chr.isdigit():
                posList = line.positions
                readCount[posList] += 1
        
    return readCount


def read_step1_result(filename):
    region_array = pd.read_csv(filename,names=['chr_name','step1_start','step1_end','type','type_num','hx_start','hx_end'],header=None,sep='\t')
    region_array = region_array[['chr_name','step1_start','step1_end','type','hx_start','hx_end']]
    region_array_len = region_array.shape[0]
    jz_region = region_array[['chr_name','step1_start','step1_end']].values.tolist()
    hx_region = region_array[['chr_name','hx_start','hx_end']].values.tolist()
    return region_array,region_array_len,jz_region,hx_region


def tv_smooth(y, lam):
    n = len(y)
    x = cp.Variable(n)
    obj = cp.Minimize(0.5*cp.sum_squares(x-y) + lam*cp.norm(cp.diff(x), 1))
    prob = cp.Problem(obj)
    prob.solve(solver='OSQP')
    return x.value


def find_hx_pos(arr,threshold):
    grad = np.diff(arr)
    indices = np.where(np.abs(grad)>threshold)[0]+1
    return indices


def find_nearest_boundary(candidates, left_boundary, right_boundary):
    # if the boundary of the candidate CNV region is null，use the boundary of the approximate CNV region
    if len(candidates) == 0:
        return left_boundary, right_boundary

    # find the nearest left boundary
    if left_boundary in candidates:
        cnv_left_boundary = left_boundary
    else:
        left_pointers = [left_boundary - i for i in range(1, len(candidates) + 1)]
        right_pointers = [left_boundary + i for i in range(1, len(candidates) + 1)]
        for l, r in zip(left_pointers, right_pointers):
            if l in candidates:
                cnv_left_boundary = l
                break
            elif r in candidates:
                cnv_left_boundary = r
                break
        else:
            cnv_left_boundary = left_boundary

    # find the nearest right boundary
    if right_boundary in candidates:
        cnv_right_boundary = right_boundary
    else:
        left_pointers = [right_boundary - i for i in range(1, len(candidates) + 1)]
        right_pointers = [right_boundary + i for i in range(1, len(candidates) + 1)]
        for l, r in zip(left_pointers, right_pointers):
            if r in candidates:
                cnv_right_boundary = r
                break
            elif l in candidates:
                cnv_right_boundary = l
                break
        else:
            cnv_right_boundary = right_boundary

    return cnv_left_boundary, cnv_right_boundary


def main(bam,hx_result_path,refpath,chrnum,bam_path,outpath):
    hx_result = hx_result_path + bam + "_hx_result.txt"
    region_array,region_array_len,jz_region,hx_region = read_step1_result(hx_result)
    
    
    binLen = 10
    ref = refpath + 'chr'+str(chrnum)+'.fa'
    chrName = ref.split('.')[0]
    rdFile = bam_path + bam
    outputFile = outpath + bam + '_step2_result.txt'
    
    treeNum = 256
    treeSampleNum = 256
    alpha = 0.25
    
    
    errorNum = 0.005
    seq = readFasta(ref)
    #The length of seq
    seqlen = len(seq)
    print("seqlen:"+str(seqlen))
    rd = readRd(rdFile, seqlen)
    print("rd already finished")
    #fillin n and N position
    rd_mean = np.mean(rd)
    for i in range(seqlen):
        if seq[i] not in ['a', 'A', 't', 'T', 'c', 'C', 'g', 'G']:
            rd[i] = rd_mean
    rd = rd.astype(np.float64)
    print("preprocessing has finished!")
    
    #hx_rd
    result = [[elem for elem in rd[hx_region[i][1]:hx_region[i][2]]] for i in range(len(hx_region))]
    print("rd_result is saved!")
    
    #hx_pos
    result_tv = []
    hx_pos = []
    for i in range(len(result)):
        y = np.array(result[i])
        x_smooth = tv_smooth(y, alpha)
        start = int(jz_region[i][1])
        end = int(jz_region[i][2])
        start_left = start-(end-start+1)
        end_left = start-1
        rd_left = np.mean(rd[start_left:end_left+1])
        rd_jz = np.mean(rd[start:end])
        threshold = np.abs(np.diff(np.array([rd_left,rd_jz])))[0]
        print(threshold)
        hx_pos_i = find_hx_pos(x_smooth,threshold)
        if len(hx_pos_i)==0:
            print("array is empty")
        else:
            hx_pos_i = hx_pos_i + int(hx_region[i][1]) +1
        result_tv.append(x_smooth)
        hx_pos.append(hx_pos_i)
    print("hx_pos is following....")
    print("\n")
    print(hx_pos)
    
    
    #CNVbd
    step2_result = region_array[['chr_name','type']]
    step2_result['step2_start'] = 0
    step2_result['step2_end'] = 0
    
    for j in range(len(hx_pos)):
        jz_left = int(jz_region[j][1])
        jz_right = int(jz_region[j][2])
        step2_leftbd,step2_rightbd = find_nearest_boundary(hx_pos[j], jz_left, jz_right)
        step2_result['step2_start'][j] = step2_leftbd
        step2_result['step2_end'][j] = step2_rightbd
        
    
    step2_result = step2_result[['chr_name','step2_start','step2_end','type']]
    step2_result.to_csv(outputFile,sep = '\t',index=False)
    


if __name__ == '__main__':
    for i in range(1,23):
        chrnum = i
        bam = 'chr' + str(chrnum) + '.bam'
        name = 'chr' + str(chrnum)
        genename = 'EGAR00001004802_2053_1'
        bam_path = '/EGA/'+genename+'/'
        refpath = '/hg37/'
        outpath = '/EGA/'+genename+'/result/copod_step2/'
        step1_path = '/EGA/'+genename+'/result/copod_step1/'
        merge_path = outpath
        hx_result_path = outpath
        
        segpath = outpath +str("/seg_")+name
        p_value_file = outpath + '/' + bam + '.score.txt'
        outfile = outpath + '/'+bam +".result.txt"
        
        main(bam,hx_result_path,refpath,chrnum,bam_path,outpath)




