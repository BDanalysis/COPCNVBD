# -*- coding: utf-8 -*-
"""
#traverse and sort all reads, recording the start and end positions of each read
@author: yyj
"""

import pysam
def func1(bam_path,bam,outpath):
    bam_file = pysam.AlignmentFile(bam_path+bam, "rb")
    
    query_positions = {}
    for read in bam_file:
    #for read in bam_file:
        query_positions[read.query_name] = (read.reference_start, read.reference_end)
    
    with open(outpath+bam+"_query_positions.txt", "w") as f:
        for query_name in sorted(query_positions.keys()):
            start, end = query_positions[query_name]
            f.write(f"{query_name}\t{start}\t{end}\n")
            
if __name__ == '__main__':
    
    for i in range(1,23):
    #bam = 'chr'+str(i)+'.bam'
    #name = 'chr'+str(i)
        chrnum = i
        bam = 'chr' + str(chrnum) + '.bam'
        name = 'chr' + str(chrnum)
        genename = 'EGAR00001004802_2053_1'
        bam_path = '/EGA/'+genename+'/'
        refpath = '/hg37'
        outpath = '/EGA/'+genename+'/result/copod_step2/'
        
        
        segpath = outpath +str("/seg_")+name
        p_value_file = outpath + '/' + bam + '.score.txt'
        outfile = outpath + '/'+bam +".result.txt"
        
        func1(bam_path,bam,outpath)
