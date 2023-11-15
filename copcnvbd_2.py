# -*- coding: utf-8 -*-
"""
# Merge two adjacent odd and even rows into a single line, which is the maximum and minimum values of the last four columns, respectively
@author: 24240
"""


import pandas as pd

'''
# read files
with open("query_positions.txt", "r") as f:
    lines = f.readlines()
# Split each row of data into a list
data = [line.strip().split('\t') for line in lines]
'''

def func2(inputpath,bam,outpath):
    data = pd.read_csv(inputpath+bam+"_query_positions.txt", sep='\t',header = None,names=['readID','start', 'end'])
    if len(data) % 2 == 1:
        data = data.drop(data.index[-1])
    # Merge two adjacent singular and even rows into one line, and generate a DataFrame
    df = pd.DataFrame([[data['readID'][i], data['readID'][i+1], data['start'][i],data['start'][i+1], data['end'][i], data['end'][i+1]] for i in range(0, len(data), 2)],
                      columns=['readID_odd', 'readID_even', 'start_odd', 'end_odd', 'start_even', 'end_even'])
    
   
    # Two columns are generated, the minimum and maximum values of the last four columns
    df['min_pos'] = df[['start_odd', 'end_odd', 'start_even', 'end_even']].min(axis=1)
    df['max_pos'] = df[['start_odd', 'end_odd', 'start_even', 'end_even']].max(axis=1)
    
        
    
    df.to_csv(outpath+bam+"_merged_data.txt", sep="\t", index=False,header=False)

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
        inputpath = '/EGA/'+genename+'/result/copod_step2/'
        
        
        segpath = outpath +str("/seg_")+name
        p_value_file = outpath + '/' + bam + '.score.txt'
        outfile = outpath + '/'+bam +".result.txt"
        
        func2(inputpath,bam,outpath)