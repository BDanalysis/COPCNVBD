
def func3(step1_path,bam,merge_path,outpath):

    # read the first txt file
    with open(step1_path + bam + ".result.txt") as f1:
        lines1 = f1.readlines()
    
    # read the second txt file
    with open(merge_path + bam + "_merged_data.txt") as f2:
        lines2 = f2.readlines()
    
    new_lines = []
    
    # traverse each row in the first txt file
    for line1 in lines1:
        # find the start and end position of the first txt file
        chr1, start1, end1, type1, value1 = line1.strip().split("\t")
        start1 = int(start1)
        end1 = int(end1)
        min_pos = start1
        max_pos = end1
    
        # traverse each row in the second txt file
        for line2 in lines2:
            # find the start and end position of the second txt file
            _, _, _, _, _, _, start2, end2 = line2.strip().split("\t")
            start2 = int(start2)
            end2 = int(end2)
    
            # judge whether there are overlaps between the start and end positions of the second file and those of the first file
            if start2 <= end1 and end2 >= start1:
                # find the min and the max position
                min_pos = min(min_pos, start2)
                max_pos = max(max_pos, end2)
    
        new_line = line1.strip() + "\t" + str(min_pos) + "\t" + str(max_pos) + "\n"
        new_lines.append(new_line)
    
    with open(outpath + bam + "_hx_result.txt", "w") as f:
        f.writelines(new_lines)

if __name__ == '__main__':
    for i in range(1,23):
        chrnum = i
        bam = 'chr' + str(chrnum) + '.bam'
        name = 'chr' + str(chrnum)
        genename = 'EGAR00001004802_2053_1'
        bam_path = '/EGA/'+genename+'/'
        refpath = '/hg37'
        outpath = '/EGA/'+genename+'/result/copod_step2/'
        step1_path = '/EGA/'+genename+'/result/copod_step1/'
        merge_path = outpath
        
        
        segpath = outpath +str("/seg_")+name
        p_value_file = outpath + '/' + bam + '.score.txt'
        outfile = outpath + '/'+bam +".result.txt"
        
        func3(step1_path,bam,merge_path,outpath)
 