# COPCNVBD:An Integrated Approach for Copy Number Variation and Breakpints Detection  using Next-generation Sequencing Data

##step1: COPCNV: generate CNV files by the improved COPOD algorithm
##input：bam files: /EGA/chr1.bam, ... , chry.bam
##the reference sequence files: /hg37/chr1.fa, ... , chry.fa
##output: chr*.bam.result.txt (chr_name	CNV_start	CNV_end	CNV_type	CN)
python copod_main.py

##step2: COPCNVBD: generate accurate CNV breakpoints by a PEM-based image edge detection algorithm
##intput: the approximate CNV files: chr*.bam.result.txt
##output: the final CNV files: chr*.bam_step_result.txt
python copcnvbd_1.py
python copcnvbd_2.py
python copcnvbd_3.py
python copcnvbd_4.py
