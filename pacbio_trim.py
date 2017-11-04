# encoding: utf-8
# author: Bo Yan, released in 2017

"""
Include functions used for treating the SMRT-cappable-seq data to remove the adapter and polyA, polyC sequence.
Based on python2.7
Preinstall regex module
"""

try:
    import regex
    import string
    import os
    import time
    import os
    import sys
    import argparse
except ImportError:
    print "python module import error."
    quit()

##-------1. functions

def RC_seq(sequence):
    '''convert input seq to reverse complementary seq'''
    ls = list(sequence)
    ls_new = []
    for letter in ls[::-1]:
        if letter == "A":
            data = "T"
        elif letter == "T":
            data = "A"
        elif letter == "U":
            data = "A"
        elif letter == "G":
            data = "C"
        elif letter == "C":
            data = "G"
        else:
            print "Sequence is not correct."
            exit()
        ls_new.append(data)

    return ''.join(ls_new)

##--------------1.1. Analyze the length of raw read of insert file

def main_insert_length(input_fastq, fa=False):

    print "============="
    print "Analyze the length of the reads."
    
    path = os.path.abspath(input_fastq) 
    dir = os.path.dirname(path) 
    os.chdir(dir)
    print "The output is saved in the dir:", os.getcwd()
    
    basename = os.path.basename(input_fastq).split('.')[0]
    
    input_fastq = open(input_fastq)
    report = open(basename+'.length.report','w') # default file name for the concatemer report
    print "Generate length report: ",  basename+'.length.report'
    
    ## creat bins to save the size of inserts, 1kb,1.5kb...4.5kb
    size = {1:0, 1.5:0, 2:0, 2.5:0, 3:0, 3.5:0, 4:0, 4.5:0}

    count = 0 # number of total reads

    if not fa:
        line = input_fastq.readline() 
        while line:
            line = input_fastq.readline().strip() 
            count +=1
            length = len(line) 
            print>>report, length
            if length<1000: 
                size[1] +=1
            elif length<1500: 
                size[1.5] +=1
            elif length<2000: 
                size[2] +=1
            elif length<2500: 
                size[2.5] +=1
            elif length<3000: 
                size[3] +=1
            elif length<3500: 
                size[3.5] +=1
            elif size<4000: 
                size[4] +=1
            else: 
                size[4.5] +=1

            line = input_fastq.readline() 
            line = input_fastq.readline() 
            line = input_fastq.readline()
    else: 
        line = input_fastq.readline() 
        while line:
            line = input_fastq.readline().strip() 
            count +=1
            length = len(line)
            print>>report, length
            if length<1000: 
                size[1] +=1
            elif length<1500: 
                size[1.5] +=1
            elif length<2000: 
                size[2] +=1
            elif length<2500:
                size[2.5] +=1
            elif length<3000: 
                size[3] +=1
            elif length<3500: 
                size[3.5] +=1
            elif size<4000: 
                size[4] +=1
            else: 
                size[4.5] +=1

            line = input_fastq.readline() 
        

    input_fastq.close()
    report.close()
    
    print "The total number of reads is %d." % count
    print "The number of reads below 1kb is %d." % size[1],
    print "{0:.0f}%".format(float(size[1])/count*100)
    print "The number of reads between 1kb and 1.5kb is %d." % size[1.5],
    print "{0:.0f}%".format(float(size[1.5])/count*100)
    print "The number of reads between 1.5kb and 2.0kb is %d." % size[2],
    print "{0:.0f}%".format(float(size[2])/count*100)
    print "The number of reads between 2.0kb and 2.5kb is %d." % size[2.5],
    print "{0:.0f}%".format(float(size[2.5])/count*100)
    print "The number of reads between 2.5kb and 3.0kb is %d." % size[3],
    print "{0:.0f}%".format(float(size[3])/count*100)
    print "The number of reads between 3kb and 3.5kb is %d." % size[3.5],
    print "{0:.0f}%".format(float(size[3.5])/count*100)
    print "The number of reads between 3.5kb and 4kb is %d." % size[4],
    print "{0:.0f}%".format(float(size[4])/count*100)
    print "The number of reads above 4kb is %d." % size[4.5],
    print "{0:.0f}%".format(float(size[4.5])/count*100)
    print "Generate length report for ploting distribution %s." % basename+'.length.report'
    print "============="
    print
    
    return 1

##--------------1.2. analyze and remove the reads with PCR concatemer

def main_filter_chimera(input_fastq, output, filter):
    print "============="
    print "Analyze the number of transcripts in each read."
    print "Only save the reads that contain one apdapterL or adpaterL_RC seq."
    print "Convert the read with adapterL_RC to its reverse complementary seq."
    
    path = os.path.abspath(input_fastq) 
    dir = os.path.dirname(path) 
    os.chdir(dir)
    print "The output is saved in the dir:", os.getcwd()
    print

    input_fastq = open(input_fastq)
    output_fasta = open(output, 'w') # file for saving the good reads
    filter_fasta = open(filter, 'w') # file for saving the filtered reads
    report = open(output+'.concatemer.report','w') # default file name for the concatemer report

    # matching criteria
    searchL = '(ACACTCTGTCGCTACGTAGATAGCGTTGAGTG){e<=3}' # adapterL seq
    searchL_RC = '(CACTCAACGCTATCTACGTAGCGACAGAGTGT){e<=3}' # adapterL_RC seq

    count_input = 0
    count_output = 0
    line1 = input_fastq.readline().strip() 
    line2 = input_fastq.readline().strip() 
    input_fastq.readline() 
    input_fastq.readline().strip() 

    while line1:
        count_input +=1

        count_L = 0
        count_L_RC = 0
        for item in regex.finditer(searchL, line2, overlapped=False):
            count_L +=1
        for item in regex.finditer(searchL_RC, line2, overlapped=False):
            count_L_RC +=1

        print>>report, count_L + count_L_RC

        if count_L + count_L_RC == 1:
            count_output +=1
            if count_L == 1: # Type1 read
                print>>output_fasta, line1.replace('@','>') # change the ID line to standard fasta format.
                print>>output_fasta, line2
            else: # Type2 read, convert to RC seq into output_fasta
                print>>output_fasta, line1.replace('@','>')
                print>>output_fasta, RC_seq(line2)

        else:
            print>>filter_fasta, line1
            print>>filter_fasta, line2
        
        line1 = input_fastq.readline().strip() 
        line2 = input_fastq.readline().strip() 
        input_fastq.readline() 
        input_fastq.readline().strip() 
        
    input_fastq.close()
    output_fasta.close()
    filter_fasta.close()
    report.close()

    print "There are %d reads in the input fastq file." % count_input
    print "There are %d reads in the output fasta file." % count_output
    print "============="
    print
    return 1


##--------------1.3. extract UID and trim the polyA tail at 3'end
def trimpolyA(seq):
    '''
    find index for the start of polyA tail;
    matching criteria: 'A{5,}[ACTG]{6,11}(AGATCGGAAGAGCACACGTCTGAACTCCAGTCA){e<=3}'
    as well as find the index for UID;
    '''
    if next(regex.finditer('A{5,}[ACTG]{6,11}(AGATCGGAAGAGCACACGTCTGAACTCCAGTCA){e<=3}', seq, overlapped=True), None):
        for item in regex.finditer('A{5,}[ACTG]{6,11}(AGATCGGAAGAGCACACGTCTGAACTCCAGTCA){e<=3}', seq, overlapped=True):
            polyAindex = item.start() 
            break
            
        for item in regex.finditer('A{5,}', seq[polyAindex:], overlapped=True): 
            UIDindex = item.end() 
            
        return polyAindex, seq[polyAindex:][UIDindex+4:UIDindex+10] 
    else:
        return None

def trimpolyC(seq):
    '''
    find index for the end of polyC tail;
    matching criteria: '(ACACTCTGTCGCTACGTAGATAGCGTTGAGTG){e<=3}C{5,}'
    '''
   
    if next(regex.finditer('(ACACTCTGTCGCTACGTAGATAGCGTTGAGTG){e<=3}C{5,}', seq, overlapped=True), None):
        for item in regex.finditer('(ACACTCTGTCGCTACGTAGATAGCGTTGAGTG){e<=3}C{5,}', seq, overlapped=True):
                polyCindex = item.end()
        return polyCindex 
    else:
        return None

    
def main_trim(input_fasta, output, filter):
    print "============="
    localtime = time.asctime(time.localtime(time.time()))
    print "Start at :", localtime
    print "The UID is added into the ID line."
        
    path = os.path.abspath(input_fasta) 
    dir = os.path.dirname(path) 
    os.chdir(dir)
    print "The output is saved in the dir:", os.getcwd()
    print
   
    filter_fasta = open(filter, 'w')

    input_fasta = open(input_fasta)
    output_fasta = open(output,'w')
    
    count_input = 0
    count_output = 0

    line1 = input_fasta.readline().strip() 
    line2 = input_fasta.readline().strip() 

    while line1:
        count_input +=1
        polyA = trimpolyA(line2) 
        polyC = trimpolyC(line2) 
        if polyA and polyC: 
            count_output+=1
            print>>output_fasta, ''.join([line1,'/',polyA[1]]) 
            print>>output_fasta, line2[polyC:polyA[0]]
        else:
            print>>filter_fasta, line1
            print>>filter_fasta, line2

        line1 = input_fasta.readline().strip() 
        line2 = input_fasta.readline().strip() 

    input_fasta.close()
    output_fasta.close()
    filter_fasta.close()

    print "There are %d reads in the input file." % count_input
    print "There are %d reads in the output file, which contain UID." % count_output
    localtime = time.asctime(time.localtime(time.time()))
    print "End at :", localtime
    print "============="
    print 
   
    return 1

##-------2. mainbody
   
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help='sub-command help', dest = 'mode')
    
    ##--------------1.1. main_insert_length()
    parser_1 = subparsers.add_parser('length', help='analyze the length of reads')
    parser_1.add_argument('--input', help='input fastq file', dest='input_fastq')
    parser_1.add_argument('--fa', action = 'store_true', help='add this tag for fasta file', dest='fa')

    ##--------------1.2. main_filter_chimera()
    parser_2 = subparsers.add_parser('filter', help='analyze and filter chimera')
    parser_2.add_argument('--input', '-i', help='input fastq file', dest='input_fastq')
    parser_2.add_argument('--output', '-o', help='output fasta file', dest='output')
    parser_2.add_argument('--filter', '-f', help='output fasta file', dest='filter')

    ##--------------1.3.main_trim(input_fasta)
    parser_3 = subparsers.add_parser('poly', help='trim polyA and polyC tails, extract UID')
    parser_3.add_argument('--input', '-i', help='input fasta file', dest='input_fasta')
    parser_3.add_argument('--output', '-o', help='output fasta file', dest='output')
    parser_3.add_argument('--filter', '-f', help='output fasta file', dest='filter')

    
    args = parser.parse_args()

    if args.mode =='length':
        main_insert_length(args.input_fastq, args.fa)
        
    elif args.mode == 'filter':
        main_filter_chimera(args.input_fastq, args.output, args.filter)

    elif args.mode == 'poly':
        main_trim(args.input_fasta, args.output, args.filter)
