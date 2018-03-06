# encoding: utf-8
# author Bo Yan, released in 2017.

"""
Contain functions for SMRT cappable-seq TSS analysis.
Based on Python 2.7.

Usage:
python TSS_analysis.py count --input inputbed --output count.output
python TSS_analysis.py cluster --input count.output --output --control --cutoff (default5)

Note:
Currently this script is only designed to count the reads mapped to only one genome, 
in other words 'cut -f1 inputbed | sort -u | uniq -c' == 1

"""
try: 
    from collections import deque
    from subprocess import check_call
    from subprocess import check_output
    import os
    import bisect
    import argparse
    import time
except:
    print "Import module error."
    quit()

##-------1. functions

##--------------1. count the number of reads with the same 5'end (TSS) and 3'end
def countreads(input_bed, output_bed, i):
    if os.path.isfile('countreads.temp.input'):
        print "countreads.temp.input file exists!"
        quit()
     
    ## sort based on required column
    if i==1: # positive strand
        command = 'sort -k2,2n ' + input_bed + '>countreads.temp.input'
        check_call(command, shell = True)
        pos5=1 # python index for TSS
        pos3=2 # python index for 3'end
        strand = "+"
    
    elif i==2: # negative strand
        command = 'sort -k3,3n ' + input_bed + '>countreads.temp.input'
        check_call(command, shell = True)
        pos5=2 # python index for TSS
        pos3=1 # python index for 3'end
        strand = '-'
    else:
        print "Error! Choose 1 for positive strand or 2 for negative strand."
        quit()
    
    
    input_bed = open('countreads.temp.input')
    output= open(output_bed, 'w') 
    output_largest = open(output_bed + '.largest', 'w') 
    q = deque(maxlen=2)
    line = input_bed.readline()
    q.append(line.strip().split("\t"))

    count = 1 
    dic = dict([(q[0][pos3],1)]) 
    for line in input_bed:
        
        q.append(line.strip().split("\t"))
        if q[1][pos5]==q[0][pos5]: # both the chr and TSS are the same

            count +=1
            if dic.has_key(q[1][pos3]):
                dic[q[1][pos3]]+=1 
            else:
                dic[q[1][pos3]]=1
        else:
            
            if i==1:
                largest=str(max(int(item) for item in dic.keys())) 
            else:
                largest=str(min(int(item) for item in dic.keys())) 
            
            print>>output_largest, '\t'.join(str(x) for x in [q[0][0], q[0][pos5], largest , 'max', '.', strand, count, dic[largest]])
            
            for item in dic.keys():
                print>>output, '\t'.join(str(x) for x in [q[0][0], q[0][pos5], item, largest, '.', strand, count, dic[item]])

            dic = dict([(q[1][pos3],1)])
            count = 1

    if i==1:
        largest=str(max(int(item) for item in dic.keys()))
    else:
        largest=str(min(int(item) for item in dic.keys()))
    print>>output_largest, '\t'.join(str(x) for x in [q[1][0], q[1][pos5], largest , 'max', '.', strand, count, dic[largest]])
    for item in dic.keys():
        print>>output, '\t'.join(str(x) for x in [q[1][0], q[1][pos5], item, largest, '.', strand, count, dic[item]])
            
    input_bed.close()
    output.close()
    output_largest.close()
    
    os.remove('countreads.temp.input')

    return 1

def main_countreads(input_bed, output_bed):
    '''
    Use to count the number of reads at the same TSS(colum2 for +; column3 for -);
    and the number of reads at the same TSS as well as 3'end(column3 for -; column2 for +). 
    '''
    print "============="
    localtime = time.asctime(time.localtime(time.time()))
    print "Start at :", localtime
    print "Counting the number of reads with the same 5'end and 3'end."
       
    path = os.path.abspath(output_bed)
    dir = os.path.dirname(path) 
    os.chdir(dir)
    print "The output is saved in the dir:", os.getcwd()
    
    with open(input_bed) as f:
        chr = [line.strip().split('\t')[0] for line in f]
        
    if len(set(chr)) >1:
        print "There are more than one genome used for mapping, please extract the reads that are only mapped to the targeting bacterial genome!"
        print "Quit."
        quit()
    
    # forward strand 
    print "Counting positive strand."
    command = "awk -v OFS=\'\\t\' \'{if($1==\"" + chr[0] + "\" && $6==\"+\"){print $0}}\' " + input_bed + '>countread.positive.bed'
    check_call(command, shell = True)
    countreads('countread.positive.bed', 'countread.positive.bed.output', 1) # also generate countread.positive.bed.output.largest by default


    # reverse strand
    print "Counting negative strand."
    command = "awk -v OFS=\'\\t\' \'{if($1==\"" + chr[0] + "\" && $6==\"-\"){print $0}}\' " + input_bed + '>countread.negative.bed'
    check_call(command, shell = True)
    countreads('countread.negative.bed', 'countread.negative.count.bed', 2) # also generate countread.negative.count.bed.largest by default
   
        
    ## rearragne the column2 and column3 for negative strand, the final output is as the standard bed file.
    command = "awk -v OFS=\'\\t\' \'{print $1,$3,$2,$4,$5,$6,$7,$8}\' countread.negative.count.bed > countread.negative.bed.output"
    check_call(command, shell = True)
    command = "awk -v OFS=\'\\t\' \'{print $1,$3,$2,$4,$5,$6,$7,$8}\' countread.negative.count.bed.largest > countread.negative.bed.output.largest"
    check_call(command, shell = True)
    
    
    ## combine two strands
    command = "cat countread.negative.bed.output >" + output_bed
    check_call(command, shell = True)
    command = "cat countread.positive.bed.output >>" + output_bed
    check_call(command, shell = True)
    
    command = "cat countread.negative.bed.output.largest >" + output_bed + ".largest"
    check_call(command, shell = True)
    command = "cat countread.positive.bed.output.largest >>" + output_bed + ".largest"
    check_call(command, shell = True)
    
    print "Also generate a file to save the number for the longest read for each TSS."
    ## delete tmp
   
    print "Removing the tempory files."
    ls = ['countread.positive.bed', 'countread.negative.bed', 
    'countread.negative.count.bed', 'countread.negative.count.bed.largest',
    'countread.positive.bed.output', 'countread.negative.bed.output',
    'countread.positive.bed.output.largest', 'countread.negative.bed.output.largest']
    for item in ls:
        os.remove(item)
    
    
    localtime = time.asctime(time.localtime(time.time()))
    print "End at :", localtime
    print "============="
    return 1

        
##--------------1.2. cluster TSS
def find_most(ls):
    ls_number = [int(item[6]) for item in ls]
    return ls[ls_number.index(max(ls_number))]
    
 
def cluster_TSS(input_file, output_file, cutoff, i):
    if i==1: # positive strand
        command = 'sort -k2,2n ' + input_file + '>cluster.temp.input' 
        check_call(command, shell = True)
        pos5=1 # python index for TSS
        pos3=2 # python index for 3'end
        strand = '+'
    
    elif i==2: # negative strand
        command = 'sort -k3,3n ' + input_file + '>cluster.temp.input'
        check_call(command, shell = True)
        pos5=2 # python index for TSS
        pos3=1 # python index for 3'end
        strand = '-'

    input_bed = open('cluster.temp.input')
    output_control = open(output_file+'.cluster_control','w') 
    output_file = open(output_file,'w')
    
    q = deque(maxlen=2)
    line = input_bed.readline()
    q.append(line.strip().split("\t")) 
    temp = [] 
    temp.append(q[0]) 
    count_control=0
    
    for line in input_bed:
        q.append(line.strip().split("\t"))
        
        if int(q[1][pos5]) - int(q[0][pos5]) <=cutoff: 
            temp.append(q[1])
            
        else:
            TSS = find_most(temp)[pos5] 
            TSS_original = list(set([str(x[pos5]) for x in temp]))
            print>>output_control, TSS, '\t', ';'.join(TSS_original), '\t',strand
            for item in temp:
                item[pos5] = TSS 
                times = int(item[7]) 
                for i in range(0, times):
                    print>>output_file, "\t".join(item[:6]) 
            temp = []
            temp.append(q[1])
          
    TSS = find_most(temp)[pos5] 
    TSS_original = list(set([str(x[pos5]) for x in temp]))
    print>>output_control, TSS, '\t', ';'.join(TSS_original), '\t',strand
    for item in temp:
        item[pos5] = TSS 
        times = int(item[7]) 
        for i in range(0, times):
            print>>output_file, "\t".join(item[:6])

    output_file.close()
    input_bed.close()
    output_control.close()
    os.remove('cluster.temp.input')
    return 1
   
def main_cluster_TSS(input_file, output_file, output_control, cutoff):
    '''
    Cluster the nearby TSSs [start at T1, T2 ..Tn] if abs(T2-T1)<=cutoff ... abs(Tn -Tn-1)<=cutoff;
    in other words, the distance between TSSs in the Cluster <=cutoff;
    then change all the position of TSS in the Cluster to the TSS with maiximum number of reads; 
    output all the original reads with the modified TSS and the original 3'end;
    So there are the same number of reads in the output file with modified TSS as in the original bed file before counting.
    Then the number of new TSS can be re counted using countread.py to count the number of TSS and the same reads after cluster.
    '''
    
    print "============="
   
    path = os.path.abspath(input_file) 
    dir = os.path.dirname(path) 
    os.chdir(dir)
    print "The current working dir is:", os.getcwd()
    print
    print "Making the standard bed8 file for clutering for both strands."
    input_file = open(input_file)
    positive_temp = open('cluster_positive.temp','w')
    negative_temp = open('cluster_negative.temp','w')
    
    for line in input_file:
        line = line.strip().split('\t')
        
        if line[0]=='NC_000913.3' and line[5] == "-": 
            temp = line[0:3] 
            temp.append('.')
            temp.append('.')
            temp.append('-') 
            temp.append(line[6]) 
            temp.append(line[7]) 
            print>>negative_temp, "\t".join(temp)
        elif line[0]=='NC_000913.3' and line[5] == "+": 
            temp = line[0:3] 
            temp.append('.')
            temp.append('.')
            temp.append('+') 
            temp.append(line[6]) 
            temp.append(line[7]) 
            print>>positive_temp, "\t".join(temp)
    input_file.close()
    positive_temp.close()
    negative_temp.close()
    
    ## cluster
    print "Clustering for the positive strand."
    cluster_TSS('cluster_positive.temp', 'cluster_positive.output', cutoff, i=1)
    print "Clustering for the negative strand."
    cluster_TSS('cluster_negative.temp', 'cluster_negative.output', cutoff, i=2)
   
    ## combine two strands
    command = "cat cluster_negative.output >" + output_file
    check_call(command, shell = True)
    command = "cat cluster_positive.output >>" + output_file
    check_call(command, shell = True)
    
    ## combine two strands for output file saving the number of reads for control at certain TSS
    command = "cat cluster_negative.output.cluster_control >" + output_control
    check_call(command, shell = True)
    command = "cat cluster_positive.output.cluster_control >>" + output_control
    check_call(command, shell = True)
   
    print "Removing the tempory files."
    ls = ['cluster_positive.temp', 'cluster_negative.temp', 
    'cluster_positive.output', 'cluster_negative.output',
    'cluster_negative.output.cluster_control', 'cluster_positive.output.cluster_control']
    for item in ls:
        os.remove(item)
    
    print "Done."
    print "============="


##--------------Parser
parser = argparse.ArgumentParser()

subparsers = parser.add_subparsers(help='sub-command help', dest = 'mode')

##--------------1.1. main_countreads(input_bed, output_bed)
parser_1 = subparsers.add_parser('count', help='count the number of reads with same TSS and TTS')
parser_1.add_argument('--input','-i', help='input bed file from enrich', dest='input_file')
parser_1.add_argument('--output','-o', help='output bed file', dest='output_file')

##--------------1.2. main_cluster_TSS(input_file, output_file, output_control, cutoff)
parser_2 = subparsers.add_parser('cluster', help='cluster the nearby TSS')
parser_2.add_argument('--input','-i', help='input bed file from enrich', dest='input_file')
parser_2.add_argument('--output','-o', help='output bed file', dest='output_file')
parser_2.add_argument('--control','-con', help='output saving the origianl TSS in the cluster', dest='output_control')
parser_2.add_argument('--cutoff','-c', help='cutoff for clustering', type=int, default=5, dest='cutoff')

if __name__ == '__main__':
    args = parser.parse_args()

    if args.mode =='count':
        main_countreads(args.input_file, args.output_file)
    elif args.mode == 'cluster':
        main_cluster_TSS(args.input_file, args.output_file, args.output_control, args.cutoff)
