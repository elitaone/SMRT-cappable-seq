# encoding: utf-8
# author Bo Yan, released in 2017

"""
Use to adjust the 3'end for reads in the bed file after alignment;
to remove the matching due to the incomplete removal of the addition of polyA tail in the previous filter step.
Here all the input reads are mapped to only one bacterial genome.

Based on python2.7

Need bedtools pre-installed.

Usage:
python adjust_3end.py --information file.txt

--information: is a file containing all the information requried for this example. The file.txt is as following:
Input:/Users/yan/Documents/script/python_script/bioinformatics_function/pacbio_tool/Github/enrich_161020.bed
Output:enrich_161020.adjust.bed
Bedtools Path:/usr/local/bin/
Reference:/Users/yan/Documents/script/python_script/bioinformatics_function/pacbio_tool/reference/NC_000913.3/GCF_000005845.2_ASM584v2_genomic.fna

See explaination of this function in README.txt
"""

try:
    from subprocess import check_call, check_output 
    import argparse
    import re
    import os
except ImportError:
    print "python import module Error!"


##-------1. create a script for bedtools getfasta
def getfasta(input_bed, reference, bedtoolspath):
    '''
    creat sh script for using bedtools getfasta
    reference if the E coli genome fasta file.
    '''
    with open('my_bedtoolscript.sh', 'w') as script:
        print>>script, '#!/bin/sh'
        print>>script, '{}bedtools getfasta -s -fi {} -bed {} -fo {}'.format(bedtoolspath, reference, input_bed, 'adjust.temp.fasta')
    check_call('sh my_bedtoolscript.sh', shell=True)
    os.remove('my_bedtoolscript.sh')
    return 1
    
##-------2. modify the bed file based on the fasta    

def modify(input_bed, output, input_fasta):
    '''
    Use to adjust the 3'end for reads in the bed file after alignment;
    to remove the unsure matching due to the addition of polyA tail.
    '''
        
    input_fasta = open(input_fasta)
    input_bed = open(input_bed)
    output = open(output,'w')

    input_fasta.readline()
    line_fasta = input_fasta.readline().strip() 

    for line in input_bed:
        line = line.strip().split('\t')
        if line[-1] == "+":
            if re.findall('A{1,}$',line_fasta):
                line[2] = int(line[2]) - len(re.findall('A{1,}$',line_fasta)[0])
                line[2] = str(line[2])
                line[4] = 'mod3end' 
                print>>output, '\t'.join(line)
            else:
                print>>output, '\t'.join(line)
        elif line[-1] == "-":
            if re.findall('A{1,}$',line_fasta): 
                line[1] = int(line[1]) + len(re.findall('A{1,}$',line_fasta)[0]) 
                line[1] = str(line[1])
                line[4] = 'mod3end' 
                print>>output, '\t'.join(line)
            else:
                print>>output, '\t'.join(line)
                    
        input_fasta.readline()
        line_fasta = input_fasta.readline().strip()

    input_fasta.close()
    input_bed.close()
    output.close() 
    
    return 1
    

##-------Parser
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--information', help='information file', dest='information')
    args = parser.parse_args()
    
    print "======================"
    print "Adjust the 3'end if the matching is due to the addition of A tail."
    
    print "Check the information."

    with open(args.information) as f:
        input_file = f.readline().strip().split(':')[-1].strip()
        print "input file:", input_file
        output_file = f.readline().strip().split(':')[-1].strip()
        print "output file:", output_file
        bedtoolspath = f.readline().strip().split(':')[-1].strip()
        print "bedtools path:", bedtoolspath
        reference = f.readline().strip().split(':')[-1].strip()
        print "reference genome:", reference
      
    path = os.path.abspath(output_file)
    dir = os.path.dirname(path) 
    os.chdir(dir)
    print "The output is saved in the dir:", os.getcwd()
    
    try:
        with open(input_file) as f:
            chr = [line.strip().split('\t')[0] for line in f]
    except:
        print "Input file does not exist."
        print "Quit."
        quit()
    
    if len(set(chr)) >1:
        print "input file contains reads mapped to more than one genome!"
        print "Quit."
        quit()
        
    command = "awk -v OFS=\'\\t\' \'{if($1==\"" + chr[0] + "\"){print $0}}\' " + input_file + '>adjust.temp.bed'
    print command
    check_call(command, shell = True)

    print "Get the fasta seq.\n"
    try:
        getfasta('adjust.temp.bed', reference, bedtoolspath)
    except:
        print "Error generated from using bedtools getfasta."
        print "Double check the reference file and bedtools path."
        quit()
    
    
    print "Add mod3end tag in column5 for the modified reads."
    modify('adjust.temp.bed', output_file, 'adjust.temp.fasta')
    
    os.remove('adjust.temp.fasta')
    os.remove('adjust.temp.bed')
    print "Done."    
    print "======================"

