#!/usr/local/bin/python3.6
import os
import sys
import glob
import re


py_script = sys.argv[0]
sample_file = sys.argv[1]


print(sys.argv, len(sys.argv))

def main():
    with open(sample_file, "r") as f:
        sample=f.read().rstrip()
    print(sample)
    files="data/P53*trimmed.fastq"
    print(files)
    
    output_name=get_output_name(sample)
    print(output_name)
    IDR1R2=combine_ID_read(files,output_name)
    for item in output_name:
        output = 'data/' + output_name[str(item)] 
        
        with open(output,"w") as f:
            for l in range(len(IDR1R2)):
                if str(item) == IDR1R2[l][0]:
                    f.write("%s\t%s\t%s\n" % (IDR1R2[l][1],IDR1R2[l][2],IDR1R2[l][3]))
    
####### create a dictionary with barcodes as keys and output file names as value.
output_name={}
def get_output_name(sample):
    sample_item=sample.split('\n')
    for item in sample_item:
        ID=item.split('\t')
        output_name[ID[0]]=ID[1]
    return output_name

##### format(ID, read1, read2)
IDR1R2=[]
def combine_ID_read(files,output_name):
    for name in glob.glob(files):
        inputfilename = str(name)
        if "R1" in inputfilename:
            with open(inputfilename, 'r') as f_R1:
                data_R1 = f_R1.read().rstrip().split('\n')
                print(len(data_R1))
        elif "R2" in inputfilename:
            with open(inputfilename, 'r') as f_R2: 
                data_R2 = f_R2.read().rstrip().split('\n') 
                print(len(data_R2))
        else:
            print("file is not found")    
    for i in range(0,(len(data_R1)),4):
        match = re.search(r"^@([^:]+):([0-9]+):([^:]+):([0-9]+):([0-9]+):([0-9]+):([0-9]+)", data_R1[i])
        if match:
            read1_ID = match.group()
            read1_seq = data_R1[i+1]
            if read1_ID in data_R2[i]:
                read2_seq = data_R2[i+1]
                barcode_seq = "1" + read1_seq[0:8] +"2" +read2_seq[0:8]
                if barcode_seq in output_name:
                    IDR1R2.append([barcode_seq,read1_ID,read1_seq,read2_seq])
            else:
                print("read1 does not match with read2")
        else:
            print ("no ID")
    print(len(IDR1R2))
    return IDR1R2
    
main()    
