#!/usr/local/bin/python3.6

import glob
from operator import itemgetter

# check format of input files
# check if P53_Miseq_screen.txt is pre-exist
# check the filter for mutation rate

exon2_F = "ATGACTGCCATGGAGGAGTCACAGTCGGATATCAGCCTCGAGCTCCCTCTGAGCCAGGAGACATTTTCAGGCTTATGGAAACT"
exon2_R = "AGTTTCCATAAGCCTGAAAATGTCTCCTGGCTCAGAGGGAGCTCGAGGCTGATATCCGACTGTGACTCCTCCATGGCAGTCAT"

exon3_F = "GACTTCCTCCAGAAGATATCCTG"
exon3_R = "CAGGATATCTTCTGGAGGAAGTC"

exon4_F = "GCCATCACCTCACTGCATGGACGATCTGTTGCTGCCCCAGGATGTTGAGGAGTTTTTTGAAGGCCCAAGTGAAGCCCTCCGAGTGTCAGGAGCTCCTGCAGCACAGGACCCTGTCACCGAGACCCCTGGGCCAGTGGCCCCTGCCCCAGCCACTCCATGGCCCCTGTCATCTTTTGTCCCTTCTCAAAAAACTTACCAGGGCAACTATG"
exon4_R = "ACCGTGCACATAACAGACTTGGCTGTCCCAGACTGCAGGAAGCCCAGGTGGAAGCCATAGTTGCCCTGGTAAGTTTTTTGAGAAGGGACAAAAGATGACAGGGGCCATGGAGTGGCTGGGGCAGGGGCCACTGGCCCAGGGGTCTCGGTGACAGGGTCCTGTGCTGCAGGAGCTCCTG"                              

exon5_F ="AGTACTCTCCTCCCCTCAATAAGCTATTCTGCCAGCTGGCGAAGACGTGCCCTGTGCAGTTGTGGGTCAGCGCCACACCTCCAGCTGGGAGCCGTGTCCGCGCCATGGCCATCTACAAGAAGTCACAGCACATGACGGAGGTCGTGAGACGCTGCCCCCACCATGAGCGCTGCTC"
exon5_R ="ACCATCACCATCGGAGCAGCGCTCATGGTGGGGGCAGCGTCTCACGACCTCCGTCATGTGCTGTGACTTCTTGTAGATGGCCATGGCGCGGACACGGCTCCCAGCTGGAGGTGTGGCGCTGACCCACAACTGCACAGGGCACG" 

exon6_F = "GCCTGGCTCCTCCCCAGCATCTTATCCGGGTGGAAGGAAATTTGTATCCCGAGTATCTGGAAGACAGGCAGACTTTTCGCCACAGCGTGGTGGTACCTTATGAGCCACCCGAG"
exon6_R = "CTCGGGTGGCTCATAAGGTACCACCACGCTGTGGCGAAAAGTCTGCCTGTCTTCCAGATACTCGGGATACAAATTTCCTTCCACCCGGATAAGATGCTGGGGAGGAGCCAGGC"

exon7_F = "GCCGGCTCTGAGTATACCACCATCCACTACAAGTACATGTGTAATAGCTCCTGCATGGGGGGCATGAACCGCCGACCTATCCTTACCATCATCACACTGGAAGACTCCAG"
exon7_R = "CTGGAGTCTTCCAGTGTGATGATGGTAAGGATAGGTCGGCGGTTCATGCCCCCCATGCAGGAGCTATTACACATGTACTTGTAGTGGATGGTGGTATACTCAGAGCCGGC"

exon8_F = "TGGGAACCTTCTGGGACGGGACAGCTTTGAGGTTCGTGTTTGTGCCTGCCCTGGGAGAGACCGCCGTACAGAAGAAGAAAATTTCCGCAAAAAGGAAGTCCTTTGCCCTGAACTGCCCCCAGGGAGCGCAAAGAGAG"
exon8_R = "CTCTCTTTGCGCTCCCTGGGGGCAGTTCAGGGCAAAGGACTTCCTTTTTGCGGAAATTTTCTTCTTCTGTACGGCGGTCTCTCCCAGGGCAGGCACAAACACGAACCTCAAAGCTGTCCCGTCCCAGAAGGTTCCCA"

exon9_F =  "CGCTGCCCACCTGCACAAGCGCCTCTCCCCCGCAAAAGAAAAAACCACTTGATGGAGAGTATTTCACCCTCAAG"
exon9_R =  "CTTGAGGGTGAAATACTCTCCATCAAGTGGTTTTTTCTTTTGCGGGGGAGAGGCGCTTGTGCAGGTGGGCAGCG"

exon10_F = "ATCCGCGGGCGTAAACGCTTCGAGATGTTCCGGGAGCTGAATGAGGCCTTAGAGTTAAAGGATGCCCATGCTACAGAGGAGTCTGGAGACAGCAGGGCTCACTCCAG"
exon10_R = "CTGGAGTGAGCCCTGCTGTCTCCAGACTCCTCTGTAGCATGGGCATCCTTTAACTCTAAGGCCTCATTCAGCTCCCGGAACATCTCGAAGCGTTTACGCCCGCGGAT"

exon11_F = "CTACCTGAAGACCAAGAAGGGCCAGTCTACTTCCCGCCATAAAAAAACAATGGTCAAGAAAGTGGGGCCTGACTCAGACT"
exon11_R = "AGTCTGAGTCAGGCCCCACTTTCTTGACCATTGTTTTTTTATGGCGGGAAGTAGACTGGCCCTTCTTGGTCTTCAGGTAG"

# If some exons have low coverage, perform this function to check coverage
exon=input("If some exons have low coverage, enter exon name and use comma to seperate each exon, e.g Exon4,Exon8: " )

def main():
    
    for name in glob.glob('data/1*2*.txt'):                            
        inputfilename = str(name)
        print(inputfilename)
        outputfilename1 = "data/Top100_seq" + str(name[20:-4])+".txt"
        outputfilename2 = "data/Check_seq_all.txt"
        
        # input one sample each time
        with open(inputfilename, 'r') as f:                         
            data_in = f.read().rstrip().split('\n')
        
        # at the beginning, counted_line have 3 columns [seq, count_read, Percentage in total reads]                                      
        counted_line = get_unique(data_in)
        
        # start a new dictionary.
        exon_reads = {"Exon2":0, "Exon3":0, "Exon4":0, "Exon5":0, "Exon6":0, "Exon7":0, "Exon8":0, "Exon9":0, "Exon10":0, "Exon11":0, "unknown":0} 
        
        # renewed counted_line have 4 columns [seq, count_read, Percentage in total reads,Exon]
        counted_line, exon_reads  = count_exon_reads(counted_line,exon_reads)
        print(len(counted_line))
        
        # renewed counted_line have 5 columns [seq, count_read, Percentage in total reads, Exon, Percentage in each exon]
        counted_line = filter_reads(counted_line, exon_reads)
        print(len(counted_line))
            
        # output a summary file all samples including top 100 reads for each sample
        with open("data/P53_Miseq_screen_genotyping.txt", 'a') as f:
            f.write("{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}\n".format("Sample ID","\t","Order","\t","Sequence","\t","Reads","\t","Percentage in total reads","\t","Exon","\t","Percentage in each exon","\t","Genotype"))
            for i in range(0,300):
                f.write("{}{}{}{}{}{}{:d}{}{:.2%}{}{}{}{:.2%}{}{}\n".format(inputfilename[2:],"\t",str(i+1),"\t",counted_line[i][0],"\t",counted_line[i][1],"\t",counted_line[i][2],"\t",counted_line[i][3],"\t",counted_line[i][4],"\t",counted_line[i][5]))  
        # output top 50 reads for each sample in fasta format
        with open(outputfilename1, 'w') as f:
            for i in range(0,100):
                f.write(">"+inputfilename[20:-4]+"_"+str(i+1)+"\n") 
                f.write(counted_line[i][0]+"\n")
        # output filtered reads (>1%) for each sample in fasta format. The actual mutation rate is about 4% considering the total reads including read1,read2, and many bad reads.
        with open(outputfilename2, 'a') as f:
            for i in range(0,300):
                if  counted_line[i][4] > 0.01 and counted_line[i][5] == "Check":
                    f.write(">"+inputfilename[20:-4]+"_"+str(i+1)+"\n")
                    f.write(counted_line[i][0]+"\n")
        # output reads from low-coverage exons
        if "Exon" in exon:
            ex_seq= calculate_coverage(counted_line,exon)
            with open("data/exon_low_coverage.txt","w") as f:
                f.write("{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}\n".format("Sample ID","\t","Order","\t","Sequence","\t","Reads","\t","Percentage in total reads","\t","Exon","\t","Percentage in each exon","\t","Genotype"))
                for l in ex_seq:
                    f.write("{}{}{}{}{}{}{:d}{}{:.2%}{}{}{}{:.2%}{}{}\n".format(inputfilename[2:],"\t"," ", "\t",l[0],"\t",l[1],"\t",l[2],"\t",l[3],"\t",l[4],"\t",l[5]))  
        
    with open("data/genetyping_completed.txt","w") as f: 
        f.write("genotyping complete\n")
  

# abstract read1 and read2 and put into seq_read12 list.
# create a set with only unique reads.
# count each unique read. 
# calculate percentage of each unique read in total reads             
def get_unique(data_in):
    seq_read12 =[] 
    for i in range(len(data_in)):  
        barcoderead1read2=data_in[i].split('\t')
        seq_read12.append(barcoderead1read2[1])
        seq_read12.append(barcoderead1read2[2])
    total_read =len(seq_read12)
    unique_seq = set(seq_read12)                        ### create a set with only unique read.              
    counted_line = []
    for seq in unique_seq: 
        if len(seq) >100:                                   ### remove primer dimers
            count_read = seq_read12.count(seq)              ### count each unique read.
            percent = float(count_read)/float(total_read)
            counted_line.append([seq, count_read, percent]) ### counted_line is a list of list.
    return counted_line
    
# mark exon number for each Read
# calculate total reads for each exon
# sort reads by count_read
       
def count_exon_reads(counted_line,exon_reads):       
    for line in counted_line:
        if "ATTCTACCCTTTCCTAT" in line[0] or "TACCATGTTTGAACACTAC" in line[0]:
            line.append("Exon2")
            exon_reads["Exon2"] += float(line[1])
        elif "CATTGACTACATAGCAAGTT" in line[0] or "AAGTCCCTTTCTGCTCT" in line[0]:
            line.append("Exon3")
            exon_reads["Exon3"] += float(line[1])              
        elif "ACAGTCCTGAGGGTTCTTC" in line[0] or "TGAAAGGTCACACGAAAGACA" in line[0]:
            line.append("Exon4")
            exon_reads["Exon4"] += float(line[1])
        elif "TTAGTTCCCCACCTTGACAC" in line[0] or "CACAGGCGGTGTTGAG" in line[0]:
            line.append("Exon5")
            exon_reads["Exon5"] += float(line[1])
        elif "CTTCTGACTTATTCTTGCTCT" in line[0] or "GCTAGAAAGTCAACATCAGTC" in line[0]:
            line.append("Exon6")
            exon_reads["Exon6"] += float(line[1])
        elif "GGAATATCCCTACTCTACA" in line[0] or "AACAGGCTAACCTAACCTAC" in line[0]:
            line.append("Exon7")
            exon_reads["Exon7"] += float(line[1])
        elif "CCTAGTTTACACACAGTCAGGA" in line[0] or "TCCGCCTCCTTGGT" in line[0]:
            line.append("Exon8")
            exon_reads["Exon8"] += float(line[1])
        elif  "CACCTCTTGCTCTCTCCT" in line[0] or "AAGCTAATGTCACGGCTAGAG" in line[0]:
            line.append("Exon9")
            exon_reads["Exon9"] += float(line[1])
        elif "CAAAAACCTGTAAGTGGAGC" in line[0] or "AGGTCTGGGTAGAGCACC" in line[0]:
            line.append("Exon10")
            exon_reads["Exon10"] += float(line[1])
        elif "AGCCCAAACTGCTAGCTC" in line[0] or "CATAAGACAGCAAGGAGAG" in line[0]:
            line.append("Exon11")
            exon_reads["Exon11"] += float(line[1])
        else:
            line.append("unknown")
            exon_reads["unknown"] += float(line[1])
    print("unique reads:", len(counted_line))
    print(exon_reads)    
    return counted_line, exon_reads

# calculate percentage of each read in total reads of each exon 
def filter_reads(counted_line, exon_reads):
    for line in counted_line:
        if line[3] == "Exon2":
            percent_in_exon = float(line[1])/exon_reads["Exon2"]
            line.append(percent_in_exon)
            if exon2_F in line[0] or exon2_R in line[0]:
                line.append("WT")
            else:
                line.append("Check")
        elif line[3] == "Exon3":
            percent_in_exon = float(line[1])/exon_reads["Exon3"]
            line.append(percent_in_exon)
            if  exon3_F in line[0] or exon3_R in line[0]:
                line.append("WT")
            else:
                line.append("Check")
        elif line[3] == "Exon4":
            percent_in_exon = float(line[1])/exon_reads["Exon4"]
            line.append(percent_in_exon)
            if exon4_F in line[0] or exon4_R in line[0]:
                line.append("WT")
            else:
                line.append("Check")            
        elif line[3] == "Exon5":
            percent_in_exon = float(line[1])/exon_reads["Exon5"]
            line.append(percent_in_exon)
            if exon5_F in line[0] or exon5_R in line[0]:
                line.append("WT")
            else:
                line.append("Check")                
        elif line[3] == "Exon6":
            percent_in_exon = float(line[1])/exon_reads["Exon6"]
            line.append(percent_in_exon)
            if exon6_F in line[0] or exon6_R in line[0]:
                line.append("WT")
            else:
                line.append("Check")
        elif line[3] == "Exon7":
            percent_in_exon = float(line[1])/exon_reads["Exon7"]
            line.append(percent_in_exon)
            if exon7_F in line[0] or exon7_R in line[0]:
                line.append("WT")
            else:
                line.append("Check")
        elif line[3] == "Exon8":
            percent_in_exon = float(line[1])/exon_reads["Exon8"]
            line.append(percent_in_exon)
            if exon8_F in line[0] or exon8_R in line[0]:
                line.append("WT")
            else:
                line.append("Check")
        elif line[3] == "Exon9":
            percent_in_exon = float(line[1])/exon_reads["Exon9"]
            line.append(percent_in_exon)
            if exon9_F in line[0] or exon9_R in line[0]:
                line.append("WT")
            else:
                line.append("Check")
        elif line[3] == "Exon10":
            percent_in_exon = float(line[1])/exon_reads["Exon10"]
            line.append(percent_in_exon)
            if exon10_F in line[0] or exon10_R in line[0]:
                line.append("WT")
            else:
                line.append("Check")
        elif line[3] == "Exon11":
            percent_in_exon = float(line[1])/exon_reads["Exon11"]
            line.append(percent_in_exon)
            if exon11_F in line[0] or exon11_R in line[0]:
                line.append("WT")
            else:
                line.append("Check")
        elif line[3] == "unknown":
            percent_in_exon = float(line[1])/exon_reads["unknown"]
            line.append(percent_in_exon)
            line.append("Check")
        else:
            print("exon is unknown")

    ### sort unique sequence by count_read
    counted_line.sort(key=itemgetter(1), reverse=True)   # cannot assign counted_line.sort(key=itemgetter(1), reverse=True) to a new variable    
    return counted_line


## analyze exon with low coverage
def calculate_coverage(counted_line,exon):
    ex_list=exon.split(',')
    print(ex_list)
    ex_seq=[]
    for ex in ex_list:    
        reads=0
        for line in counted_line:             
            if line[3] == ex:
                reads += line[1]
                ex_seq.append(line)
                if reads >1000:
                    break
    print(len(ex_seq))
    return ex_seq
    
main()    
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	