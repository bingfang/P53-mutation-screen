#!/usr/local/bin/python3.8

#conda install -c bioconda biopython
#conda activate mP53
from Bio.Seq import Seq
import glob

## mutation filter is set as 10%
## R1 and R2 imbalance filter is set as 5%


def main():
    #list=["wt seq","mutated seq"]
    list1 = ['GGGAGCCGTGTCC', 'GGGAGGCGTGTCC']    #Exon5 S149R
    list2 = ['TGAGCGCTGCTC','ATGGTGATGGTA'] # exon5end
    list3 = ['GCCATGGCCATCT','GCCATGCCCATCT'] # exon5_1 A155P
    list4 = ['ACGTGCCCTGTG','ACGTGGCCTGTG'] # exon5_2 C135W
    list5 = ['CCAGTACTCTCC','CCAGTGCTCTCC'] # exon5_3 reoccur Y120C
    list6 = ['GAGGTCGTGAGAC','GAGGTCCTGAGAC'] # exon5_4n V167L
    list7 = ['ATTCTGCCAGCTG','ATTCTGGCAGCTG'] # exon5_5n C129W
    list8 = ['ATATCCTGGTA','ATATACTGGTA'] # exon3n
    list9 = ['TGAGGCCTTAG','TGAGCCCTTAG'] # exon10 A341
    list10 = ['AGATATCCTGGT','AGATATACTGGT'] # exon3 
    list11 = ['GGGCCAGTGGCCCC','GGGCCAATGGCCCC'] # exon4, GGCCCC 
    list12 = ['CATGAGCGCTG','CATGAACGCTG'] # exon5new
    list13 = ['GACCTATCCTTAC','GACCTTTCCTTAC'] # exon7n
    list14 = ['ACATGTGTAATA','ACATGAGTAAT'] # exon7-mid
    list15 = ['TAATAGCTCCTGC','TAATAGGTCCTGC'] # exon7-3

    dic={"ex5_S149R":list1, "ex5end":list2, "ex5_A155P":list3,"ex5_C129W":list4,"ex5_Y120C":list5, "ex5_4V167L":list6,"ex5_5n":list7,"ex3n":list8, "ex10_A341P":list9,"ex3":list10,"ex4-GGCCC":list11, "ex5new":list12,"ex7n":list13,"ex7-mid":list14, "ex7-3":list15}
   
    for i in dic:
        print(i, dic[i])
        outputfilename = 'data/' + str(i) + "_snp_check.txt" 
        for name in glob.glob('data/1*2*.txt'):
            inputfilename = str(name)
            with open(inputfilename, 'r') as f:
                data_in = f.read().rstrip().split('\n')
            WT = Seq(dic[i][0])
            WT_RC = WT.reverse_complement()
            mut= Seq(dic[i][1])
            mut_RC = mut.reverse_complement()
             
            wt_F, wt_R, mut_F, mut_R,F_rate,R_rate = check_SNP(data_in,str(WT),str(WT_RC),str(mut),str(mut_RC))
                       
            note = make_note(name,F_rate, R_rate)

            with open(outputfilename, "a") as f:
                f.write(inputfilename[2:]+'\n')
                f.write("{}{}{}{}\n".format("WT_Forward\t",str(WT),"\t",wt_F))
                f.write("{}{}{}{}\n".format("WT_Reward\t",str(WT_RC),"\t",wt_R))
                f.write("{}{}{}{}{}{:.2%}\n".format("mut_Forward\t",str(mut),"\t",mut_F,"\t",F_rate))
                f.write("{}{}{}{}{}{:.2%}{}{}\n".format("mut_Reward\t",str(mut_RC),"\t",mut_R,"\t",R_rate,"\t",note))
        

# calculate forward and reward mutation rate    
def check_SNP(data_in,str1,str2,str3,str4):
    wt_F = 0
    mut_F = 0
    wt_R = 0
    mut_R = 0
    for line in data_in:
        if str1 in line:
            wt_F += 1
        if str2 in line:
            wt_R += 1
        if str3 in line:
            mut_F += 1
        if str4 in line:
            mut_R += 1
    F_rate= float(mut_F)/(float(wt_F)+float(mut_F)+1)   # for exon5 end, no mutation count, add 1 to avoid error
    R_rate= float(mut_R)/(float(wt_R)+float(mut_R)+1)     
    return wt_F, wt_R, mut_F, mut_R, F_rate,R_rate
    
# annotate mutation    
def make_note(name,F_rate, R_rate):
    note=""
    if abs(F_rate - R_rate) > 0.05:
        if F_rate > 0.05 or R_rate > 0.05 :         
            note="R1R2 imbalance" 
    else:
        if F_rate > 0.05 or R_rate > 0.05:
            print("{}{}{:0.2%}{}{:.2%}{}\n".format(name[21:-4],"\t",F_rate,"\t",R_rate,"\tmutation"))
            note="MUTATION!!!"
    return note


main()

