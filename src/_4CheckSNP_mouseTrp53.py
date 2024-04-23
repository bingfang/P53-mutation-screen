#!/usr/local/bin/python3.8

#conda install -c bioconda biopython
#conda activate mP53
from Bio.Seq import Seq
import glob

## mutation filter is set as 10%
## R1 and R2 imbalance filter is set as 5%


def main():
    #list=["wt seq","mutated seq"]
    list1 = ['ATGAGGCCTTAGAG', 'ATGAGCCCTTAGA']    #Exon10n
    list2 = ['AGAGTTAAAGGAT','AGAGTAAAAGGAT'] # Exon10n1
    list3 = ['CCTGCACAAGCG','CCTGCCCAAGCG'] # exon9n
    list4 = ['CTCAAGGTACCAA','CTCAAGTTACCAA'] # exon9nend
    list5 = ['AAGAAAATTTCCG','AAGAACATTTCCG'] # exon8n
    list6 = ['AAGAAAATTTCCG','AAGAAAAGTTCCG'] # exon8n2
    list7 = ['CCCGAGGTCTGTAA','CCCGAGATCTGTAA'] # exon6_end
    list8 = ['CCCACCATGAGC','CCCCCCCACGAG'] # exon5n
    list9 = ['ACCGAGACCCCTG','ACCGAGCCCCCTG'] # exon4n
    list10 = ['TTTCCAGACTTCC','TTTCCAAACTTCC'] # exon3_start 
    list11 = ['GAGCTCCCTCT','GAGCCCCCTCT'] # exon2n    
    list12 = ['CGATGGTGATGGT','CGATGGTGATGGT'] # exon5end
    list13 = ['GAGTATACCACC','GAGTATCCCACC'] # exon7n_T224P
    #list14 = ['ACATGTGTAATA','ACATGAGTAAT'] # exon7-mid
    #list15 = ['TAATAGCTCCTGC','TAATAGGTCCTGC'] # exon7-3

    dic={"ex10n":list1, "ex10n1":list2, "ex9n":list3,"ex9nend":list4,"ex8n":list5, "ex8n2":list6,"ex6_end":list7,"ex5n":list8, "ex4n":list9,"ex3_start":list10,"ex2n":list11,"ex5end":list12,"ex7n_T224P":list13}
   
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

