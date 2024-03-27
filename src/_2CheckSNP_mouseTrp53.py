#!/usr/local/bin/python3.6
from Bio.Seq import Seq
import glob

## mutation filter is set as 10%
## R1 and R2 imbalance filter is set as 5%


def main():
    #list=["wt seq","mutated seq"]
    list1 = ['GCCCTCCGAGTGTCA', 'GCCCTCTGAGTGTCA']    #Exon4   
    list2 = ['GGGAGAGACCGCCG','GGGAGAAACCGCCG'] # exon8 D275N
    list3 = ['AGAAAATTTCCGCA','AGAAAAGTTCCGCA'] # exon8-1
    list4 = ['GTGGGAACCTTC','GTGGGACCCTTC'] # exon8-2
    list5 = ['TATACCACCATCC','TATACCCCCATCC'] # exon7_1
    list6 = ['ACCATCCACTACA','ACCAGCCACTACA'] # exon7_2 reoccur
    list7 = ['TTTTCGCCACAG','TTTTCCCCACAG'] # exon6 reoccur
    list8 = ['AACCGCCGACC','AACCCCCGACC'] # exon7_3n
    list9 = ['TCCACAGCCATC','TCCACACCCATC'] # exon4_s
    
    dic={"ex4_1":list1, "ex8_D275N":list2, "ex8_1":list3,"ex8_2":list4, "ex7_1":list5, "exon7_2R":list6,"exon6R":list7,"exon7_3n":list8,"exon4s":list9}
   
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

