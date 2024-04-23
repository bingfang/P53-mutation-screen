#!/usr/local/bin/python3.6

def main():
    with open("data/genetyping_coverage.txt", "r") as f:
        data_in = f.read().rstrip().split("\n")
        print(len(data_in))
        
    new,sample_num = re_format(data_in)

    with open("data/QC_coverage.txt", "w") as f:
        for l in new:
            #print(l)
            f.write(l+"\n")
        
def re_format(data_in):
    new=[]
    sample_num = 0
    for line in data_in:
        if "data" not in line:
            count=line.split("\t")
            if len(count)<2:
                print(line)
            else:
                count_num= float(count[1])
  
                if count_num < 200:
                   line = line + "\t low coverage"
                   
        else:
            line=line 
            sample_num +=1        


        new.append(line) 
    print(len(new),sample_num, type(new))   

    return new, sample_num
    
main()