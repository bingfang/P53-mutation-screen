# snakemake --dry-run check_target
# snakemake -np check_target
# snakemake check_target --cores 1
# snakemake -np

rule check_target:
    input: 
        "data/ex5end_snp_check.txt",
        "data/P53_S1_L001_R1_trimmed.fastq",
        "data/genetyping_completed.txt",
        "data/exon7_2R_snp_check.txt"   

        
rule trim:
    input: 
        R1="data/P53_S1_L001_R1_001.fastq.gz",
        R2="data/P53_S1_L001_R2_001.fastq.gz",
        myPath= "/Users/xubr/anaconda3/share/trimmomatic-0.39-2"
    output:
        out1="data/P53_S1_L001_R1_trimmed.fastq",
        out2="data/P53_S1_L001_R2_trimmed.fastq"
    params:
    shell:
        """
        java -jar {input.myPath}/trimmomatic.jar PE -phred33 {input.R1} {input.R2} {output.out1} P53_S1_L001_R1_unpair_trimmed.fastq {output.out2} P53_S1_L001_R2_unpair_trimmed.fastq TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        """



rule debarcode:
    input: 
        script="src/_0debarcode_Miseq.py",
        sample_file="data/sample.txt"  
    output: 
    params:
        "data/sample.txt"
    shell:
        """
        python3 {input.script} {input.sample_file}
        """


rule genotyping:
    input: 
        script="src/_1_P53_genotyped_mouseTrp53.py"   
    output:
        "data/genetyping_completed.txt", 
        "data/Check_seq_all.txt",
        "data/P53_Miseq_screen_genotyping.txt"
    params:
    shell:
        """
        python3 {input.script}
        """


rule calculate_common_mutation_rate:
    input: 
        script2="src/_3CheckSNP_mouseTrp53.py"   
    output: 
        "data/ex5end_snp_check.txt"
    params:
    shell:
        """
        python3 {input.script2}
        """
        
rule calculate_mutation_rate:
    input: 
        script1="src/_2CheckSNP_mouseTrp53.py" 
    output: 
        "data/exon7_2R_snp_check.txt"
    params:
    shell:
        """
        python3 {input.script1}
        """        
        
     
         
            