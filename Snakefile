# snakemake --dry-run check_target
# snakemake -np check_target
# snakemake check_target --cores 1
# snakemake -np
# snakemake --unlock

# check input file name
# built sample.txt file 

rule check_target:
    input: 
        "data/P53_Miseq_screen_genotyping.txt",
        "data/ex5end_snp_check.txt",
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
        R1="data/P53_S1_L001_R1_trimmed.fastq",
        R2="data/P53_S1_L001_R2_trimmed.fastq",
        script="src/_0debarcode_Miseq.py",
        file="data/sample.txt"
    output: 
        "data/debarcode_completed.txt"
    params:
    shell:
        """
        python3 {input.script}
        """


rule genotyping:
    input:
        debarcode="data/debarcode_completed.txt",
        script="src/_1_P53_genotyped_mouseTrp53_Coverage.py"   
    output:
        "data/genetyping_coverage.txt", 
        "data/Check_seq_all.txt",
        "data/P53_Miseq_screen_genotyping.txt",
        "data/genetyping_completed.txt"
    params:
    shell:
        """
        python3 {input.script}
        """


rule calculate_common_mutation_rate:
    input:
        genotype="data/genetyping_completed.txt",
        mutation="data/Check_seq_all.txt",
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
        genotype="data/genetyping_completed.txt", 
        mutation="data/Check_seq_all.txt",
        script1="src/_2CheckSNP_mouseTrp53.py" 
    output: 
        "data/exon7_2R_snp_check.txt"
    params:
    shell:
        """
        python3 {input.script1}
        """        
        
