from snakemake.remote.FTP import RemoteProvider
chromosomes=range(1,23)
populations=["CEU","YRI","LWK","TSI","JPT","CHB"]
shell.prefix("module add UHTS/Analysis/plink/1.90; ")

all_pairs = []
for i in range(len(populations)):
	for j in range(i+1,len(populations)):
		all_pairs.append((populations[i],populations[j]))


singularity: "docker://continuumio/miniconda3"
ruleorder: build_panel_all_pop > build_panel_single_pop
ruleorder: make_all_pop_prunned > make_med_single_pop_prunned 




rule all:
    input:
        expand("medeas/pop_{pop1}/pop_{pop1}.log",pop1 = populations),
        expand("medeas/pops_{pops[0]}_{pops[1]}/pops_{pops[0]}_{pops[1]}.log", pops = all_pairs),
        "medeas/pop_all/pop_all.log"


FTP = RemoteProvider()

rule initial_bcf:
    input:
        vcf=FTP.remote("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{chromosome}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"),
        exclude="sampleToExclude.lab"
    output:
        "bcf/initial_bcf.chr{chromosome}.bcf"
    shell:
        "bcftools view -S ^{input.exclude} -O b -i \"MAC > 2\" {input.vcf} -o {output}"




rule build_general_panel:
    input:
        vcf=FTP.remote("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"),
        exclude="sampleToExclude.lab"        
    output:
	    "label/sample.prunned.csv"
    run:
        shell("mkdir -p label")
        shell("grep -v -f {input.exclude} {input.vcf} > {output}")


rule build_panel_single_pop:
    input:
	    "label/sample.prunned.csv"
    output:
	    lab = "label/pop_{pop1}.lab",
	    lab_haplo = "label/pop_{pop1}_haplo.lab"
    run:
        for population in populations:
            shell("grep  {population} {input} | awk '{{print $1}}' > label/pop_{population}.ids"),
            shell("grep  {population} {input} | awk -v population={population} '{{print population}}' > label/pop_{population}.lab"),
            shell("grep  {population} {input} | awk '{{print $1; print $1}}' > label/pop_{{population}}_haplo.ids"),
            shell("grep  {population} {input} | awk -v population={population} '{{print population; print population}}' > label/pop_{population}_haplo.lab")

rule build_panel_all_pop:
    input:
	    "label/sample.prunned.csv"
    output:
        "label/pop_all.ids",        
        "label/pop_all_haplo.lab"
    run:
        shell("rm -f label/pop_all.lab")
        shell("rm -f label/pop_all_haplo.ids")
        for population in populations:
            shell("grep  {population} {input} | awk '{{print $1}}' >> label/pop_all.ids"),
            shell("grep  {population} {input} | awk -v population={population} '{{print population}}' >> label/pop_all.lab"),
            shell("grep  {population} {input} | awk '{{print $1; print $1}}' >> label/pop_all_haplo.ids"),
            shell("grep  {population} {input} | awk -v population={population} '{{print population; print population}}' >> label/pop_all_haplo.lab")
            
rule build_panel_two_pops:
    input:
	    p1 = "label/pop_{pop1}.lab",
	    p2 = "label/pop_{pop2}.lab",
	    p3 = "label/pop_{pop1}_haplo.lab",
	    p4 = "label/pop_{pop2}_haplo.lab"
    output:     
        "label/pops_{pop1}_{pop2}.lab",
        "label/pops_{pop1}_{pop2}_haplo.lab"
    run:
        for population in populations:
            shell("cat {input.p1} {input.p2} > label/pops_{wildcards.pop1}_{wildcards.pop2}.lab"),
            shell("cat {input.p3} {input.p4} > label/pops_{wildcards.pop1}_{wildcards.pop2}_haplo.lab")
            
            
rule remove_population:
    input:
        bcf="bcf/initial_bcf.chr{chromosome}.bcf",
        pop_id="label/pop_all.ids"
    output:
         "bcf/reduce.chr{chromosome}.bcf"
    shell:
        "bcftools view -S {input.pop_id} -m 1 -M 2 -v snps -O b -o {output} {input.bcf}"	
    


rule get_linked_site:
    input:
        "bcf/reduce.chr{chromosome}.bcf"
    output:
        "plink/linked.chr{chromosome}.prune.in",
    params:
        "plink/linked.chr{chromosome}"
    shell:
        "plink --bcf {input} --indep 100 10 1.01 --out {params} --allow-extra-chr"

rule remove_linked_site:
    input:
        bcf = "bcf/reduce.chr{chromosome}.bcf",
        site_list = "plink/linked.chr{chromosome}.prune.in"
    output:
        "bcf/prunned.chr{chromosome}.bcf"
    shell:
        "bcftools view {input.bcf} -i ID=@{input.site_list} -O b -o {output}"

rule merge_bcf:
    input:
    	expand("bcf/prunned.chr{chromosome}.bcf", chromosome = chromosomes) if False else expand("bcf/initial_bcf.chr{chromosome}.bcf", chromosome = chromosomes)
    output:
        "bcf/all_merged.bcf"
    shell:
        "bcftools concat {input} -O b -o {output}"


rule make_tbed:
    input:
        "bcf/all_merged.bcf"
    output:
        "plink/all_population_sorted.tped"
    shell:
        "plink --bcf  {input}  --allow-extra-chr --recode12 --transpose --out plink/all_population_sorted"

rule make_med_all_pop_sorted:
    input:
        "plink/all_population_sorted.tped"
    output:
        "med/all_pop_sorted1.med"    
    shell:
        "cut -d\" \" -f 1-4 --complement {input} > {output}"
        
rule make_med_single_pop_sorted:
    input:
        med = "med/all_pop_sorted1.med",
        pop = "label/pop_all_haplo.lab"
    output:
        "med/pop_{population}_sorted.med"    
    shell:
        "pop=\"$(nl  {input.pop} | grep {wildcards.population} | awk '{{printf(\"%s,\",$1);}}')\" && "
        "pop=${{pop::-1}} && "
        "cut -d\" \" -f${{pop}} {input.med} > {output}"
        
        
rule make_med_single_pop:
    input:
        "med/pop_{population}_sorted.med" 
    output:
        "med/pop_{population}_full.med"    
    shell:
	    "awk 'BEGIN{{srand()}}{{for (i=0;i<NF/2;i++) if (rand() > 0.5) printf \"%s %s \", $(2*i+1) ,$(2*i+2); else printf \"%s %s \", $(2*i+2), $(2*i+1);print \"\"}}' "
	    "{input} > "
	    "{output}"


rule make_med_two_pop_prunned:
    input:
        pop1 = "med/pop_{population1}_full.med" ,
        pop2 = "med/pop_{population2}_full.med" 
    output:
        "med/pops_{population1}_{population2}_prunned.med"    
    shell:
        "paste -d \"\" {input.pop1} {input.pop2} |"
        "awk '{{ for(i=1; i<=NF;i++) j+=$i; if(j < NF*2) print $0; j=0 }}'  > {output} "
        
rule make_med_single_pop_prunned:
    input:
       "med/pop_{population}_full.med"
    output:
        "med/pop_{population}_prunned.med"    
    shell:
	    "awk '{{ for(i=1; i<=NF;i++) j+=$i; if(j < NF*2) print $0; j=0 }}' {input} > {output}"
	    
rule make_all_pop_prunned:
    input:
       expand("med/pop_{population}_full.med",population = populations)
    output:
        "med/pop_all_prunned.med"    
    shell:
        "paste -d \"\" {input} | awk '{{ for(i=1; i<=NF;i++) j+=$i; if(j < NF*2) print $0; j=0 }}' > "
        "{output}"


rule run_medeas:
    input:
        med_file = "med/{pop_pattern}_prunned.med",
        lab_file = "label/{pop_pattern}_haplo.lab"        
    output:
        log = "medeas/{pop_pattern}/{pop_pattern}.log",
    shell:
        "singularity exec --bind /scratch/local/monthly/fmichau2/ ~/python3.simg python3 ~/medeas/main.py "
        "-sf {input.med_file} " 
        "-lf {input.lab_file} "
        "-of medeas/{wildcards.pop_pattern} "
        "-bws 100000 -bsn 100 "
        "> {output.log}"
        