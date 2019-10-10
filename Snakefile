from snakemake.remote.FTP import RemoteProvider

report: "report/workflow.rst"

configfile: "config.yaml"


populations = config["populations"] if "populations" in config.keys() else ["CEU","YRI","LWK"]
MAC = config["MAC"] if "MAC" in config.keys() else 0
chromosomes = config["chromosomes"] if "chromosomes" in config.keys() else [21, 22]

remove_linkage = config["remove_linkage"] if "remove_linkage" in config.keys() else True
use_plink_to_prune = config["use_plink_to_prune"] if "use_plink_to_prune" in config.keys() else True
if use_plink_to_prune:
	linkage_command = config["linkage_command"] if "linkage_command" in config.keys() else " --indep-pairwise 50 10 0.1"
else:
	reduction_factor = config["reduction_factor"] if "reduction_factor" in config.keys() else 10


all_pairs = []

for i in range(len(populations)):
	for j in range(i+1,len(populations)):
		all_pairs.append((populations[i],populations[j]))

ruleorder: run_medeas_single_pop > run_medeas


rule all:
    input:
        expand("medeas/single/{pop1}/{pop1}.log",pop1 = populations),
        expand("medeas/double/{pops[0]}_{pops[1]}/{pops[0]}_{pops[1]}.log", pops = all_pairs),
        "medeas/all/pop_all/pop_all.log"


FTP = RemoteProvider()

rule initial_bcf:
    input:
        vcf=FTP.remote("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{chromosome}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"),
        exclude="sampleToExclude.lab"
    output:
        "bcf/initial_bcf.chr{chromosome}.bcf"
    shell:
        "bcftools view -S ^{input.exclude} -O b -i \"MAC > {MAC} \" {input.vcf} -o {output}"




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
	    lab = "label/single/{pop}.lab",
	    lab_haplo = "label/single/{pop}_haplo.lab",
	    ids = "label/single/{pop}.ids",
	    ids_haplo = "label/single/{pop}_haplo.ids"
    run:
        shell("grep  {wildcards.pop} {input} | awk '{{print $1}}' > {output.ids}"),
        shell("grep  {wildcards.pop} {input} | awk -v population={wildcards.pop} '{{print population}}' > {output.lab}"),
        shell("grep  {wildcards.pop} {input} | awk '{{print $1; print $1}}' > {output.ids_haplo}"),
        shell("grep  {wildcards.pop} {input} | awk -v population={wildcards.pop} '{{print population; print population}}' > {output.lab_haplo}")

rule build_panel_all_pop:
    input:
	    "label/sample.prunned.csv"
    output:
        ids = "label/all/pop_all.ids",
        lab = "label/all/pop_all.lab",
        ids_haplo = "label/all/pop_all_haplo.ids",
        lab_haplo = "label/all/pop_all_haplo.lab"
    run:
        shell("rm -f label/pop_all.lab")
        shell("rm -f label/pop_all_haplo.ids")
        for population in populations:
            shell("grep  {population} {input} | awk '{{print $1}}' >> {output.ids}"),
            shell("grep  {population} {input} | awk -v population={population} '{{print population}}' >> {output.lab}"),
            shell("grep  {population} {input} | awk '{{print $1; print $1}}' >> {output.ids_haplo}"),
            shell("grep  {population} {input} | awk -v population={population} '{{print population; print population}}' >> {output.lab_haplo}")

rule build_panel_two_pops:
    input:
	    p1 = "label/single/{pop1}.lab",
	    p2 = "label/single/{pop2}.lab",
	    p3 = "label/single/{pop1}_haplo.lab",
	    p4 = "label/single/{pop2}_haplo.lab"
    output:
        "label/double/{pop1}_{pop2}.lab",
        "label/double/{pop1}_{pop2}_haplo.lab"
    run:
        shell("cat {input.p1} {input.p2} > label/double/{wildcards.pop1}_{wildcards.pop2}.lab"),
        shell("cat {input.p3} {input.p4} > label/double/{wildcards.pop1}_{wildcards.pop2}_haplo.lab")

rule build_panel_three_pops:
    input:
	    p1 = "label/single/{pop1}.lab",
	    p2 = "label/single/{pop2}.lab",
	    p3 = "label/single/{pop3}.lab",
	    p4 = "label/single/{pop1}_haplo.lab",
	    p5 = "label/single/{pop2}_haplo.lab",
	    p6 = "label/single/{pop3}_haplo.lab"
    output:
        lab = "label/triple/{pop1}_{pop2}_{pop3}.lab",
        lab_haplo = "label/triple/{pop1}_{pop2}_{pop3}_haplo.lab"
    run:
        shell("cat {input.p1} {input.p2} {input.p3} > {output.lab}"),
        shell("cat {input.p4} {input.p5} {input.p6} > {output.lab_haplo}")


rule remove_population:
    input:
        bcf="bcf/initial_bcf.chr{chromosome}.bcf",
        pop_id="label/all/pop_all.ids"
    output:
         "bcf/reduce.chr{chromosome}.bcf"
    shell:
        "bcftools view -S {input.pop_id} -m 1 -M 2 -v snps -O b -o {output} {input.bcf}"


rule get_unlinked_site_plink:
    input:
        "bcf/reduce.chr{chromosome}.bcf"
    output:
        "plink/linked.chr{chromosome}.prune.in",
    params:
        "plink/linked.chr{chromosome}"
    shell:
        "plink --bcf {input} {linkage_command} --out {params} --allow-extra-chr"

rule get_unlinked_site_random:
    input:
        "bcf/reduce.chr{chromosome}.bcf"
    output:
        "random_site/unlinked.chr{chromosome}.in",
    params:
        "plink/linked.chr{chromosome}"
    shell:
        "nline=\"$(bcftools view   {input} |grep -v \"^##\" | wc -l)\";nbsnp=\"$(echo \"scale=0;$nline/{reduction_factor}\" | bc)\";bcftools view {input} | grep -v \"^##\" | cut -f3 | shuf -n $nbsnp\ > {output}"

rule remove_linked_site:
    input:
        bcf = "bcf/reduce.chr{chromosome}.bcf",
        site_list = "plink/linked.chr{chromosome}.prune.in" if use_plink_to_prune == True else "random_site/unlinked.chr{chromosome}.in"
    output:
        "bcf/prunned.chr{chromosome}.bcf"
    shell:
        "bcftools view {input.bcf} -i ID=@{input.site_list} -O b -o {output}"

rule merge_bcf:
    input:
    	expand("bcf/prunned.chr{chromosome}.bcf", chromosome = chromosomes) if remove_linkage else expand("bcf/reduce.chr{chromosome}.bcf", chromosome = chromosomes)
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
        pop = "label/all/pop_all_haplo.lab"
    output:
        "med/single/{population}_sorted.med"
    shell:
        "pop=\"$(nl  {input.pop} | grep {wildcards.population} | awk '{{printf(\"%s,\",$1);}}')\" && "
        "pop=${{pop::-1}} && "
        "cut -d\" \" -f${{pop}} {input.med} > {output}"


rule make_med_single_pop:
    input:
        "med/single/{population}_sorted.med"
    output:
        "med/single/{population}_full.med"
    shell:
	    "awk 'BEGIN{{srand()}}{{for (i=0;i<NF/2;i++) if (rand() > 0.5) printf \"%s %s \", $(2*i+1) ,$(2*i+2); else printf \"%s %s \", $(2*i+2), $(2*i+1);print \"\"}}' "
	    "{input} > "
	    "{output}"


rule make_med_two_pop_prunned:
    input:
        pop1 = "med/single/{population1}_full.med" ,
        pop2 = "med/single/{population2}_full.med"
    output:
        "med/double/{population1}_{population2}_prunned.med"
    shell:
        "paste -d \"\" {input.pop1} {input.pop2} |"
        "awk '{{ for(i=1; i<=NF;i++) j+=$i; if(j < NF*2) print $0; j=0 }}'  > {output} "

rule make_med_three_pop_prunned:
    input:
        pop1 = "med/single/{population1}_full.med" ,
        pop2 = "med/single/{population2}_full.med" ,
        pop3 = "med/single/{population3}_full.med"
    output:
        "med/triple/{population1}_{population2}_{population3}_prunned.med"
    shell:
        "paste -d \"\" {input.pop1} {input.pop2} {input.pop3}  |"
        "awk '{{ for(i=1; i<=NF;i++) j+=$i; if(j < NF*2) print $0; j=0 }}'  > {output} "

rule make_med_single_pop_prunned:
    input:
       "med/single/{population}_full.med"
    output:
        "med/single/{population}_prunned.med"
    shell:
	    "awk '{{ for(i=1; i<=NF;i++) j+=$i; if(j < NF*2) print $0; j=0 }}' {input} > {output}"

rule make_all_pop_prunned:
    input:
       expand("med/single/{population}_full.med",population = populations)
    output:
        "med/all/pop_all_prunned.med"
    shell:
        "paste -d \"\" {input} | awk '{{ for(i=1; i<=NF;i++) j+=$i; if(j < NF*2) print $0; j=0 }}' > "
        "{output}"

rule run_medeas:
    input:
        med_file = "med/{folder}/{pop_pattern}_prunned.med",
        lab_file = "label/{folder}/{pop_pattern}_haplo.lab"
    output:
        log = "medeas/{folder}/{pop_pattern}/{pop_pattern}.log",
    shell:
        "python ~/medeas/main.py "
        "-sf {input.med_file} "
        "-lf {input.lab_file} "
        "-of medeas/{wildcards.folder}/{wildcards.pop_pattern} "
        "-bws 100000 -bsn 100 "
        "--ncpus 10 "
        "> {output.log}"

rule run_medeas_single_pop:
    input:
        med_file = "med/single/{pop_pattern}_prunned.med",
        lab_file = "label/single/{pop_pattern}_haplo.lab"
    output:
        log = "medeas/single/{pop_pattern}/{pop_pattern}.log",
    shell:
        "python ~/medeas/main.py "
        "-sf {input.med_file} "
        "-lf {input.lab_file} "
        "-of medeas/single/{wildcards.pop_pattern} "
        "-bws 100000 -bsn 100 "
        "--ncpus 10 "
        "> {output.log} "
        "|| true"
