# mds_1000_genome

This repository contain a snakemake pipeline to download the 1000 genome project, filter the data and run medeas on them. 

It needs a configuration name `config.yam` as well as a list of the sample to be excluded name `sampleToExclude.lab`. The configuration file can be empty, or the following option can be specified:

 `populations`: A list of populations to use (default ["CEU","YRI","LWK"])
 `MAC`: A number specifying a global MAC (Max Allele Count) cutoff (default 0)
`chromosomes`: A list of chromosome to use (default [21 22])
`remove_linkage`: A boolean specifying if linkage should be removed (default True) 
`use_plink_to_prune` A boolean specifying if to use plink to prune the data for linkage (True) or to use random removal of site (False) (default True)
`linkage_command` The plink command to remove linkage (default `--indep-pairwise 50 10 0.1`)
`reduction_factor` One over the proportion of site to randomly remove (default 10, i.e. we randomly select 1/10 sites)

Some configuration files are provided as examples into to config_file directory. Names as self-explanatory. 

If the configuration file and the sampleToExclude.lab file are place in the same directory as the snakemake file, it can be directly run using:
`snakemake --cores 4`


Alternatively, if the configuration file is in a different directory as the snakemake file, it can be run using: 
`snakemake --directory path/to/configurationfile --cores 4`

In this case, to generate all the data without any pruning, simply run: 
`snakemake --directory config_file/no_prunning --cores 4`.

A representation of the pipeline can be seen using (the directory need to be specified if needed):

`snakemake --dag | dot -Tpdf > dag.pdf`

A conda environment in which the pipeline works can be created using:

`conda env create -n mds_1000_genome -f conda_envs.yaml`

`source activate mds_1000_genome`

It can be run on a slurm cluster using the cluster.json on a slurm cluster using:

`Snakemake -j 80 --directory path/to/configurationfile --cluster-config cluster.json --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -A {cluster.account} --cpus-per-task={cluster.cpus} --job-name={cluster.name} --error={cluster.error} --output={cluster.output}"`
